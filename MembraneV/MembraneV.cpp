/*
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SC_PlugIn.h"
#include "Membrane_shape.h"
#include "assert.h"

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

// twiddle-ables
#define SHAPE_SZ 1 // diameter
//#define FEEDFORWARD
#define SELF_LOOP        // for control over tension
#define RIMGUIDES // extra self-loops around edge, which invert signal
#define RIMFILTER
#define AUDIO_INPUT
#define TRIGGER_DURATION 1024 /* number of samples worth of white noise to inject */
#define EUCLIDEAN_POS
// some constants from Brook Eaton's roto-drum
// http://www-ccrma.stanford.edu/~be/drum/drum.htm

#define DELTA 6.0f // distance between junctions
#define GAMMA 8.0f // wave speed

// A unit delay
typedef struct {
  float a;
  float b;
  float c;
  int invert; // whether signal should be inverted
} t_delay;

// a scatter junction with up to 6 input and output delays, plus a
// self loop
typedef struct {
  int ins, outs;
  t_delay *in[6];
  t_delay *out[6];
#ifdef SELF_LOOP
  t_delay *self_loop;
#endif
  float last_val;
#ifdef EUCLIDEAN_POS
  float x,y,input_fac;
#endif
} t_junction;


// supercollider stuff starts here...

// InterfaceTable contains pointers to functions in the host (server).
static InterfaceTable *ft;

// declare struct to hold unit generator state
struct Membrane : public Unit
{
  float yj; // junction admittence, calculated from tension parameter
#ifndef AUDIO_INPUT
  int triggered; // flag
  int excite;    // number of samples left in a triggered excitation
#endif
  t_shape *shape;
  t_junction *junctions;
  t_delay *delays;
  float loss;
  int delay_n;   // number of delays in mesh including self loops etc
  float a1;
#ifdef EUCLIDEAN_POS
  float saved_input_pos,saved_input_width;
  int input_pos;
  float Size;
#endif
};

// declare unit generator functions
extern "C"
{
  void Membrane_next_a(Membrane *unit, int inNumSamples);
  void MembraneCircle_Ctor(Membrane* unit);
  void MembraneHexagon_Ctor(Membrane* unit);
  void Membrane_Dtor(Membrane* unit);
};

////////////////////////////////////////////////////////////////////
//calculate input_fac for each junction
int calc_input_fac(Membrane *unit) {
  float radius = unit->Size*0.5f;
  float input_radius = sc_max(1.0f,unit->saved_input_width*radius);
  //Print("input_radius %f\n",input_radius);
  float center_x = unit->saved_input_pos * radius;
  float center_y = 0.0f;
  t_junction *junctions = unit->junctions;
  double mindist = unit->Size;
  int mindist_i = 0;
  //distances to input
  for (int i = 0; i < unit->shape->points_n; ++i) {
      t_junction *junction = &junctions[i];
	  double fx = junction->x - center_x;
	  double fy = junction->y;//-center_y
	  double distance = sqrt((fx * fx) + (fy * fy));
	  if (distance < mindist)  {
		  mindist = distance;
		  mindist_i = i;
	  }
	  float fac = sc_max(0.0f, 1.0f - distance / input_radius);
      fac = 0.5 - cos(fac * M_PI) * 0.5;
      junction->input_fac = fac;
  }
  return mindist_i;
}
// execute one sample cycle over the mesh
float cycle(Membrane *unit, float input, float yj_r,int input_width, int input_pos) {
  t_delay *delays = unit->delays;
  t_junction *junctions = unit->junctions;

  int i;

  int middle = (int) (unit->shape->points_n / 2);
  float result;
  //float input1 = input/((float)input_width+1);
  float input1 = input;
  int input_pos_end = input_pos + input_width;
#ifdef EUCLIDEAN_POS
  input_pos = unit->input_pos;
#endif
  
  // update all the junctions and waveguides
  for (i = 0; i < unit->shape->points_n; ++i) {
    t_junction *junction = &junctions[i];
    int j;
    float total = 0;

    float yc = unit->yj - junction->ins;

    for (j = 0; j < junction->ins; ++j) {
      total += junction->in[j]->b;
    }

#ifdef SELF_LOOP
    total = 2.0f * (total + (yc * junction->self_loop->b)) * yj_r;
#else
    total *= (2.0f / ((float) junction->ins));
#endif
/*
    if (i < middle) {
      total += (input / middle);
    }
*/
    if (i == input_pos) {//(i == 0) {
      result = total;
    }
#ifdef EUCLIDEAN_POS
	total += junction->input_fac * input1/((float) junction->ins);
#else
	if(i>=input_pos && i <= input_pos_end){
		total += input1/((float) junction->ins);
	}
#endif
	
	/*
	if(i <= input_width){
		total += input1/((float) junction->ins);
	}*/
    total *= unit->loss;
	//one pole
	total = total*(1.f - unit->a1) + unit->a1 * junction->last_val;
	junction->last_val = total;

    for (j = 0; j < junction->outs; ++j) {
      junction->out[j]->a = total - junction->in[j]->b;
    }
#ifdef SELF_LOOP
    junction->self_loop->a = total - junction->self_loop->b;
#endif
/*
    if (i == 0) {
      result = total;
    }
	*/
  }

  // circulate the unit delays
  for (i = 0; i < unit->delay_n; ++i) {
    t_delay *delay = &delays[i];
    if (delay->invert) {
#ifdef RIMFILTER
	//one zero filter
      delay->b = ((0.0f - delay->a) + delay->c) * 0.5f;
	  delay->c = (0.0f - delay->a);
	  
	  //with variable a1
	  //delay->b = - delay->a + unit->a1*(delay->c + delay->a);
      //delay->c = (0.0f - delay->a);
	  
	  //one pole
	  //delay->c = delay->a * (1.f - unit->a1) + unit->a1 * delay->c;
	  //delay->b = - delay->c;
#else
      delay->b = 0.f - delay->a;
#endif
    }
    else {
      delay->b = delay->a;
    }
  }
  return(result);
}

////////////////////////////////////////////////////////////////////

void Membrane_init(Membrane* unit, int shape_type)
{
  t_shape *shape;
  int d = 0;
  int i = 0;
  int j = 0;
  // 1. Set calculation function -- always run at audio rate
  SETCALC(Membrane_next_a);

  int SIZ = IN0(5);
#ifdef EUCLIDEAN_POS
  unit->Size = (SIZ * 4) + 8;
#endif
  int doprint = IN0(7);
  
  // 2. Initialise state
#ifndef AUDIO_INPUT
  unit->triggered = 0;
  unit->excite = 0;
#endif
  unit->yj = 0;
  shape = unit->shape = getShape(shape_type, SIZ);//SHAPE_SZ);

  unit->delay_n = (shape->lines_n * 2)
#ifdef RIMGUIDES
    + shape->edge_n
#endif
#ifdef SELF_LOOP
    + shape->points_n
#endif
    ;
///////////////////////////////////////////printf
	if(doprint > 0){
		Print("Made shape with %d lines, %d points and %d edge points.\n",shape->lines_n,shape->points_n,shape->edge_n);
		Print("size: %d numinputs %d\n",SIZ,unit->mNumInputs);
		Print("edges are\n");
		for(int ii=0;ii<shape->points_n;ii++){
			if (shape->points[ii]->is_edge){
				Print("edge %d id:%d, x:%d, y:%d\n",ii,shape->points[ii]->id,shape->points[ii]->x,shape->points[ii]->y);
			}else{
				Print("center %d id:%d, x:%d, y:%d\n",ii,shape->points[ii]->id,shape->points[ii]->x,shape->points[ii]->y);
			}
		}
	}
///////////////////////////////
  // make all the junctions and delays
  unit->delays =
    (t_delay *) RTAlloc(unit->mWorld, unit->delay_n * sizeof(t_delay));
  memset((void *) unit->delays, 0, unit->delay_n * sizeof(t_delay));

  unit->junctions =
    (t_junction *) RTAlloc(unit->mWorld,
			   unit->shape->points_n * sizeof(t_junction)
			   );
  memset((void *) unit->junctions, 0,
	 unit->shape->points_n * sizeof(t_junction)
	 );

  for (i = 0; i < shape->lines_n; ++i) {
    t_line *line = shape->lines[i];
    t_junction *from, *to;
    t_delay *delay;

    from = &unit->junctions[line->a->id];
    to = &unit->junctions[line->b->id];

    // leftward delay
    delay = &unit->delays[d++];
    from->out[from->outs++] = delay;
    to->in[to->ins++] = delay;

    // rightward delay
    delay = &unit->delays[d++];
    from->in[from->ins++] = delay;
    to->out[to->outs++] = delay;
  }
#ifdef EUCLIDEAN_POS
  double euclidean_fac = sqrt(3);
#endif
  for (i = 0; i < shape->points_n; ++i) {
    t_point *point = shape->points[i];
    t_junction *junction = &unit->junctions[i];
#ifdef EUCLIDEAN_POS
    junction->x = point->x;
    junction->y = point->y*euclidean_fac;
#endif

#ifdef SELF_LOOP
    t_delay *delay = &unit->delays[d++];
    junction->self_loop = delay;
#endif

#ifdef RIMGUIDES
    assert((junction->ins < 6) == point->is_edge);
    if (point->is_edge) {
      t_delay *delay = &unit->delays[d++];
      delay->invert = 1;
      junction->out[junction->outs++] = delay;
      junction->in[junction->ins++] = delay;
    }
#endif
  }

  if(unit->mWorld->mVerbosity > 0){
    printf("%d delays initialised.\n", unit->delay_n);
  }

  // 3. Calculate one sample of output.
  // (why do this?)
  Membrane_next_a(unit, 1);
}


////////////////////////////////////////////////////////////////////

void Membrane_next_a(Membrane *unit, int inNumSamples) {
  // get the pointer to the output buffer
  float *out = OUT(0);
  int input_n = 0;
  // get the control rate input
#ifdef AUDIO_INPUT
  float *in = IN(input_n++);
#else
  float trigger = IN0(input_n++);
#endif
  float tension = IN0(input_n++);
  float loss = IN0(input_n++);
  
  float input_width_fac = sc_clip(IN0(3),0,1);
  float input_pos_fac = sc_clip(IN0(4),0,1);
  int input_pos = 0;
  int input_width = 0;
#ifdef EUCLIDEAN_POS
  if (input_pos_fac != unit->saved_input_pos || input_width_fac != unit->saved_input_width){
      unit->saved_input_pos = input_pos_fac;
	  unit->saved_input_width = input_width_fac;
	  input_pos = calc_input_fac(unit);
	  unit->input_pos = input_pos;
  }
#else
  input_width = ((float)unit->shape->points_n-1)*input_width_fac;
  input_pos = ((float)unit->shape->points_n-1)*input_pos_fac;
#endif
  
  unit->a1 = sc_clip(IN0(6),0,1);
  
  if (tension == 0) {
    // default tension
    tension =  0.0001;
  }
  unit->yj = 2.f * DELTA * DELTA / (tension * tension * GAMMA * GAMMA);
  float yj_r = 1.0f / unit->yj;

  if (loss > 1) {
    loss = 1.f;//0.99999;
  }
  unit->loss = loss;

#ifndef AUDIO_INPUT
  if (trigger >= 0.5 && (! unit->triggered)) {
    unit->triggered = 1;
    unit->excite = TRIGGER_DURATION;

  }
  else if (trigger < 0.5 && unit->triggered) {
    unit->triggered = 0;
  }
#endif

  for (int k=0; k < inNumSamples; ++k) {
    float input = 0.0;
#ifdef AUDIO_INPUT
    input = in[k];
#else
    if (unit->excite > 0) {
      input = (0.01 - (((float) rand() / RAND_MAX) * 0.02));
      unit->excite--;
    }
#endif
    out[k] = cycle(unit, input, yj_r,input_width,input_pos);
  }
}

void MembraneCircle_Ctor(Membrane* unit) {
  Membrane_init(unit, 0);
}

void MembraneHexagon_Ctor(Membrane* unit) {
  Membrane_init(unit, 1);
}

void Membrane_Dtor(Membrane* unit) {
  free_shape(unit->shape);
  RTFree(unit->mWorld, unit->delays);
  RTFree(unit->mWorld, unit->junctions);
}

////////////////////////////////////////////////////////////////////

// the load function is called by the host when the plug-in is loaded
PluginLoad(Membrane)
{
  ft = inTable;
  (*ft->fDefineUnit)("MembraneCircleV",
		     sizeof(Membrane),
		     (UnitCtorFunc)&MembraneCircle_Ctor,
		     (UnitDtorFunc)&Membrane_Dtor,
		     0);
  (*ft->fDefineUnit)("MembraneHexagonV",
		     sizeof(Membrane),
		     (UnitCtorFunc)&MembraneHexagon_Ctor,
		     (UnitDtorFunc)&Membrane_Dtor,
		     0);
}

////////////////////////////////////////////////////////////////////
