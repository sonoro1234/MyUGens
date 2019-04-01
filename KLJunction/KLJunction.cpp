/*
	SuperCollider real time audio synthesis system
 Copyright (c) 2002 James McCartney. All rights reserved.
	http://www.audiosynth.com

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */

#include "SC_PlugIn.h"
#include "DWG.cpp"
#include <cstdio>
#define SCWrapClassINI(classname) extern "C" {\
	void classname##_Ctor(classname* unit);\
	void classname##_next(classname *unit, int inNumSamples);\
	void classname##_Dtor(classname *unit);\
}\
void classname##_Ctor(classname* unit){new(unit) classname(unit);SETCALC(classname##_next);(unit->mCalcFunc)(unit, 1);}\
void classname##_Dtor(classname* unit){unit->~classname();}
InterfaceTable *ft;

////////////////////////////////////////////////////////////
struct HumanVdel : public Unit
{
	int numtubes;
	int numdels;
    float Hgo,Hlo,Hgi,Hli;
    float rg,rl,loss;
	//convenience variables for copying particular input data
	float * areas;
	//float * scattering;
	float * dels;
	TUBE ** tubes;
    KL_Junction ** kl_j;
	HumanVdel(HumanVdel* unit);
	~HumanVdel();
};
SCWrapClassINI(HumanVdel);
/*
extern "C" {

	void HumanVdel_next(HumanVdel *unit, int inNumSamples);
	void HumanVdel_Ctor(HumanVdel* unit);
	void HumanVdel_Dtor(HumanVdel* unit);

}
*/


////output= HumanV.ar(input, loss,rg,rl,dels, areas );
HumanVdel::HumanVdel(HumanVdel* unit){
//void HumanVdel_Ctor(HumanVdel* unit) {

	int i,j;
    unit->Hgo = unit->Hlo = unit->Hgi = unit->Hli =0.0;

	unit->loss= ZIN0(1);
    unit->rg= ZIN0(2);
    unit->rl= ZIN0(3);

	unit->numdels = ZIN0(4);
	int numtubes = unit->numtubes = ZIN0(5 + unit->numdels);
	
	if(numtubes<2) {
		printf("too few tubes! only %d \n", numtubes);
		return;
	}
	
	unit->tubes = (TUBE **)RTAlloc(unit->mWorld, (unit->numtubes)*sizeof(TUBE*));
    for (i=0; i<(unit->numtubes); ++i){
        unit->tubes[i] = new(unit) TUBE(unit);
    }
	
    unit->kl_j = (KL_Junction **)RTAlloc(unit->mWorld, (numtubes-1)*sizeof(KL_Junction*));
    for (i=0; i<(numtubes-1); ++i){
        unit->kl_j[i] = new(unit) KL_Junction(unit);
    }
    
	unit->dels= (float*)RTAlloc(unit->mWorld, (unit->numdels) * sizeof(float));
	for(i=0;i<unit->numdels;i++)
		unit->dels[i] = sc_max(ZIN0(5 + i),1);
	
    unit->areas= (float*)RTAlloc(unit->mWorld, (numtubes) * sizeof(float));
	for(i=0;i<numtubes;i++)
        unit->areas[i] = ZIN0(6 + unit->numdels + i);
	
	
	//SETCALC(HumanVdel_next);
}

HumanVdel::~HumanVdel(){
//void HumanVdel_Dtor(HumanVdel* unit) {

	int i;
	for (i=0; i<(numtubes); ++i){
        delete tubes[i];
    }
    RTFree(mWorld,tubes);
	
    for (i=0; i<(numtubes-1); ++i){
        delete kl_j[i];
    }
    RTFree(mWorld,kl_j);
   
    RTFree(mWorld, areas);
	RTFree(mWorld, dels);
}

void Calc_Scattering(HumanVdel *unit){
    //float *k = unit->scattering;
    float *areas = unit->areas;
    for(int i=0;i<unit->numtubes-1;i++){
        //k[i] = getK(areas[i],areas[i+1]);
        unit->kl_j[i]->set_areas(areas[i],areas[i+1]);
    }
}
void HumanVdel_next(HumanVdel *unit, int inNumSamples) {

	int i,j;
    KL_Junction **kl_j = unit->kl_j;
	TUBE **tubes = unit->tubes;
	int numtubes= unit->numtubes;
	int numdels = unit->numdels;
	
	//value to store
	float * in= IN(0);
	float * out= OUT(0);

	unit->loss= ZIN0(1);
    unit->rg= ZIN0(2);
    unit->rl= ZIN0(3);

	for(i=0;i<numtubes;i++)
        unit->areas[i] = ZIN0(6 + unit->numdels + i);
    
	for(i=0;i<numdels;i++)
		unit->dels[i] = sc_max(ZIN0(5 + i),1);
	
	 //set tubes
    for (j=0; j<numtubes; ++j) {
            unit->tubes[j]->del = unit->dels[j] - 1.0;
            unit->tubes[j]->set_area(unit->areas[j], unit->loss);
	}
	
    Calc_Scattering(unit);


	for (i=0; i<inNumSamples; ++i) {

        //update ins
        unit->Hgi = in[i] + unit->rg * tubes[0]->outL;
		
        //fill tubes
        tubes[0]->inR = unit->Hgo;
        tubes[0]->inL = kl_j[0]->outL;
        
        for (j=1; j<numtubes-1; ++j) {
            tubes[j]->inR = kl_j[j-1]->outR;
            tubes[j]->inL = kl_j[j]->outL;
		}
		
		tubes[numtubes-1]->inR = kl_j[unit->numtubes-2]->outR;
        tubes[numtubes-1]->inL = unit->Hlo*unit->rl;
		
		//advance tubes
		for (j=0; j<numtubes; ++j) {
            tubes[j]->go();
		}
		
		 //fill junctions area1
        for (j=0; j< numtubes-1; ++j) {
            kl_j[j]->inR = tubes[j]->outR;
            kl_j[j]->inL = tubes[j + 1]->outL;
        }
		
        //kl_j[numtubes-2]->inR = tubes[numtubes-2]->outR;
       // kl_j[numtubes-2]->inL = unit->Hlo*unit->rl;
        
        unit->Hli = tubes[numtubes-1]->outR;
        //compute outs
        for (j=0; j<(numtubes-1); ++j) {
            kl_j[j]->go();
		}
        unit->Hgo = unit->Hgi;
        unit->Hlo = unit->Hli;
        
        out[i] = unit->Hlo;
	}
}

////////////////////////////////////////////////////////////
struct HumanV : public Unit
{
	int numtubes;
    float Hgo,Hlo,Hgi,Hli;
    float rg,rl,loss;
	//convenience variables for copying particular input data
	float * areas;
	float * scattering;

    KL_Junction ** kl_j;
};


extern "C" {

	void HumanV_next(HumanV *unit, int inNumSamples);
	void HumanV_Ctor(HumanV* unit);
	void HumanV_Dtor(HumanV* unit);

}


////output= HumanV.ar(input, loss,rg,rl, areas );
void HumanV_Ctor(HumanV* unit) {

	int i,j;
    unit->Hgo = unit->Hlo = unit->Hgi = unit->Hli =0.0;
	int numinputs = unit->mNumInputs;
	int numtubes= numinputs-5;
	unit->numtubes= numtubes;

	if(numtubes<2) {
		printf("too few tubes! only %d \n", numtubes);
		return;
	}

    unit->kl_j = (KL_Junction **)RTAlloc(unit->mWorld, (numtubes-1)*sizeof(KL_Junction*));
    for (i=0; i<(numtubes-1); ++i){
        unit->kl_j[i] = new(unit) KL_Junction(unit);
    }
    
	unit->loss= ZIN0(1);
    unit->rg= ZIN0(2);
    unit->rl= ZIN0(3);
    
    unit->areas= (float*)RTAlloc(unit->mWorld, (numtubes) * sizeof(float));
	unit->scattering= (float*)RTAlloc(unit->mWorld, (numtubes-1) * sizeof(float));
	for(i=0;i<numtubes;i++)
        unit->areas[i] = ZIN0(i + 4);
    //for(i=0;i<numtubes-1;i++)
        //unit->scattering[i] = ZIN0(i + 5);
	SETCALC(HumanV_next);
	(unit->mCalcFunc)(unit, 1);
}

void HumanV_Dtor(HumanV* unit) {

	int i;

    for (i=0; i<(unit->numtubes-1); ++i){
        delete unit->kl_j[i];
    }
    RTFree(unit->mWorld,unit->kl_j);
   
	RTFree(unit->mWorld, unit->scattering);
    RTFree(unit->mWorld, unit->areas);
}

void Calc_Scattering(HumanV *unit){
    //float *k = unit->scattering;
    float *areas = unit->areas;
    for(int i=0;i<unit->numtubes-1;i++){
        //k[i] = getK(areas[i],areas[i+1]);
        unit->kl_j[i]->set_areas(areas[i],areas[i+1]);
    }
}
void HumanV_next(HumanV *unit, int inNumSamples) {

	int i,j;
    KL_Junction **kl_j = unit->kl_j;
	int numtubes= unit->numtubes;

	//value to store
	float * in= IN(0);
	float * out= OUT(0);



	unit->loss= ZIN0(1);
    unit->rg= ZIN0(2);
    unit->rl= ZIN0(3);

	for(i=0;i<numtubes;i++)
        unit->areas[i] = ZIN0(i + 4);
    
    Calc_Scattering(unit);

    /*
    for (i=0; i<(numtubes-1); ++i) {
		kl_j[i]->k = unit->scattering[i];
		kl_j[i]->lossL = unit->loss;
        kl_j[i]->lossR = unit->loss;
	}
    */
	for (i=0; i<inNumSamples; ++i) {

		
        //update ins
        unit->Hgi = in[i] + unit->rg * kl_j[0]->outL;
        
        kl_j[0]->inR = unit->Hgo;//in[i] + losses[0]*leftouts[0];;
        kl_j[0]->inL = kl_j[1]->outL;
        
        for (j=1; j<(numtubes-2); ++j) {
            kl_j[j]->inR = kl_j[j-1]->outR;
            kl_j[j]->inL = kl_j[j+1]->outL;
		}
        kl_j[numtubes-2]->inR = kl_j[numtubes-3]->outR;
        kl_j[numtubes-2]->inL = unit->Hlo*unit->rl;
        
        unit->Hli = kl_j[numtubes-2]->outR;
        //compute outs
        for (j=0; j<(numtubes-1); ++j) {
            kl_j[j]->go();
		}
        unit->Hgo = unit->Hgi;
        unit->Hlo = unit->Hli;
        
        out[i] = unit->Hlo;
	}
}
////////////////////////////////////////////////////////////
struct HumanVN : public Unit
{
	int numtubes;
    float Hgo,Hlo,Hgi,Hli,Hno,Hni;
    float rg,rl,rn,loss,lmix,nmix;
    int area1len,areasNlen;
	//convenience variables for copying particular input data
	float * areas;
	float * scattering;
    float *areasN;
	float *scatteringN;
    KL_Junction ** kl_j;
    KL_Junction ** kl_jN;
    NJunction<3> NoseJunction;//();//{};
};


extern "C" {

	void HumanVN_next(HumanVN *unit, int inNumSamples);
	void HumanVN_Ctor(HumanVN* unit);
	void HumanVN_Dtor(HumanVN* unit);

}


////output= HumanVN.ar(input, loss,rg,rl,rn,area1len, areas,areasN );
//(input, loss,rg,rl,rn,lmix,nmix,area1len,numtubes, areas,areasNlen,areasN );
#define AREASZIN 9
void HumanVN_Ctor(HumanVN* unit) {

	int i,j;
    unit->Hgo = unit->Hlo = unit->Hgi = unit->Hli = unit->Hno = unit->Hni=0.0;

	unit->loss= ZIN0(1);
    unit->rg= ZIN0(2);
    unit->rl= ZIN0(3);
    unit->rn= ZIN0(4);
    unit->lmix= ZIN0(5);
    unit->nmix= ZIN0(6);
    unit->area1len = ZIN0(7);
    unit->numtubes = ZIN0(8);
    if(unit->numtubes<2) {
		printf("too few tubes! only %d \n", unit->numtubes);
		return;
	}
    
    unit->areas= (float*)RTAlloc(unit->mWorld, (unit->numtubes) * sizeof(float));
	unit->scattering= (float*)RTAlloc(unit->mWorld, (unit->numtubes-2) * sizeof(float));
    
    unit->kl_j = (KL_Junction **)RTAlloc(unit->mWorld, (unit->numtubes-2)*sizeof(KL_Junction*));
    for (i=0; i<(unit->numtubes-2); ++i){
        unit->kl_j[i] = new(unit) KL_Junction(unit);
    }
    
    unit->areasNlen = ZIN0(unit->numtubes + AREASZIN);
    
    unit->areasN = (float*)RTAlloc(unit->mWorld, (unit->areasNlen) * sizeof(float));
	unit->scatteringN= (float*)RTAlloc(unit->mWorld, (unit->areasNlen-1) * sizeof(float));
    
    unit->kl_jN = (KL_Junction **)RTAlloc(unit->mWorld, (unit->areasNlen-1)*sizeof(KL_Junction*));
    for (i=0; i<(unit->areasNlen-1); ++i){
        unit->kl_jN[i] = new(unit) KL_Junction(unit);
    }
    unit->NoseJunction.Init();
	SETCALC(HumanVN_next);
	(unit->mCalcFunc)(unit, 1);
}

void HumanVN_Dtor(HumanVN* unit) {

	int i;

    for (i=0; i<(unit->numtubes-2); ++i){
        delete unit->kl_j[i];
    }
    for (i=0; i<(unit->areasNlen-1); ++i){
        delete unit->kl_jN[i];
    }
    RTFree(unit->mWorld,unit->kl_j);
    RTFree(unit->mWorld,unit->kl_jN);
    
	RTFree(unit->mWorld, unit->scattering);
    RTFree(unit->mWorld, unit->areas);
    
    RTFree(unit->mWorld, unit->scatteringN);
    RTFree(unit->mWorld, unit->areasN);
}

void Calc_Scattering(HumanVN *unit){
    //float *k = unit->scattering;
    float *areas = unit->areas;
    KL_Junction **kl_j = unit->kl_j;
    
    for(int i = 0;i<unit->area1len-1;i++){
        //k[i] = getK(areas[i],areas[i+1]);
        kl_j[i]->set_areas(areas[i],areas[i+1]);
    }
    for(int i=unit->area1len;i<unit->numtubes-1;i++){
        //k[i-1] = getK(areas[i],areas[i+1]);
        kl_j[i-1]->set_areas(areas[i],areas[i+1]);
    }
    

    //float *kN = unit->scatteringN;
    float *areasN = unit->areasN;
    KL_Junction **kl_jN = unit->kl_jN;
    for(int i=0;i<unit->areasNlen -1;i++){
        //kN[i] = getK(areasN[i],areasN[i+1]);//(areasN[i] - areasN[i+1])/(areasN[i] + areasN[i+1]);
        kl_jN[i]->set_areas(areasN[i],areasN[i+1]);
    }
    

    float NJareas[3];
    NJareas[0] = areas[unit->area1len-1];
    NJareas[1] = areas[unit->area1len];
    NJareas[2] = areasN[0];
    unit->NoseJunction.calc_alpha(NJareas);
}
void HumanVN_next(HumanVN *unit, int inNumSamples) {

	int i,j;
    KL_Junction **kl_j = unit->kl_j;
    KL_Junction **kl_jN = unit->kl_jN;
	int numtubes= unit->numtubes;

	//value to store
	float * in= IN(0);
	float * out= OUT(0);


	unit->loss= ZIN0(1);
    unit->rg= ZIN0(2);
    unit->rl= ZIN0(3);
    unit->rn= ZIN0(4);
    unit->lmix= ZIN0(5);
    unit->nmix= ZIN0(6);
    unit->area1len = ZIN0(7);
    
    for(i=0;i<numtubes;i++)
        unit->areas[i] = sc_max(ZIN0(i + AREASZIN),0.0);
    
    for(i=0;i<unit->areasNlen;i++)
        unit->areasN[i] = sc_max(ZIN0(i + unit->numtubes + AREASZIN + 1),0.0);
    
    Calc_Scattering(unit);

    
	for (i=0; i<inNumSamples; ++i) {

		
        //update ins
        //area1
        unit->Hgi = in[i] + unit->rg * kl_j[0]->outL;
        
        kl_j[0]->inR = unit->Hgo;//in[i] + losses[0]*leftouts[0];;
        kl_j[0]->inL = kl_j[1]->outL;
        
        for (j=1; j<(unit->area1len-2); ++j) {
            kl_j[j]->inR = kl_j[j-1]->outR;
            kl_j[j]->inL = kl_j[j+1]->outL;
		}
        
        kl_j[unit->area1len-2]->inR = kl_j[unit->area1len-3]->outR;
        kl_j[unit->area1len-2]->inL = unit->NoseJunction.out[0];
        //nose
        unit->NoseJunction.in[0] = kl_j[unit->area1len-2]->outR;
        unit->NoseJunction.in[1] = kl_j[unit->area1len-1]->outL;
        unit->NoseJunction.in[2] = kl_jN[0]->outL;
        //area2
        kl_j[unit->area1len-1]->inR = unit->NoseJunction.out[1];
        kl_j[unit->area1len-1]->inL = kl_j[unit->area1len]->outL;
        
        for (j=unit->area1len; j<(unit->numtubes-3); ++j) {
            kl_j[j]->inR = kl_j[j-1]->outR;
            kl_j[j]->inL = kl_j[j+1]->outL;
		}
        
        kl_j[unit->numtubes-3]->inR = kl_j[unit->numtubes-4]->outR;
        kl_j[unit->numtubes-3]->inL = unit->Hlo*unit->rl;
        
        unit->Hli = kl_j[numtubes-3]->outR;
        //noseareas
        kl_jN[0]->inR = unit->NoseJunction.out[2];
        kl_jN[0]->inL = kl_jN[1]->outL;
        for (j=1; j<(unit->areasNlen-2); ++j) {
            kl_jN[j]->inR = kl_jN[j-1]->outR;
            kl_jN[j]->inL = kl_jN[j+1]->outL;
		}
        kl_jN[unit->areasNlen-2]->inR = kl_jN[unit->areasNlen-3]->outR;
        kl_jN[unit->areasNlen-2]->inL = unit->Hno*unit->rn;
        
        unit->Hni = kl_jN[unit->areasNlen-2]->outR;
        //compute outs
        for (j=0; j<(numtubes-2); ++j) {
            kl_j[j]->go();
		}
        for (j=0; j<(unit->areasNlen-1); ++j) {
            kl_jN[j]->go();
		}
        unit->NoseJunction.go();
        
        unit->Hgo = unit->Hgi;
        unit->Hlo = unit->Hli;
        unit->Hno = unit->Hni;
        
        out[i] = unit->Hlo*unit->lmix + unit->Hno*unit->nmix;
	}
}
////////////////////////////////////////////////////////////
struct HumanVNdel : public Unit
{
	HumanVNdel(HumanVNdel *unit);
    ~HumanVNdel();
    int numtubes,numdels;
    float Hgo,Hlo,Hgi,Hli,Hno,Hni;
    float rg,rl,rn,loss,lmix,nmix;
    int area1len,areasNlen;
	//convenience variables for copying particular input data
	float * areas;
    float *dels;
	//float * scattering;
    float *areasN;
	float *scatteringN;
    //float del;
    KL_Junction ** kl_j;
    KL_Junction ** kl_jN;
    TUBE ** tubes;
    TUBE ** tubesN;
    NJunction<3> NoseJunction;//();//{};
};

SCWrapClassINI(HumanVNdel);
/*
extern "C" {

	void HumanVNdel_next(HumanVNdel *unit, int inNumSamples);
	void HumanVNdel_Ctor(HumanVNdel* unit);
	void HumanVNdel_Dtor(HumanVNdel* unit);

}
*/

////output= HumanVNdel.ar(input, loss,rg,rl,rn,area1len, areas,areasN );
//(input, loss,rg,rl,rn,lmix,nmix,area1len,numtubes, areas,areasNlen,areasN );
//#define AREASZINdel 12
HumanVNdel::HumanVNdel(HumanVNdel *unit){
//void HumanVNdel_Ctor(HumanVNdel* unit) {

	int i,j;
    unit->Hgo = unit->Hlo = unit->Hgi = unit->Hli = unit->Hno = unit->Hni=0.0;

	unit->loss= ZIN0(1);
    unit->rg= ZIN0(2);
    unit->rl= ZIN0(3);
    unit->rn= ZIN0(4);
    unit->lmix= ZIN0(5);
    unit->nmix= ZIN0(6);
    unit->area1len = ZIN0(7);
    unit->numdels = ZIN0(10);
    unit->numtubes = ZIN0(11 + unit->numdels);
    if(unit->numtubes<2) {
		printf("too few tubes! only %d \n", unit->numtubes);
		return;
	}
    unit->dels= (float*)RTAlloc(unit->mWorld, (unit->numdels) * sizeof(float));
    
    unit->areas= (float*)RTAlloc(unit->mWorld, (unit->numtubes) * sizeof(float));
	//unit->scattering= (float*)RTAlloc(unit->mWorld, (unit->numtubes-2) * sizeof(float));
    
    unit->tubes = (TUBE **)RTAlloc(unit->mWorld, (unit->numtubes)*sizeof(TUBE*));
    for (i=0; i<(unit->numtubes); ++i){
        unit->tubes[i] = new(unit) TUBE(unit);
    }
    
    unit->kl_j = (KL_Junction **)RTAlloc(unit->mWorld, (unit->numtubes-2)*sizeof(KL_Junction*));
    for (i=0; i<(unit->numtubes-2); ++i){
        unit->kl_j[i] = new(unit) KL_Junction(unit);
    }
    
    unit->areasNlen = ZIN0(12 + unit->numtubes + unit->numdels);
    
    unit->areasN = (float*)RTAlloc(unit->mWorld, (unit->areasNlen) * sizeof(float));
	unit->scatteringN= (float*)RTAlloc(unit->mWorld, (unit->areasNlen-1) * sizeof(float));
    
    unit->tubesN = (TUBE **)RTAlloc(unit->mWorld, (unit->areasNlen)*sizeof(TUBE*));
    for (i=0; i<unit->areasNlen; ++i){
        unit->tubesN[i] = new(unit) TUBE(unit);
    }
    
    unit->kl_jN = (KL_Junction **)RTAlloc(unit->mWorld, (unit->areasNlen-1)*sizeof(KL_Junction*));
    for (i=0; i<(unit->areasNlen-1); ++i){
        unit->kl_jN[i] = new(unit) KL_Junction(unit);
    }
    
    //Print("%d,%d,%d\n",unit->numdels,unit->numtubes,unit->areasNlen);
    unit->NoseJunction.Init();
	//SETCALC(HumanVNdel_next);
	//(unit->mCalcFunc)(unit, 1);
}

//void HumanVNdel_Dtor(HumanVNdel* unit) {
HumanVNdel::~HumanVNdel(){
	int i;

    for (i=0; i<(numtubes-2); ++i){
        delete kl_j[i];
    }
    for (i=0; i<(areasNlen-1); ++i){
        delete kl_jN[i];
    }
    RTFree(mWorld,kl_j);
    RTFree(mWorld,kl_jN);
    
    for (i=0; i<(numtubes); ++i){
        delete tubes[i];
    }
    for (i=0; i<(areasNlen); ++i){
        delete tubesN[i];
    }
    
    RTFree(mWorld,tubes);
    RTFree(mWorld,tubesN);
    
	//RTFree(mWorld, scattering);
    RTFree(mWorld, areas);
    RTFree(mWorld, dels);
    
    RTFree(mWorld, scatteringN);
    RTFree(mWorld, areasN);
}

void Calc_Scattering(HumanVNdel *unit){
    //float *k = unit->scattering;
    float *areas = unit->areas;
    KL_Junction **kl_j = unit->kl_j;
    
    for(int i = 0;i<unit->area1len-1;i++){
        //k[i] = getK(areas[i],areas[i+1]);
        kl_j[i]->set_areas(areas[i],areas[i+1]);
    }
    for(int i=unit->area1len;i<unit->numtubes-1;i++){
        //k[i-1] = getK(areas[i],areas[i+1]);
        kl_j[i-1]->set_areas(areas[i],areas[i+1]);
    }
    
    //float *kN = unit->scatteringN;
    float *areasN = unit->areasN;
    KL_Junction **kl_jN = unit->kl_jN;
    for(int i=0;i<unit->areasNlen -1;i++){
        //kN[i] = getK(areasN[i],areasN[i+1]);//(areasN[i] - areasN[i+1])/(areasN[i] + areasN[i+1]);
        kl_jN[i]->set_areas(areasN[i],areasN[i+1]);
    }
    

    float NJareas[3];
    NJareas[0] = areas[unit->area1len-1];
    NJareas[1] = areas[unit->area1len];
    NJareas[2] = areasN[0];
    unit->NoseJunction.calc_alpha(NJareas);
}
void HumanVNdel_next(HumanVNdel *unit, int inNumSamples) {

	int i,j;
    KL_Junction **kl_j = unit->kl_j;
    KL_Junction **kl_jN = unit->kl_jN;
	int numtubes= unit->numtubes;

	//value to store
	float * in= IN(0);
	float * out= OUT(0);


	unit->loss= ZIN0(1);
    unit->rg= ZIN0(2);
    unit->rl= ZIN0(3);
    unit->rn= ZIN0(4);
    unit->lmix= ZIN0(5);
    unit->nmix= ZIN0(6);
    unit->area1len = ZIN0(7);
    //unit->del = sc_max(ZIN0(8),1.0);
    
    float * in_noise = IN(8);
    int noiseloc = ZIN0(9);
    
    //get dels
    for(i=0;i<unit->numdels;i++)
        unit->dels[i] = sc_max(ZIN0(11 + i),1.0);
    
    //get areas
    for(i=0;i<numtubes;i++)
        unit->areas[i] = sc_max(ZIN0(12 + unit->numdels +i),0.0);
    
    for(i=0;i<unit->areasNlen;i++)
        unit->areasN[i] = sc_max(ZIN0(13 + unit->numtubes + unit->numdels + i),0.0);
    
    //set tubes
    for (j=0; j<numtubes; ++j) {
            unit->tubes[j]->del = unit->dels[j] - 1.0;
            unit->tubes[j]->set_area(unit->areas[j], unit->loss);
	}
    for (j=0; j<(unit->areasNlen); ++j) {
            unit->tubesN[j]->del = unit->dels[0] - 1.0;
            unit->tubesN[j]->set_area(unit->areasN[j], unit->loss);
	}
    Calc_Scattering(unit);

    
	for (i=0; i<inNumSamples; ++i) {

		
        //update ins
        unit->Hgi = in[i] + unit->rg * unit->tubes[0]->outL;
        //fill tubes area1
        unit->tubes[0]->inR = unit->Hgo;
        unit->tubes[0]->inL = kl_j[0]->outL;        
        for (j=1; j<(unit->area1len-1); ++j) {
            unit->tubes[j]->inR = kl_j[j-1]->outR;
            unit->tubes[j]->inL = kl_j[j]->outL;
		}
        unit->tubes[unit->area1len-1]->inR = kl_j[unit->area1len-2]->outR;
        unit->tubes[unit->area1len-1]->inL = unit->NoseJunction.out[0];
        //fill tubes area2
        unit->tubes[unit->area1len]->inR = unit->NoseJunction.out[1];
        unit->tubes[unit->area1len]->inL = kl_j[unit->area1len -1]->outL;
        for (j=unit->area1len + 1; j<(unit->numtubes-1); ++j) {
            unit->tubes[j]->inR = kl_j[j-2]->outR;
            unit->tubes[j]->inL = kl_j[j-1]->outL;
		}
        unit->tubes[unit->numtubes-1]->inR = kl_j[unit->numtubes-3]->outR;
        unit->tubes[unit->numtubes-1]->inL = unit->Hlo*unit->rl;
        
        //fill noise
        if(noiseloc > 0 && noiseloc <= numtubes)
            unit->tubes[noiseloc-1]->inR += in_noise[i];
        
        //advance tubes area1 and 2
        for (j=0; j<numtubes; ++j) {
            unit->tubes[j]->go();
		}
        //fill tubesN
        unit->tubesN[0]->inR = unit->NoseJunction.out[2];
        unit->tubesN[0]->inL = kl_jN[0]->outL;
        for (j=1; j<unit->areasNlen-1; ++j) {
            unit->tubesN[j]->inR= kl_jN[j-1]->outR;
            unit->tubesN[j]->inL= kl_jN[j]->outL;
		}
        unit->tubesN[unit->areasNlen-1]->inR = kl_jN[unit->areasNlen-2]->outR;
        unit->tubesN[unit->areasNlen-1]->inL = unit->Hno*unit->rn;
        
        //advance tubesN
        for (j=0; j<(unit->areasNlen); ++j) {
            unit->tubesN[j]->go();
		}
        //fill junctions area1
        for (j=0; j<(unit->area1len-1); ++j) {
            kl_j[j]->inR = unit->tubes[j]->outR;
            kl_j[j]->inL = unit->tubes[j + 1]->outL;
        }
        //fill nose
        unit->NoseJunction.in[0] = unit->tubes[unit->area1len-1]->outR;
        unit->NoseJunction.in[1] = unit->tubes[unit->area1len]->outL;
        unit->NoseJunction.in[2] = unit->tubesN[0]->outL;
        
        //fill junctions area2
        for (j=unit->area1len-1; j<(unit->numtubes-2); ++j) {
            kl_j[j]->inR = unit->tubes[j + 1]->outR;
            kl_j[j]->inL = unit->tubes[j + 2]->outL;
        }
        
        unit->Hli = unit->tubes[numtubes-1]->outR;
        ///////noseareas
        
        //fill junctionsN
        for (j=0; j<(unit->areasNlen-1); ++j) {
            kl_jN[j]->inR = unit->tubesN[j]->outR;
            kl_jN[j]->inL = unit->tubesN[j + 1]->outL;
		}
        
        unit->Hni = unit->tubesN[unit->areasNlen-1]->outR;
        //compute outs

        
        for (j=0; j<(numtubes-2); ++j) {
            kl_j[j]->go();
		}
        for (j=0; j<(unit->areasNlen-1); ++j) {
            kl_jN[j]->go();
		}
        unit->NoseJunction.go();
        
        unit->Hgo = unit->Hgi;
        unit->Hlo = unit->Hli;
        unit->Hno = unit->Hni;
        
        out[i] = unit->Hlo*unit->lmix + unit->Hno*unit->nmix;
	}
}
//////////////////////////////////////////
struct HumanVNdelO2 : public HumanVNdel{
    HumanVNdelO2(HumanVNdelO2 *unit);
    ~HumanVNdelO2();
    float lastin,lastin_noise;
};
SCWrapClassINI(HumanVNdelO2);
HumanVNdelO2::HumanVNdelO2(HumanVNdelO2 *unit):HumanVNdel(unit){
    lastin_noise = lastin = 0.0;
    //SETCALC(HumanVNdelO2_next);
	//(unit->mCalcFunc)(unit, 1);
}
HumanVNdelO2::~HumanVNdelO2(){
}
void HumanVNdelO2_next(HumanVNdelO2 *unit, int inNumSamples) {

	int i,j;
    KL_Junction **kl_j = unit->kl_j;
    KL_Junction **kl_jN = unit->kl_jN;
	int numtubes= unit->numtubes;

	//value to store
	float * in= IN(0);
	float * out= OUT(0);


	unit->loss= ZIN0(1);
    unit->rg= ZIN0(2);
    unit->rl= ZIN0(3);
    unit->rn= ZIN0(4);
    unit->lmix= ZIN0(5);
    unit->nmix= ZIN0(6);
    unit->area1len = ZIN0(7);
    //unit->del = sc_max(ZIN0(8),1.0);
    
    float * in_noise = IN(8);
    int noiseloc = sc_min(ZIN0(9),numtubes);
    
    //get dels
    for(i=0;i<unit->numdels;i++)
        unit->dels[i] = sc_max(ZIN0(11 + i),1.0);
    
    //get areas
    for(i=0;i<numtubes;i++)
        unit->areas[i] = sc_max(ZIN0(12 + unit->numdels +i),0.0);
    
    for(i=0;i<unit->areasNlen;i++)
        unit->areasN[i] = sc_max(ZIN0(13 + unit->numtubes + unit->numdels + i),0.0);
    
    //set tubes
    for (j=0; j<numtubes; ++j) {
            unit->tubes[j]->del = unit->dels[j] - 1.0;
            unit->tubes[j]->set_area(unit->areas[j], unit->loss);
	}
    for (j=0; j<(unit->areasNlen); ++j) {
            unit->tubesN[j]->del = unit->dels[0] - 1.0;
            unit->tubesN[j]->set_area(unit->areasN[j], unit->loss);
	}
    Calc_Scattering(unit);

    float inp;
    float inpnoise;
	for (i=0; i<inNumSamples; ++i) {
        for(int ii=0;ii<2;ii++){
            if (ii==0){
                inp = (unit->lastin + in[i])*0.5;
                inpnoise = (unit->lastin_noise + in_noise[i])*0.5;
            }else{
                unit->lastin = inp = in[i];
                unit->lastin_noise = inpnoise = in_noise[i];
            }
            
        //update ins
        unit->Hgi = inp + unit->rg * unit->tubes[0]->outL;
        //fill tubes area1
        unit->tubes[0]->inR = unit->Hgo;
        unit->tubes[0]->inL = kl_j[0]->outL;        
        for (j=1; j<(unit->area1len-1); ++j) {
            unit->tubes[j]->inR = kl_j[j-1]->outR;
            unit->tubes[j]->inL = kl_j[j]->outL;
		}
        unit->tubes[unit->area1len-1]->inR = kl_j[unit->area1len-2]->outR;
        unit->tubes[unit->area1len-1]->inL = unit->NoseJunction.out[0];
        //fill tubes area2
        unit->tubes[unit->area1len]->inR = unit->NoseJunction.out[1];
        unit->tubes[unit->area1len]->inL = kl_j[unit->area1len -1]->outL;
        for (j=unit->area1len + 1; j<(unit->numtubes-1); ++j) {
            unit->tubes[j]->inR = kl_j[j-2]->outR;
            unit->tubes[j]->inL = kl_j[j-1]->outL;
		}
        unit->tubes[unit->numtubes-1]->inR = kl_j[unit->numtubes-3]->outR;
        unit->tubes[unit->numtubes-1]->inL = unit->Hlo*unit->rl;
        
        //fill noise
        if(noiseloc > 0){// && noiseloc <= numtubes){
            unit->tubes[noiseloc-1]->inR += inpnoise;
            unit->tubes[noiseloc-1]->inL += inpnoise;
            //unit->tubes[noiseloc]->inR += inpnoise;
            //unit->tubes[noiseloc]->inL += inpnoise;
        }
        //advance tubes area1 and 2
        for (j=0; j<numtubes; ++j) {
            unit->tubes[j]->go();
		}
        //fill tubesN
        unit->tubesN[0]->inR = unit->NoseJunction.out[2];
        unit->tubesN[0]->inL = kl_jN[0]->outL;
        for (j=1; j<unit->areasNlen-1; ++j) {
            unit->tubesN[j]->inR= kl_jN[j-1]->outR;
            unit->tubesN[j]->inL= kl_jN[j]->outL;
		}
        unit->tubesN[unit->areasNlen-1]->inR = kl_jN[unit->areasNlen-2]->outR;
        unit->tubesN[unit->areasNlen-1]->inL = unit->Hno*unit->rn;
        
        //advance tubesN
        for (j=0; j<(unit->areasNlen); ++j) {
            unit->tubesN[j]->go();
		}
        //fill junctions area1
        for (j=0; j<(unit->area1len-1); ++j) {
            kl_j[j]->inR = unit->tubes[j]->outR;
            kl_j[j]->inL = unit->tubes[j + 1]->outL;
        }
        //fill nose
        unit->NoseJunction.in[0] = unit->tubes[unit->area1len-1]->outR;
        unit->NoseJunction.in[1] = unit->tubes[unit->area1len]->outL;
        unit->NoseJunction.in[2] = unit->tubesN[0]->outL;
        
        //fill junctions area2
        for (j=unit->area1len-1; j<(unit->numtubes-2); ++j) {
            kl_j[j]->inR = unit->tubes[j + 1]->outR;
            kl_j[j]->inL = unit->tubes[j + 2]->outL;
        }
        
        unit->Hli = unit->tubes[numtubes-1]->outR;
        ///////noseareas
        
        //fill junctionsN
        for (j=0; j<(unit->areasNlen-1); ++j) {
            kl_jN[j]->inR = unit->tubesN[j]->outR;
            kl_jN[j]->inL = unit->tubesN[j + 1]->outL;
		}
        
        unit->Hni = unit->tubesN[unit->areasNlen-1]->outR;
        //compute outs

        
        for (j=0; j<(numtubes-2); ++j) {
            kl_j[j]->go();
		}
        for (j=0; j<(unit->areasNlen-1); ++j) {
            kl_jN[j]->go();
		}
        unit->NoseJunction.go();
        
        unit->Hgo = unit->Hgi;
        unit->Hlo = unit->Hli;
        unit->Hno = unit->Hni;
        }
        out[i] = unit->Hlo*unit->lmix + unit->Hno*unit->nmix;
	}
}
////////////////////////////////////////////////////////////
struct HumanVNdelU : public Unit
{
	int numtubes;
    float Hgo,Hlo,Hgi,Hli,Hno,Hni;
    float rg,rl,rn,loss,lmix,nmix;
    int area1len,areasNlen;
	//convenience variables for copying particular input data
	float * areas;
	float * scattering;
    float *areasN;
	float *scatteringN;
    float del;
    KL_JunctionU ** kl_j;
    KL_JunctionU ** kl_jN;
    TUBE ** tubes;
    TUBE ** tubesN;
    NJunctionU<3> NoseJunction;//();//{};
};


extern "C" {

	void HumanVNdelU_next(HumanVNdelU *unit, int inNumSamples);
	void HumanVNdelU_Ctor(HumanVNdelU* unit);
	void HumanVNdelU_Dtor(HumanVNdelU* unit);

}


////output= HumanVNdelU.ar(input, loss,rg,rl,rn,area1len, areas,areasN );
//(input, loss,rg,rl,rn,lmix,nmix,area1len,numtubes, areas,areasNlen,areasN );
#define AREASZINdel 12
void HumanVNdelU_Ctor(HumanVNdelU* unit) {

	int i,j;
    unit->Hgo = unit->Hlo = unit->Hgi = unit->Hli = unit->Hno = unit->Hni=0.0;
	//int numinputs = unit->mNumInputs;
	//int numtubes= numinputs-5;
	//unit->numtubes= numtubes;


    
	unit->loss= ZIN0(1);
    unit->rg= ZIN0(2);
    unit->rl= ZIN0(3);
    unit->rn= ZIN0(4);
    unit->lmix= ZIN0(5);
    unit->nmix= ZIN0(6);
    unit->area1len = ZIN0(7);
    unit->del = ZIN0(8);
    unit->numtubes = ZIN0(11);
    if(unit->numtubes<2) {
		printf("too few tubes! only %d \n", unit->numtubes);
		return;
	}
    
    unit->areas= (float*)RTAlloc(unit->mWorld, (unit->numtubes) * sizeof(float));
	unit->scattering= (float*)RTAlloc(unit->mWorld, (unit->numtubes-2) * sizeof(float));
    
    unit->tubes = (TUBE **)RTAlloc(unit->mWorld, (unit->numtubes)*sizeof(TUBE*));
    for (i=0; i<(unit->numtubes); ++i){
        unit->tubes[i] = new(unit) TUBE(unit);
    }
    
    unit->kl_j = (KL_JunctionU **)RTAlloc(unit->mWorld, (unit->numtubes-2)*sizeof(KL_JunctionU*));
    for (i=0; i<(unit->numtubes-2); ++i){
        unit->kl_j[i] = new(unit) KL_JunctionU(unit);
    }
    
    unit->areasNlen = ZIN0(unit->numtubes + AREASZINdel);
    
    unit->areasN = (float*)RTAlloc(unit->mWorld, (unit->areasNlen) * sizeof(float));
	unit->scatteringN= (float*)RTAlloc(unit->mWorld, (unit->areasNlen-1) * sizeof(float));
    
    unit->tubesN = (TUBE **)RTAlloc(unit->mWorld, (unit->areasNlen)*sizeof(TUBE*));
    for (i=0; i<unit->areasNlen; ++i){
        unit->tubesN[i] = new(unit) TUBE(unit);
    }
    
    unit->kl_jN = (KL_JunctionU **)RTAlloc(unit->mWorld, (unit->areasNlen-1)*sizeof(KL_JunctionU*));
    for (i=0; i<(unit->areasNlen-1); ++i){
        unit->kl_jN[i] = new(unit) KL_JunctionU(unit);
    }
    unit->NoseJunction.Init();
	SETCALC(HumanVNdelU_next);
	(unit->mCalcFunc)(unit, 1);
}

void HumanVNdelU_Dtor(HumanVNdelU* unit) {

	int i;

    for (i=0; i<(unit->numtubes-2); ++i){
        delete unit->kl_j[i];
    }
    for (i=0; i<(unit->areasNlen-1); ++i){
        delete unit->kl_jN[i];
    }
    RTFree(unit->mWorld,unit->kl_j);
    RTFree(unit->mWorld,unit->kl_jN);
    
    for (i=0; i<(unit->numtubes); ++i){
        delete unit->tubes[i];
    }
    for (i=0; i<(unit->areasNlen); ++i){
        delete unit->tubesN[i];
    }
    
    RTFree(unit->mWorld,unit->tubes);
    RTFree(unit->mWorld,unit->tubesN);
    
	RTFree(unit->mWorld, unit->scattering);
    RTFree(unit->mWorld, unit->areas);
    
    RTFree(unit->mWorld, unit->scatteringN);
    RTFree(unit->mWorld, unit->areasN);
}

void Calc_Scattering(HumanVNdelU *unit){
    //float *k = unit->scattering;
    float *areas = unit->areas;
    KL_JunctionU **kl_j = unit->kl_j;
    
    for(int i = 0;i<unit->area1len-1;i++){
        //k[i] = getK(areas[i],areas[i+1]);
        kl_j[i]->set_areas(areas[i],areas[i+1]);
    }
    for(int i=unit->area1len;i<unit->numtubes-1;i++){
        //k[i-1] = getK(areas[i],areas[i+1]);
        kl_j[i-1]->set_areas(areas[i],areas[i+1]);
    }
    
    //float *kN = unit->scatteringN;
    float *areasN = unit->areasN;
    KL_JunctionU **kl_jN = unit->kl_jN;
    for(int i=0;i<unit->areasNlen -1;i++){
        //kN[i] = getK(areasN[i],areasN[i+1]);//(areasN[i] - areasN[i+1])/(areasN[i] + areasN[i+1]);
        kl_jN[i]->set_areas(areasN[i],areasN[i+1]);
    }
    

    float NJareas[3];
    NJareas[0] = areas[unit->area1len-1];
    NJareas[1] = areas[unit->area1len];
    NJareas[2] = areasN[0];
    unit->NoseJunction.calc_alpha(NJareas);
}
void HumanVNdelU_next(HumanVNdelU *unit, int inNumSamples) {

	int i,j;
    KL_JunctionU **kl_j = unit->kl_j;
    KL_JunctionU **kl_jN = unit->kl_jN;
	int numtubes= unit->numtubes;

	//value to store
	float * in= IN(0);
	float * out= OUT(0);


	unit->loss= ZIN0(1);
    unit->rg= ZIN0(2);
    unit->rl= ZIN0(3);
    unit->rn= ZIN0(4);
    unit->lmix= ZIN0(5);
    unit->nmix= ZIN0(6);
    unit->area1len = ZIN0(7);
    unit->del = sc_max(ZIN0(8),1.0);
    
    float * in_noise = IN(9);
    int noiseloc = ZIN0(10);
    
    //get areas
    for(i=0;i<numtubes;i++)
        unit->areas[i] = sc_max(ZIN0(i + AREASZINdel),0.0);
    
    for(i=0;i<unit->areasNlen;i++)
        unit->areasN[i] = sc_max(ZIN0(i + unit->numtubes + AREASZINdel + 1),0.0);
    
    //set tubes
    for (j=0; j<numtubes; ++j) {
            unit->tubes[j]->del = unit->del - 1.0;
            unit->tubes[j]->set_area(unit->areas[j], unit->loss);
	}
    for (j=0; j<(unit->areasNlen); ++j) {
            unit->tubesN[j]->del = unit->del - 1.0;
            unit->tubesN[j]->set_area(unit->areasN[j], unit->loss);
	}
    Calc_Scattering(unit);

    
	for (i=0; i<inNumSamples; ++i) {

		
        //update ins
        unit->Hgi = in[i] + unit->rg * unit->tubes[0]->outL;
        //fill tubes area1
        unit->tubes[0]->inR = unit->Hgo;
        unit->tubes[0]->inL = kl_j[0]->outL;        
        for (j=1; j<(unit->area1len-1); ++j) {
            unit->tubes[j]->inR = kl_j[j-1]->outR;
            unit->tubes[j]->inL = kl_j[j]->outL;
		}
        unit->tubes[unit->area1len-1]->inR = kl_j[unit->area1len-2]->outR;
        unit->tubes[unit->area1len-1]->inL = unit->NoseJunction.out[0];
        //fill tubes area2
        unit->tubes[unit->area1len]->inR = unit->NoseJunction.out[1];
        unit->tubes[unit->area1len]->inL = kl_j[unit->area1len -1]->outL;
        for (j=unit->area1len + 1; j<(unit->numtubes-1); ++j) {
            unit->tubes[j]->inR = kl_j[j-2]->outR;
            unit->tubes[j]->inL = kl_j[j-1]->outL;
		}
        unit->tubes[unit->numtubes-1]->inR = kl_j[unit->numtubes-3]->outR;
        unit->tubes[unit->numtubes-1]->inL = unit->Hlo*unit->rl;
        
        //fill noise
        if(noiseloc > 0 && noiseloc <= numtubes)
            unit->tubes[noiseloc-1]->inR += in_noise[i];
        
        //advance tubes area1 and 2
        for (j=0; j<numtubes; ++j) {
            unit->tubes[j]->go();
		}
        //fill tubesN
        unit->tubesN[0]->inR = unit->NoseJunction.out[2];
        unit->tubesN[0]->inL = kl_jN[0]->outL;
        for (j=1; j<unit->areasNlen-1; ++j) {
            unit->tubesN[j]->inR= kl_jN[j-1]->outR;
            unit->tubesN[j]->inL= kl_jN[j]->outL;
		}
        unit->tubesN[unit->areasNlen-1]->inR = kl_jN[unit->areasNlen-2]->outR;
        unit->tubesN[unit->areasNlen-1]->inL = unit->Hno*unit->rn;
        
        //advance tubesN
        for (j=0; j<(unit->areasNlen); ++j) {
            unit->tubesN[j]->go();
		}
        //fill junctions area1
        for (j=0; j<(unit->area1len-1); ++j) {
            kl_j[j]->inR = unit->tubes[j]->outR;
            kl_j[j]->inL = unit->tubes[j + 1]->outL;
        }
        //fill nose
        unit->NoseJunction.in[0] = unit->tubes[unit->area1len-1]->outR;
        unit->NoseJunction.in[1] = unit->tubes[unit->area1len]->outL;
        unit->NoseJunction.in[2] = unit->tubesN[0]->outL;
        
        //fill junctions area2
        for (j=unit->area1len-1; j<(unit->numtubes-2); ++j) {
            kl_j[j]->inR = unit->tubes[j + 1]->outR;
            kl_j[j]->inL = unit->tubes[j + 2]->outL;
        }
        
        unit->Hli = unit->tubes[numtubes-1]->outR;
        ///////noseareas
        
        //fill junctionsN
        for (j=0; j<(unit->areasNlen-1); ++j) {
            kl_jN[j]->inR = unit->tubesN[j]->outR;
            kl_jN[j]->inL = unit->tubesN[j + 1]->outL;
		}
        
        unit->Hni = unit->tubesN[unit->areasNlen-1]->outR;
        //compute outs

        
        for (j=0; j<(numtubes-2); ++j) {
            kl_j[j]->go();
		}
        for (j=0; j<(unit->areasNlen-1); ++j) {
            kl_jN[j]->go();
		}
        unit->NoseJunction.go();
        
        unit->Hgo = unit->Hgi;
        unit->Hlo = unit->Hli;
        unit->Hno = unit->Hni;
        
        out[i] = unit->Hlo*unit->lmix + unit->Hno*unit->nmix;
	}
}
//////////////////////

struct GaussianNoise{
    int isset;
    float gset;
    RGen* rgen;
    void init(Unit* unit){
        rgen = unit->mParent->mRGen;
        isset = 0;
    }
    float get(){
        float fac,rsq,v1,v2;
        if(isset == 0){
            do{
                v1 = 2.0*rgen->frand() - 1.0;
                v2 = 2.0*rgen->frand() - 1.0;
                rsq = v1*v1 + v2*v2;
            }while(rsq >= 1.0 || rsq == 0.0);
            fac = sqrt(-2.0*log(rsq)/rsq);
            gset = v1*fac;
            isset = 1;
            return v2*fac;
        }else{
            isset = 0;
            return gset;
        }
    }
};
//gives hanning window centered at 0 and 1 width
float hanning(float t){
    if (fabs(t) > 0.5)
        return 0.0;
    return 0.5*(1.0 + cos(2.0*M_PI*t));
}
inline float S1(float t){
    return t - floor(t);
}

float hpB[] = {0.7786600804279, -1.5573201608558, 0.7786600804279};
float hpA[] = { -1.5077255607275, 0.60691476098406};
float deB[] = {1.0};
float deA[] = { -0.9};
struct LFglottal:public Unit
{
	LFglottal(Unit* unit);
	float freq,phase,alpha,Te,Tp,Ta,mFreqMul;
    GaussianNoise GN;
    LTIT<hpB,NumElements(hpB),hpA,NumElements(hpA)> hpfilter;
    LTIT<deB,NumElements(deB),deA,NumElements(deA)> deemph;
};
SCWrapClassINI(LFglottal);
LFglottal::LFglottal(Unit* unit){
		phase = 0.0;
        alpha = ZIN0(4);
        mFreqMul = unit->mRate->mSampleDur;
		freq = ZIN0(0);
        GN.init(unit);
		//SETCALC(LFglottal_next);
}
void LFglottal_next(LFglottal* unit,int inNumSamples){
    float * out = OUT(0);
    unit->freq = ZIN0(0);
    float phaseinc = unit->freq * unit->mFreqMul;
    //float Tp = unit->Tp = sc_min(ZIN0(1),0.8);
    //float Te = unit->Te = sc_clip(ZIN0(2),Tp*1.1,1);
    float Te = unit->Te = sc_clip(ZIN0(2),0.1,1);
    float Tp = unit->Tp = sc_clip(ZIN0(1),Te*0.51,Te*0.9);
    float Ta = unit->Ta = ZIN0(3);
    float t = unit->phase;
    float alpha = unit->alpha = ZIN0(4);
    float hamp = ZIN0(5);
	float hampto1 = 1.0 - hamp;
    float hwidth = sc_clip(ZIN0(6),0.0,1.0);
    //calc Hanning constants
    float hoff = -Te + 0.5;
    double hwinv = 1.0/hwidth;
    //calc LF constants
    float ep = 1/Ta;
    float expfac = exp(-ep*(1-Te));
	float om_g = M_PI/Tp;
    float E0 = -1.0/(exp(alpha*Te)*sin(om_g*Te));
	float Ei = 1/(1.0 - expfac);
    float sig;
    float than;
    float aspir;
    for(int i=0;i<inNumSamples; i++){
        t += phaseinc;
        if (t >= 1.0)
               t -= 1.0;//t = fmod(t,1.0); //0.0;
        if(t < Te){
			sig = E0*exp(alpha*t)*sin(om_g*t);
		}else{
			sig =  -Ei*(exp(-ep*(t-Te)) - expfac);
        }
        //than = (S1(t - Te + 0.5) - 0.5)/hwidth;
        than = (S1(t + hoff) - 0.5)*hwinv;
        aspir = unit->GN.get() *hanning(than)*hamp;
        aspir = unit->hpfilter.filter(aspir);
        aspir = unit->deemph.filter(aspir);
        out[i] = hampto1*sig + aspir*hamp;
    }
    unit->phase = t;
}
//////////////////////
//Chen, G., Shue, Y.-L., Kreiman, J., Alwan, A., 2012. Estimating the voicesource in noise. In: Interspeech.
struct ChenglottalU:public Unit
{
	ChenglottalU(Unit* unit);
	float freq,phase,mFreqMul;
    float OQ,Scp,Sop,asym;
};
SCWrapClassINI(ChenglottalU);
ChenglottalU::ChenglottalU(Unit* unit){
		phase = 0.0;
        mFreqMul = unit->mRate->mSampleDur;
		freq = ZIN0(0);
		//SETCALC(ChenglottalU_next);
}
inline float F(float x,float L,float den){
    float M_PIx = M_PI*x;
    return (exp(L*x)*(L*sin(M_PIx) - M_PI*cos(M_PIx)) + M_PI)/den; //den = (M_PI*(exp(L) + 1.0));
}
void ChenglottalU_next(ChenglottalU* unit,int inNumSamples){
    float * out = OUT(0);
    unit->freq = ZIN0(0);
    float phaseinc = unit->freq * unit->mFreqMul;
    float OQ = unit->OQ = sc_clip(ZIN0(1),0.0,1.0);
    unit->asym = sc_clip(ZIN0(2),0.0,1.0);
    unit->Sop = sc_clip(ZIN0(3),0.0,1.0);
    unit->Scp = sc_clip(ZIN0(4),0.0,1.0);
    float t = unit->phase;
    //calc constants
    float to = OQ * unit->asym;
    float tc = OQ - to;
    float Lop = 12.0*(0.5 - unit->Sop);
    float Lcp = 12.0*(0.5 - unit->Scp);
    float denop = M_PI*(exp(Lop) + 1.0);
    float dencp = M_PI*(exp(Lcp) + 1.0);
    for(int i=0;i<inNumSamples; i++){
        t += phaseinc;
        if (t >= 1.0)
               t = fmod(t,1.0); //0.0;
        if(t < to)
			out[i] = F(t/to,Lop,denop);
		else if(t < OQ) //< to + tc 
			out[i] =  F((OQ - t)/tc,Lcp,dencp);
        else
            out[i] = 0.0;
    }
    unit->phase = t;
}
//Veldhuis
struct VeldhuisGlot:public Unit
{
	VeldhuisGlot(Unit* unit);
    void Calc();
	float freq,phase,te,tp,ta,mFreqMul;
    float A,tx,fff,fte;
    GaussianNoise GN;
    LTIT<hpB,NumElements(hpB),hpA,NumElements(hpA)> hpfilter;
};
SCWrapClassINI(VeldhuisGlot);
VeldhuisGlot::VeldhuisGlot(Unit* unit){
		phase = 0.0;
        mFreqMul = unit->mRate->mSampleDur;
		freq = ZIN0(0);
        te  = sc_clip(ZIN0(2),0.1,1);
        tp = sc_clip(ZIN0(1),te*0.505,te*0.9);
        ta = ZIN0(3);
        Calc();
        GN.init(unit);
		//SETCALC(VeldhuisGlot_next);
}
void VeldhuisGlot::Calc(){
    float tt = (1.0-te)/ta;
	float D = 1.0 - tt/(exp(tt) -1.0 );
    float den = 2.0*te*te-3.0*te*tp+6.0*ta*(te-tp)*D;
	//--if den >= 0 then prerror("den is >=0") end
	float num = (0.5*te*te-te*tp);
	//--if num >= 0 then prerror("num is >=0") end
	float ff = num/den;
	tx = te*(1 - ff);
	//if tx < te then prerror"tx < te" end
	fte = 4.0*te*(tp-te)*(tx-te);
	fff = exp(-(1.0-te)/ta);
	A = -1/fte;
}
void VeldhuisGlot_next(VeldhuisGlot* unit,int inNumSamples){
    float * out = OUT(0);
    unit->freq = ZIN0(0);
    float phaseinc = unit->freq * unit->mFreqMul;
    //float Tp = unit->Tp = sc_min(ZIN0(1),0.8);
    //float Te = unit->Te = sc_clip(ZIN0(2),Tp*1.1,1);
    float Te = sc_clip(ZIN0(2),0.1,1);
    float Tp =  sc_clip(ZIN0(1),Te*0.505,Te*0.9);
    float Ta  = ZIN0(3);
    float t = unit->phase;
    float hamp = ZIN0(4);
    float hwidth = sc_clip(ZIN0(5),0.0,1.0);
    //calc Hanning constants
    float hoff = -Te + 0.5;
    double hwinv = 1.0/hwidth;
    if (Te!=unit->te || Tp!=unit->tp || Ta!=unit->ta){
        unit->te = Te;
        unit->tp = Tp;
        unit->ta = Ta;
        unit->Calc();
    }

    float sig;
    float than;
    float aspir;
    float A = unit->A;
    float tx = unit->tx;
    float fff = unit->fff;
    float fte = unit->fte;
    for(int i=0;i<inNumSamples; i++){
        t += phaseinc;
        if (t >= 1.0)
               t -= 1.0;//t = fmod(t,1.0); //0.0;
        if(t < Te){
			sig = 4.0*A*t*(Tp-t)*(tx-t);
		}else{
			sig =  A*fte*(exp(-(t-Te)/Ta)-fff)/(1 - fff);
        }
        //than = (S1(t - Te + 0.5) - 0.5)/hwidth;
        than = (S1(t + hoff) - 0.5)*hwinv;
        aspir = unit->GN.get() *hanning(than)*hamp;
        aspir = unit->hpfilter.filter(aspir);
        out[i] = sig + aspir;
    }
    unit->phase = t;
}
PluginLoad(SLUGens)
{

	ft = inTable;

    DefineDtorUnit(HumanV);
	DefineDtorUnit(HumanVdel);
    DefineDtorUnit(HumanVN);
    DefineDtorUnit(HumanVNdel);
    DefineDtorUnit(HumanVNdelO2);
    DefineDtorUnit(HumanVNdelU);
    DefineDtorUnit(LFglottal);
    DefineDtorUnit(ChenglottalU);
    DefineDtorUnit(VeldhuisGlot);
}
