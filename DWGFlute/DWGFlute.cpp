

/*
 *
 *    Copyright (C) 2013 Victor Bombi
 *
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program; if not, write to the Free Software
 *    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
#include "DWG.cpp"
#define MAXDELAY 1024
InterfaceTable *ft;

float repair_input(float val,float min,float max){
	if( std::fpclassify(val)== FP_INFINITE ||  std::fpclassify(val)== FP_NAN )
		return min;
	else
		return sc_clip(val,min,max);
	
}

//////////////////////////////////////////////////////////
//[,-0.42471514412277,-0.2264060561713,] / [,-0.34889059582329,1,]
float reflB[] = {-0.2264060561713,-0.42471514412277};
float reflA[] = {-0.34889059582329};
struct DWGFlute:public Unit
{
	LagrangeT<MAXDELAY> boreDelay;
	LagrangeT<MAXDELAY> jetDelay;
	LTIT<reflB,2,reflA,1> reflfilt;
	DCBlocker dcblocker;
	DCBlocker dcblockerout;
	DWGFlute(Unit* unit);
	void Release(float trig,float *out,int NumSamples);
	float m_trig;
	int relcount;
	float rellevel;
	float rellevelstep;
};
SCWrapClass(DWGFlute);
DWGFlute::DWGFlute(Unit* unit){
		m_trig = 0.0;
		float release = ZIN0(6);
		relcount = SAMPLERATE * release;
		rellevel = 1.0;
		rellevelstep = 1.0/(float)relcount;
		SETCALC(DWGFlute_next);
}

struct DWGFlute2:public DWGFlute
{
	float lastN;
	DWGFlute2(Unit* unit);
};
SCWrapClass(DWGFlute2);
DWGFlute2::DWGFlute2(Unit* unit):DWGFlute(unit),lastN(0){SETCALC(DWGFlute2_next);};

void DWGFlute::Release(float trig,float *out,int NumSamples){
	
	if(this->m_trig <=0 && trig > 0){
		this->m_trig = trig;
	}
	if((this->m_trig >0 && trig <= 0)){
		int relcount = this->relcount;
		float rellevel = this->rellevel;
		float rellevelstep = this->rellevelstep;
		
		for(int i=0; i<NumSamples; i++){
			if(relcount > 0){
				rellevel -= rellevelstep;
				relcount--;
			}
			out[i] *=rellevel;
		}
		if(relcount <=0)
			DoneAction(2,this);
			
		this->relcount = relcount;
		this->rellevel = rellevel;
	}
}


/////////////////////////////////////////////////////////
inline float JetTable( float input )
{
  // Perform "table lookup" using a polynomial
  // calculation (x^3 - x), which approximates
  // the jet sigmoid behavior.
  float lastFrame = input * (input * input - 1.0);

  // Saturate at +/- 1.0.
  if ( lastFrame > 1.0 ) lastFrame = 1.0;
  if ( lastFrame < -1.0 ) lastFrame = -1.0; 
  return lastFrame;
}
float VortexNoise(float N,float R){
	return -0.6 + 0.1*R -0.6*N +2.0*N*N;
}
float EdgeSimul(float X,float N){

	//return b0 + (b1+c1*N)*X + (b2+c2*N)*X*X + (b3+c3*N)*X*X*X;
	float lastFrame  = (-1.0 -0.015*N)*X + (1.0 + 0.015*N)*X*X*X;
	
	// Saturate at +/- 1.0.
  if ( lastFrame > 1.0 ) lastFrame = 1.0;
  if ( lastFrame < -1.0 ) lastFrame = -1.0; 
  return lastFrame;
}
void DWGFlute_next(DWGFlute *unit, int inNumSamples)
{

	float *out = OUT(0);
	float freq = ZIN0(0)*0.6666666666666666;
	float *Pm = IN(1);
	float endReflection = ZIN0(2);
	float jetReflection = ZIN0(3);
	float jetRatio = ZIN0(4);
	float trig = ZIN0(5);


	float omega = 2*M_PI*freq/SAMPLERATE;
	//float refldelay = (unit->reflfilt.phasedelay(freq,SAMPLERATE)*omega + M_PI)/omega;
	float refdelay = (M_PI-unit->reflfilt.phase(freq,SAMPLERATE))/omega;
	
	float deltot = SAMPLERATE/freq;
	//float del1 = (deltot - lossdelay)*0.5 - 1;
	float del1 = (deltot  - refdelay) - 1;
	float deljet = del1 * jetRatio;
	
	float pressureDiff;
	float PMENOS;
	for (int i=0; i < inNumSamples; ++i)
	{
		PMENOS = unit->boreDelay.delay(del1);
		PMENOS = unit->reflfilt.filter(PMENOS);
		PMENOS = unit->dcblocker.filter(PMENOS);
		
		pressureDiff = Pm[i] - (jetReflection * PMENOS);
		
		unit->jetDelay.push(pressureDiff);
		pressureDiff = unit->jetDelay.delay(deljet);
		pressureDiff = JetTable( pressureDiff ) + (endReflection * PMENOS);
		
		unit->boreDelay.push(pressureDiff);
		

		out[i] = unit->dcblockerout.filter(pressureDiff);//PMAS;// + PMAS2;

	}
	unit->Release(trig,out,inNumSamples);
}

void DWGFlute2_next(DWGFlute2 *unit, int inNumSamples)
{

	float *out = OUT(0);
	float freq = ZIN0(0)*0.6666666666666666;
	float *Pm = IN(1);
	float endReflection = ZIN0(2);
	float jetReflection = ZIN0(3);
	float jetRatio = ZIN0(4);
	float trig = ZIN0(5);


	float omega = 2*M_PI*freq/SAMPLERATE;
	//float refldelay = (unit->reflfilt.phasedelay(freq,SAMPLERATE)*omega + M_PI)/omega;
	float refdelay = (M_PI-unit->reflfilt.phase(freq,SAMPLERATE))/omega;
	
	float deltot = SAMPLERATE/freq;
	//float del1 = (deltot - lossdelay)*0.5 - 1;
	float del1 = (deltot  - refdelay) - 1;
	float deljet = del1 * jetRatio;
	
	float pressureDiff;
	float PMENOS;
	for (int i=0; i < inNumSamples; ++i)
	{
		PMENOS = unit->boreDelay.delay(del1);
		PMENOS = unit->reflfilt.filter(PMENOS);
		PMENOS = unit->dcblocker.filter(PMENOS);
		
		pressureDiff = Pm[i] - (jetReflection * PMENOS);
		
		unit->jetDelay.push(pressureDiff);
		pressureDiff = unit->jetDelay.delay(deljet);
		
		unit->lastN = VortexNoise(unit->lastN,PMENOS);
		pressureDiff = EdgeSimul(pressureDiff,unit->lastN) + (endReflection * PMENOS);
		//pressureDiff = JetTable( pressureDiff ) + (endReflection * PMENOS);
		
		unit->boreDelay.push(pressureDiff);
		

		out[i] = unit->dcblockerout.filter(pressureDiff);//PMAS;// + PMAS2;

	}
	unit->Release(trig,out,inNumSamples);
}

/////////////////////////////////////////
PluginLoad(DWGFlute)
{
	ft = inTable;
	DefineDtorUnit(DWGFlute);
	DefineDtorUnit(DWGFlute2);

}
