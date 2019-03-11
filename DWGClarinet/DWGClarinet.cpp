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
//[,-0.0510756103675,-0.059692489429428,-0.11161650972268,] / [,0.27334492044681,-1.0509645904693,1,]
float reflB[] = {-0.11161650972268,-0.059692489429428,-0.0510756103675};
float reflA[] = {-1.0509645904693,0.27334492044681};
//////////////////////////////////////////////////////////
struct DWGClarinet:public Unit
{
	LagrangeT<MAXDELAY> DWGF;
	FilterC1C3 Loss;
	DWGClarinet(Unit* unit);
	void Release(float trig,float *out,int NumSamples);
	float m_trig;
	int relcount;
	float rellevel;
	float rellevelstep;
};
SCWrapClass(DWGClarinet);
DWGClarinet::DWGClarinet(Unit* unit){
		m_trig = 0.0;
		float release = ZIN0(5);
		relcount = SAMPLERATE * release;
		rellevel = 1.0;
		rellevelstep = 1.0/(float)relcount;
		SETCALC(DWGClarinet_next);
}
struct DWGClarinet2:public DWGClarinet
{
	LTIT<reflB,3,reflA,2> reflfilt;
	DWGClarinet2(Unit* unit);
};
SCWrapClass(DWGClarinet2);
DWGClarinet2::DWGClarinet2(Unit* unit):DWGClarinet(unit){SETCALC(DWGClarinet2_next);};

struct DWGClarinet3:public DWGClarinet
{
	LTIT<reflB,3,reflA,2> reflfilt;
	DWGClarinet3(Unit* unit);
};
SCWrapClass(DWGClarinet3);
DWGClarinet3::DWGClarinet3(Unit* unit):DWGClarinet(unit){SETCALC(DWGClarinet3_next);};

void DWGClarinet::Release(float trig,float *out,int NumSamples){
	
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
float reflection1(float pd,float pc,float m)
{
	if(pd >= pc)
		return 1.;
	//return sc_clip(1 + m*(pd - pc),-1,1);
	return 1 + m*(pd - pc);
}
float reflection2(float pd,float pc,float m)
{
	//if(pd <= pc)
	//	return 1.;
	return sc_clip(pc - m*(pd),-1,1);
}
float Ureed(float pdelta,float Pc){
	float w = 1;
	float y0 = 1;
	float rho = 1;
	float sgn = (pdelta > 0)?1 :-1;
	if (pdelta > Pc) 
		return 0; 
	return w*y0*(1-pdelta/Pc)*sqrt(2*fabs(pdelta)/rho)*sgn;
}
float RefTable(float pdelta,float Pc){
	if(pdelta == 0)
		return -1;
	float aa = Ureed(pdelta,Pc)/pdelta;
	return (1-aa)/(1+aa);
}
float minPm(float m,float Pc){
	float raiz = sqrt(pow(m*Pc,2) + 4.0f);
	return -(raiz -m*Pc -2)/m ;//,(raiz +m*Pc +2)/m
}

/////////////////////////////////////////////////////////

void DWGClarinet_next(DWGClarinet *unit, int inNumSamples)
{

	float *out = OUT(0);
	float freq = ZIN0(0);
	float *Pm = IN(1);
	float Pc = ZIN0(2);
	float m = ZIN0(3);
	float trig = ZIN0(4);

	float c1 = ZIN0(6);
	float c3 = std::max(ZIN0(7),(float)1e-9);

	float minpm = 0;//minPm(m,Pc);
	//Print("minpm %g\n",minpm);
	unit->Loss.setcoeffs(freq,c1,c3);
	float lossdelay = unit->Loss.phasedelay(freq,SAMPLERATE);
	float deltot = 0.5*SAMPLERATE/freq;
	//float del1 = (deltot - lossdelay)*0.5 - 1;
	float del1 = (deltot - lossdelay ) - 1;

	float PMAS,PMAS2;
	float PMENOS;
	for (int i=0; i < inNumSamples; ++i)
	{
		PMENOS = unit->DWGF.delay(del1);
		PMENOS = -unit->Loss.filter(PMENOS);

		float pmdiv2 = (minpm + Pm[i])/2;
		float pdeltap = pmdiv2 - PMENOS;
		PMENOS = -reflection1(pdeltap,Pc,m)*pdeltap + pmdiv2;
		
		//float pdeltap = PMENOS - Pm[i]/2;
		//PMENOS = reflection2(pdeltap,Pc,m)*pdeltap + Pm[i]/2;
		
		//float pdeltap = Pm[i] - (unit->DWGF.get(0) + PMENOS);
		//PMENOS = RefTable(pdeltap*Pc,Pc)*(PMENOS - Pm[i]/2) + Pm[i]/2;

		unit->DWGF.push(PMENOS);

		out[i] = PMENOS ;//PMAS;// + PMAS2;

	}
	unit->Release(trig,out,inNumSamples);
}

void DWGClarinet3_next(DWGClarinet3 *unit, int inNumSamples)
{

	float *out = OUT(0);
	float freq = ZIN0(0);
	float *Pm = IN(1);
	float Pc = ZIN0(2);
	float m = ZIN0(3);
	float trig = ZIN0(4);

	float c1 = ZIN0(6);
	float c3 = std::max(ZIN0(7),(float)1e-9);

	float minpm = 0;//minPm(m,Pc);
	//Print("minpm %g\n",minpm);
	unit->Loss.setcoeffs(freq,c1,c3);
	float lossdelay = unit->Loss.phasedelay(freq,SAMPLERATE);
	
	float omega = 2*M_PI*freq/SAMPLERATE;
	//float refldelay = (unit->reflfilt.phasedelay(freq,SAMPLERATE)*omega + M_PI)/omega;
	float refldelay = (M_PI-unit->reflfilt.phase(freq,SAMPLERATE))/omega;
	
	float deltot = 0.5*SAMPLERATE/freq;
	//float del1 = (deltot - lossdelay)*0.5 - 1;
	float del1 = (deltot - lossdelay -refldelay) - 1;
//Print("freq %g,refldelay %g refldelay2 %g deltot %g\n",freq,refldelay,refldelay2,deltot);
	float PMAS,PMAS2,OUTv;
	float PMENOS;
	for (int i=0; i < inNumSamples; ++i)
	{
		PMENOS = unit->DWGF.delay(del1);
		PMENOS = unit->Loss.filter(PMENOS);
		OUTv = PMENOS;
		PMENOS = unit->reflfilt.filter(PMENOS);
		//OUTv += PMENOS;
		float pmdiv2 = (minpm + Pm[i])/2;
		float pdeltap = pmdiv2 - PMENOS;
		PMENOS = -reflection1(pdeltap,Pc,m)*pdeltap + pmdiv2;
		
		//float pdeltap = PMENOS - Pm[i]/2;
		//PMENOS = reflection2(pdeltap,Pc,m)*pdeltap + Pm[i]/2;
		
		//float pdeltap = Pm[i] - (unit->DWGF.get(0) + PMENOS);
		//PMENOS = RefTable(pdeltap*Pc,Pc)*(PMENOS - Pm[i]/2) + Pm[i]/2;

		unit->DWGF.push(PMENOS);

		out[i] = OUTv;//PMENOS ;//PMAS;// + PMAS2;

	}
	unit->Release(trig,out,inNumSamples);
}
void DWGClarinet2_next(DWGClarinet2 *unit, int inNumSamples)
{

	float *out = OUT(0);
	float freq = ZIN0(0);
	float *Pm = IN(1);
	float Pc = ZIN0(2);
	float m = ZIN0(3);
	float trig = ZIN0(4);

	float c1 = ZIN0(6);
	float c3 = std::max(ZIN0(7),(float)1e-9);

	
	unit->Loss.setcoeffs(freq,c1,c3);
	float lossdelay = unit->Loss.phasedelay(freq,SAMPLERATE);
	float refldelay = unit->reflfilt.groupdelay(freq,SAMPLERATE);
	float deltot = 0.5*SAMPLERATE/freq;
	//float del1 = (deltot - lossdelay)*0.5 - 1;
	float del1 = (deltot - lossdelay - refldelay) - 1;

	float PMAS,PMAS2;
	float PMENOS,OUTv;
	for (int i=0; i < inNumSamples; ++i)
	{
		PMENOS = unit->DWGF.delay(del1);
		//PMENOS = -unit->Loss.filter(PMENOS);
		PMENOS = unit->Loss.filter(PMENOS);
		OUTv = PMENOS;
		PMENOS = unit->reflfilt.filter(PMENOS);
		OUTv += PMENOS;

		float pdeltap = Pm[i] - (unit->DWGF.get(0) + PMENOS);
		//float pdeltap = Pm[i] -  PMENOS;
		PMENOS = RefTable(pdeltap*Pc,Pc)*(PMENOS - Pm[i]/2) + Pm[i]/2;

		unit->DWGF.push(PMENOS);

		out[i] = OUTv ;//PMAS;// + PMAS2;

	}
	unit->Release(trig,out,inNumSamples);
}
/////////////////////////////////////////
PluginLoad(DWGClarinet)
{
	ft = inTable;
	DefineDtorUnit(DWGClarinet);
	DefineDtorUnit(DWGClarinet2);
	DefineDtorUnit(DWGClarinet3);
}
