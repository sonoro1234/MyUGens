////////////////////////////////////////////////////////////////////////////////////////////////
//Based on 
//THE PLUCKSYNTH TOUCH STRING
//Fredrik Eckerholm, Gianpaolo Evangelista
//Proc. of the 11th Int. Conference on Digital Audio Effects (DAFx-08), Espoo, Finland, September 1-4, 2008
///////////////////////////////////////////////////////////////////////////////////////////////

#include "SC_PlugIn.h"
#include <float.h>
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define MAXDELAY 1024
static InterfaceTable *ft;
#include "DWG.cpp"

float repair_input(float val,float min,float max){
	if( std::fpclassify(val)== FP_INFINITE ||  std::fpclassify(val)== FP_NAN )
		return min;
	else
		return sc_clip(val,min,max);
	
}
//////////////plucksccatering
struct PluckScattering
{
	LTITv<1,2> Am1;
	//float L,r,rho;
	//float mM,mK,mR;
	float X,rhoL;
	float T;
	PluckScattering(){

	}
	void setParams(float f,float fs,float M,float K,float Rc,float L,float r,float rho){
		//L = 0.65;
		//r = 0.001;
		//rho = 7850;
		rhoL = M_PI*r*r*rho;
		T = (2*L*f)*(2*L*f)*rhoL;
		//float c = sqrt(T*rhoL)
		float c = 2*L*f;
			//Z = sqrt(T*rhoL);
		X = c/fs;
		float R = Rc*T/c;
		float pho = R/(2*sqrt(rhoL*T));
		float kappa = K*X/T;
		float p1 = 1 + M/(rhoL*X) + pho;
		float p0 = 2*M/(rhoL*X) - kappa;
		float pm1 = 1 + M/(rhoL*X) - pho;
		
		float Afil[2];
		Afil[0] = -p0/p1;
		Afil[1] = -(1 - pm1)/p1;
		float Numfil[1] = {1/p1};
		Am1.setKernel(Numfil,Afil);
	}
	void dotick(float F0,float inM,float inm,float &outM,float &outm){
		float prefil = inM + inm + F0 * X / T;
		float postfil = Am1.pushConvol(prefil);
		float other = inm - inM;
		outM = 0.5*(postfil - other);
		outm = 0.5*(postfil + other);
	}
};
////////////////////////////////////////////////////
float coefDCB[2] = {1,-1};
float coefDCA[1] = {-0.995};
struct PluckSynth : public Unit
{

	LagrangeT<MAXDELAY> DWGF[2];
	LTITv<1,1> Loss;
	PluckScattering pscat;
	//LTIT<coefDCB,NumElements(coefDCB),coefDCA,NumElements(coefDCA)> DCblocker;
	PluckSynth(Unit* unit){//Print("PluckSynth constructor\n");
	}
	~PluckSynth(){//Print("PluckSynth destructor\n");
	}
	float m_freq;
	float m_amp;
	float m_trig;
	float DCdelay;
	int relcount;
	float rellevel;
	float rellevelstep;
};

extern "C"
{
	void PluckSynth_Ctor(PluckSynth* unit);
	void PluckSynth_next(PluckSynth *unit, int inNumSamples);
	void PluckSynth_Dtor(PluckSynth *unit);
}

void getcoeffs(float freq,float c1,float c3,float BB[1],float AA[1])
{
	float g = 1.0 - c1/freq; 
	float b = 4.0*c3+freq;
	float a1 = (-b+sqrt(b*b-16.0*c3*c3))/(4.0*c3);
	BB[0] = g*(1+a1);
	AA[0] = a1;
}
//////////////////////////////////////////
void PluckSynth_Ctor(PluckSynth* unit)
{

	new(unit) PluckSynth(unit);

	unit->m_freq = repair_input(ZIN0(0),1.f,20000);
	
	
	unit->m_amp = ZIN0(1);
	unit->m_trig = 0.0;
	
	/*
	float c1 = ZIN0(4);
	float c3 = std::max(ZIN0(5),(float)1e-9);
	float BB[1];
	float AA[1];
	getcoeffs(unit->m_freq,c1,c3,BB,AA);
	unit->Loss.setKernel(BB,AA);
	float lossdelay = unit->Loss.groupdelay(unit->m_freq,SAMPLERATE);
	
	float deltot = SAMPLERATE/unit->m_freq;
	float del1 = deltot*0.5 - lossdelay;
	float del2 = deltot*0.5 ;
	
	//del1 = (int) (del1 + 0.5);
	//del2 = (int) (del2 + 0.5);
	//unit->DCdelay = unit->DCblocker.groupdelay(unit->m_freq,SAMPLERATE);
	*/
	float release = ZIN0(6);
	unit->relcount = SAMPLERATE * release;
	unit->rellevel = 1.0;
	unit->rellevelstep = 1.0/(float)unit->relcount;
	
	//Print("freq %g BB %g AA %g lossdelay %g\n",unit->m_freq,BB[0],AA[0],lossdelay);
	//Print("deltot %g del1, %g del2 %g\n",deltot,del1,del2);
    SETCALC(PluckSynth_next);

}

void PluckSynth_Dtor(PluckSynth* unit)
{
	unit->~PluckSynth ();
}
#define PRINT(A) Print(#A ":%f\n",A);
#define DUMPONNAN(val) \
do{ if( (std::fpclassify(val)== FP_INFINITE ||  std::fpclassify(val)== FP_NAN ) && unit->doprint){ \
PRINT(val) PRINT(P0) PRINT(Flip) PRINT(Slip) PRINT(Uacoust) PRINT(factAreas) PRINT(unit->m_flip)\
PRINT(P1) PRINT(Ulip) PRINT(m) PRINT(outgi) PRINT(convol) PRINT(gate)\
Print("X.x :%f\n",X.x); Print("X.y :%f\n",X.y); \
if(abierto) Print("abierto\n"); else Print("cerrado\n");\
unit->doprint = false; \
}}while(0)
void PluckSynth_next(PluckSynth *unit, int inNumSamples)
{

	float *out = OUT(0);
	float freq = ZIN0(0);
	float amp = ZIN0(1);
	float trig = ZIN0(2);
	float pos = ZIN0(3);

	float c1 = ZIN0(4);
	float c3 = std::max(ZIN0(5),(float)1e-9);
	float *F0 = IN(7);
	float *M = IN(8);
	float *K = IN(9);
	float *R = IN(10);
	
	float L = ZIN0(11);
	float r = ZIN0(12);
	float rho = ZIN0(13);
	
	float BB[1];
	float AA[1];
	getcoeffs(freq,c1,c3,BB,AA);
	unit->Loss.setKernel(BB,AA);
	float lossdelay = unit->Loss.groupdelay(freq,SAMPLERATE);
	float deltot = SAMPLERATE/freq;
	float del1 = (deltot - lossdelay )*0.5 - 1;


	if(unit->m_trig <=0 && trig > 0){
		unit->m_trig = trig;
	}

	float PMAS,PMAS2;
	float PMENOS,OUT1,OUT2;
	for (int i=0; i < inNumSamples; ++i)
	{
		unit->pscat.setParams(freq,SAMPLERATE,M[i],K[i],R[i],L,r,rho);
		float inM = unit->DWGF[0].get(pos*del1);
		float inm = unit->DWGF[1].get(del1*(1-pos));
		float outM,outm;
		unit->pscat.dotick(F0[i],inM,inm,outM,outm);
		
		unit->DWGF[0].set(outM,pos*del1);
		unit->DWGF[1].set(outm,del1*(1-pos));
		
		PMAS = unit->DWGF[0].delay(del1);
		PMAS2 = unit->Loss.pushConvol(PMAS);
		PMENOS = unit->DWGF[1].delay(del1);
		
		unit->DWGF[1].push(-PMAS2);
		unit->DWGF[0].push(-PMENOS);

		out[i] =  PMAS + PMAS2;
		
	}
	if((unit->m_trig >0 && trig <= 0)){
		int relcount = unit->relcount;
		float rellevel = unit->rellevel;
		float rellevelstep = unit->rellevelstep;
		
		for(int i=0; i<inNumSamples; i++){
			if(relcount > 0){
				rellevel -= rellevelstep;
				relcount--;
			}
			out[i] *=rellevel;
		}
		if(relcount <=0)
			DoneAction(2,unit);
			
		unit->relcount = relcount;
		unit->rellevel = rellevel;
	}
	//unit->m_flip = Flip;
	//unit->m_p0 = P0;
	//unit->m_delay = delay;
	//assert(Uacoust == unit->Uacoust);
}

/////////////////////////////////////////
PluginLoad(PluckSynth)
{
	ft = inTable;
	DefineDtorUnit(PluckSynth);
}
