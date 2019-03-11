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
//#include "SC_PlugIn.h"
//#define MAXDELAY 1024
InterfaceTable *ft;


////////////////////////////////////////
void makewindow(float *A,int N)
{
int i;
float sumW=0.;
float TWOPI = 2.*M_PI;

 	for(i=0;i<N;i++){
		A[i]=0.54 - 0.46*cos(TWOPI*i/(N-1));
		sumW+=A[i]*A[i];
	}
	sumW/=N;
	sumW=sqrt(sumW);
	for(i=0;i<N;i++)
		A[i]/=sumW;
	 

}
/////////////////////////////////////
void makesynwindow(float *A,int N,int D)
{
int i;
float sumW=0.;
float TWOPI = 2.*M_PI;

 	for(i=0;i<N;i++){
		A[i]=0.54 - 0.46*cos(TWOPI*i/(N-1));
		sumW+=A[i]*A[i];
	}
	sumW/=N;
	sumW=sqrt(sumW);
	for(i=0;i<N;i++)
		A[i]/=sumW;

	sumW=0.;
	for(i=0;i<N;i +=D)
	    sumW +=A[i]*A[i];

	for(i=0;i<N;i++)
	    A[i] /=sumW; 

}

void firstdiff(float *B,int N)
{
	for(int a=N-1;a>0;a--)
		B[a]=(B[a]-B[a-1]);
		
}

void invfirstdiff(float *A,int N)
{

	for(int a=1;a<N;a++)
		A[a]=A[a]+A[a-1];

}
// Ordinary ClearUnitOutputs outputs zero, potentially telling the IFFT (+ PV UGens) to act on buffer zero, so let's skip that:
void LPC_ClearUnitOutputs(Unit *unit, int wrongNumSamples)
{
	ZOUT0(0) = -1;
}
/*
#define SCWrapClass(classname) extern "C" {\
	void classname##_Ctor(classname* unit);\
	void classname##_next(classname *unit, int inNumSamples);\
	void classname##_Dtor(classname *unit);\
}\
void classname##_Ctor(classname* unit){new(unit) classname(unit);}\
void classname##_Dtor(classname* unit){unit->~classname();}
*/
static SndBuf * LPCGetBuffer(Unit * unit, uint32 bufnum, const char * ugenName, int inNumSamples)
{
	Print("getting buffer %d\n", bufnum);
	SndBuf *buf;
	World *world = unit->mWorld;

	if (bufnum >= world->mNumSndBufs) {
		int localBufNum = bufnum - world->mNumSndBufs;
		Graph *parent = unit->mParent;
		if (localBufNum <= parent->localMaxBufNum) {
			buf = parent->mLocalSndBufs + localBufNum;
		} else {
			if (unit->mWorld->mVerbosity > -1)
				Print("%s: invalid buffer number (%d).\n", ugenName, bufnum);
			goto handle_failure;
		}
	} else {
		buf = world->mSndBufs + bufnum;
	}

	if (buf->data == NULL) {
		if (unit->mWorld->mVerbosity > -1)
			Print("%s: uninitialized buffer (%i).\n", ugenName, bufnum);
		goto handle_failure;
	}

	return buf;

handle_failure:
	SETCALC(*LPC_ClearUnitOutputs);
	LPC_ClearUnitOutputs(unit, inNumSamples);
	unit->mDone = true;
	return NULL;
}

#define MAX_POLES 80
float lpcD(float *coef,int M,float *Kaes,float *R)
{
//int i,j,c;
float E,K,temp;
//static float *R,*coeftemp;
	float coeftemp[MAX_POLES + 2];
	
   Kaes[0] = coef[0] = coeftemp[0] = -1.;
   Kaes[1] = coef[1] = coeftemp[1] = K = R[0]!=0.? R[1]/R[0]:0.f;
   E=(1-K*K)*R[0];

   for(int i=2;i<=M;i++){
     	temp=0.;
     	for(int j=1;j<i;j++)
      		temp+=coef[j]*R[i-j];
     	temp=R[i]-temp;
     	if(E!=0.)
      		K=temp/E;
     	else
      		K=0.;
	   
		if(K>1.)
			K=1.;
		else if(K<-1.)
			K=-1.;
		 

     	Kaes[i]=coeftemp[i]=K;

     	for(int c=1;c<i;c++)							  //calculo de alfas
	   		coeftemp[c]=coef[c]-K*coef[i-c];
     	memcpy(coef,coeftemp,(M+1)*sizeof(float));

     	E=(1-K*K)*E;
	}
	return(E);
}

void calcerrorD(float *S,float *error,int N,float *coef,int M)
{
	float tempo;
	int a,b;
	int limite;
	int c;
	for(a=0;a<N;a++){
   		tempo=0.;
		limite=__min(a,M);
   		for(b=1,c=a-1;b<=limite;b++,c--){
     		tempo+=coef[b]*S[c];
   		}
   		error[a]=S[a]-tempo;
	}
}

void calculalfa(float *coef,float *Kaes,int M)
{
	float coeftemp[MAX_POLES + 2];
	
	coef[0]=coeftemp[0]=Kaes[0];
	coef[1]=coeftemp[1]=Kaes[1];

 for(int i=2;i<=M;i++){
     coeftemp[i]=Kaes[i];
     for(int c=1;c<i;c++)							  //calculo de alfas
	   coeftemp[c]=coef[c]-Kaes[i]*coef[i-c];
     memcpy(coef,coeftemp,(i+1)*sizeof(float));
  }

}
void sintesisD(float *O,float *E,int N,float *CO,int M,float G)
{

	float tempo;

 	for(int a=0;a<N;a++){
  		O[a]=E[a]*G;
 	}

	int limite;
	//int c;
	for(int a=1;a<N;a++){
		tempo=0.;
		limite=__min(a,M);
   		for(int c=a-1,b=1;b<=limite;b++,c--){	
    	//	O[a]+=CO[b]*O[a-b];
		//	tempo+=CO[b]*O[a-b];
			tempo+=CO[b]*O[c];
   		}
		O[a]+=tempo;
	}
	   
}

/////////////////////////////////////////////////////////////////////////////////////////////
struct SonLPCBase:public Unit
{
	SndBuf *m_coefsndbuf;
	int m_pos, m_coefbufsize, m_audiosize; // "fullbufsize" includes any zero-padding, "audiosize" does not.
	uint32 m_coefbufnum;
	scfft* m_scfft;
	int m_hopsize, m_shuntsize; // These add up to m_audiosize

	int m_numSamples;
	SonLPCBase(Unit *unit);
	World *mWorld;
};

SCWrapClass(SonLPCBase);

SonLPCBase::SonLPCBase(Unit *unit)
{
	mWorld = unit->mWorld;
	m_coefbufnum = (uint32)ZIN0(0);
	m_coefsndbuf = LPCGetBuffer(unit, m_coefbufnum, "SonLPCBase", 1);
	if(!m_coefsndbuf)
		return;
	//m_bufsize = m_coefsndbuf->samples;
	//m_audiosize = m_fullbufsize - (MAX_POLES +1);
	//m_log2n_full  = LOG2CEIL(m_fullbufsize);
	m_pos = 0;
	ZOUT0(0) = ZIN0(0);
	return ;
}


//////////////////////////////////////////////////////////
struct SonLPC:public Unit
{
	SndBuf *m_coefsndbuf;
	int m_pos, m_coefbufsize; // "fullbufsize" includes any zero-padding, "audiosize" does not.
	uint32 m_coefbufnum;
	uint32 m_maxpoles,m_poles;
	int m_hopsize, m_shuntsize; // These add up to m_audiosize

	int m_numSamples;
	World *mWorld;
	
	SonLPC(Unit* unit);
	~SonLPC();
	float * m_inbuf;
	//int m_inbufpos;
	int m_audiosize;
	
	scfft* m_scfft;
	scfft* m_iscfft;
	float m_coef[MAX_POLES + 2];
	float m_kaes[MAX_POLES + 2];
	float *m_fftbuf;
	int m_fftbufsize;
};

SCWrapClass(SonLPC);

SonLPC::SonLPC(Unit* unit):m_inbuf(NULL),m_pos(0){

	mWorld = unit->mWorld;
	m_coefbufnum = (uint32)ZIN0(0);
	m_coefsndbuf = LPCGetBuffer(unit, m_coefbufnum, "SonLPCBase", 1);
	if(!m_coefsndbuf)
		return;
	m_maxpoles = sc_min(MAX_POLES,m_coefsndbuf->samples);
	int poles = (uint32)sc_min(m_maxpoles,sc_max(4, ZIN0(3)));
	m_pos = 0;
	
	m_audiosize = (int)ZIN0(1);
	m_inbuf = (float *)RTAlloc(unit->mWorld,m_audiosize * sizeof(float));
	memset(m_inbuf, 0, m_audiosize * sizeof(float));
	
	m_fftbufsize = NEXTPOWEROFTWO(m_audiosize + m_maxpoles + 2);
	m_fftbuf = (float *)RTAlloc(unit->mWorld,m_fftbufsize * sizeof(float));

	
	int hopsize = (int)(sc_max(sc_min(ZIN0(2), 1.f), 0.f) * m_audiosize);
	if (hopsize < unit->mWorld->mFullRate.mBufLength) {
		Print("SonLPC_Ctor: hopsize smaller than SC's block size (%i) - automatically corrected.\n", hopsize, unit->mWorld->mFullRate.mBufLength);
		hopsize = unit->mWorld->mFullRate.mBufLength;
	} else if (((int)(hopsize / unit->mWorld->mFullRate.mBufLength)) * unit->mWorld->mFullRate.mBufLength
				!= hopsize) {
		Print("SonLPC_Ctor: hopsize (%i) not an exact multiple of SC's block size (%i) - automatically corrected.\n", hopsize, unit->mWorld->mFullRate.mBufLength);
		hopsize = ((int)(hopsize / unit->mWorld->mFullRate.mBufLength)) * unit->mWorld->mFullRate.mBufLength;
	}
	
	m_hopsize = hopsize;
	m_shuntsize = m_audiosize - m_hopsize;
	Print("anal hopsize %d\n",hopsize);
	if (INRATE(4) == calc_FullRate) {
		m_numSamples = unit->mWorld->mFullRate.mBufLength;
	} else {
		m_numSamples = 1;
	}
	
	SCWorld_Allocator alloc(ft, unit->mWorld);
	m_scfft = scfft_create(m_fftbufsize, m_audiosize,(SCFFT_WindowFunction) -1, m_inbuf,m_fftbuf, kForward, alloc);
	m_iscfft = scfft_create(m_fftbufsize, m_audiosize, (SCFFT_WindowFunction) -1, m_fftbuf, m_fftbuf, kBackward, alloc);
	SETCALC(SonLPC_next);
	
	ZOUT0(0) = ZIN0(0);
	OUT(0)[1] = 1;
	OUT(0)[2] = poles;
	return ;
}

SonLPC::~SonLPC(){
	//SCWorld_Allocator alloc(ft, unit->mWorld);
	//if(unit->m_scfft)
	//	scfft_destroy(unit->m_scfft, alloc);
	if(m_inbuf)
		RTFree(mWorld, m_inbuf);
	if(m_fftbuf)
		RTFree(mWorld, m_fftbuf);

}

/////////////////////////////////////////////////////////

void SonLPC_next(SonLPC *unit, int wronginNumSamples)
{
	//float *out = OUT(0);
	float *in = IN(4);
	int poles = (uint32)sc_min(unit->m_maxpoles,sc_max(4, ZIN0(3)));
	int numSamples = unit->m_numSamples;


	float *out = unit->m_inbuf + unit->m_pos + unit->m_shuntsize;
	// copy input
	memcpy(out, in, numSamples * sizeof(float));
	unit->m_pos += numSamples;
	
	if (unit->m_pos != unit->m_hopsize){// || !unit->m_fftsndbuf->data || unit->m_fftsndbuf->samples != unit->m_fullbufsize) {
		if(unit->m_pos == unit->m_hopsize)
			unit->m_pos = 0;
		//Print("SonLPC m_pos %u\n",unit->m_pos);
		ZOUT0(0) = -1.f;
	} else {

		unit->m_pos = 0;
		//SndBuf * bbb = unit->m_fftsndbuf;
		//LOCK_SNDBUF(bbb); 

		//perform autocorrelation
		scfft_dofft(unit->m_scfft);
		
		int numbins = (unit->m_fftbufsize - 2) >> 1;
		float * data = unit->m_fftbuf;
		data[0] = data[0]*data[0];
		data[1] = data[1]*data[1];
		for(int i=1 ; i < numbins ; i++) {
			int re = i*2;
			int im = re +1;
			data[re] = data[re]*data[re] + data[im]*data[im];
			data[im] = 0.0f;
		}
		scfft_doifft(unit->m_iscfft);
		
		float R[MAX_POLES + 2];
	
		float reescalar = float(unit->m_fftbufsize)/float(unit->m_audiosize);

		for(int a=0;a<=poles + 1;a++)
			R[a]=data[a]*reescalar;
			
		float G2 = lpcD(unit->m_coef,poles,unit->m_coefsndbuf->data,R);
		G2=sqrt(G2);
		
		/////////////////////////////////////////////////////
		memmove(unit->m_inbuf, unit->m_inbuf + unit->m_hopsize, unit->m_shuntsize * sizeof(float));
		//Print("SonLPC buffnum %u\n",unit->m_coefbufnum);
		ZOUT0(0) = unit->m_coefbufnum;
		OUT(0)[1] = G2;
		OUT(0)[2] = poles;
	}
}

//#include "../dwgugens/dwglib/dwg.cpp"

////////////////////////////////////////////////////
struct SonLPCSynth : public Unit
{
	float m_coef[MAX_POLES + 2];
	float m_kaes[MAX_POLES + 2];
	uint32 m_poles;
	float m_G;
	SndBuf * m_coefsndbuf;
	uint32 m_coefbufnum;
	LTIv *filter;
	SonLPCSynth(Unit* unit);
};

SCWrapClassT(SonLPCSynth);

SonLPCSynth::SonLPCSynth(Unit *unit)
{
	m_coefbufnum = (uint32)ZIN0(0);
	float *out = OUT(0);
	float G = IN(0)[1];
	float poles = IN(0)[2];
	
	m_G = G;
	m_poles = poles;
	
	m_coefsndbuf = LPCGetBuffer(unit, m_coefbufnum, "SonLPCSynth", 1);
	if(!m_coefsndbuf)
		return;
	//m_maxpoles = sc_min(MAX_POLES,m_coefsndbuf->samples);
	Print("poles %g data %g\n",poles,(float)m_coefsndbuf->samples);
	calculalfa(m_coef,m_coefsndbuf->data,poles);
	
	int nB = 1;//(int)IN0(1);
	int nA = poles ;//(int)IN0(2);
	filter = new(unit) LTIv(unit,nB,nA);
	for(int i=0;i<nB;i++)
		this->filter->KernelB[i] = 0;//G;//IN0(3+i);
	for(int i=0;i<nA;i++)
		this->filter->KernelA[i] = 0;//-m_coef[i + 1];//IN0(3+nB+i);
	SETCALC(SonLPCSynth_next);
}


void SonLPCSynth_next(SonLPCSynth *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(1);
	
	float fbufnum = ZIN0(0);
	float G = IN(0)[1];
	float poles = IN(0)[2];
	if (fbufnum >= 0.f){
		unit->m_poles = poles;
		unit->m_G = G;
		//Print("fbufnum %f\n",fbufnum);
		//for(int i=0;i<=unit->m_poles;i++){
		//	Print("%d es %f\n",i,unit->m_coefsndbuf->data[i]);
		//}
		calculalfa(unit->m_coef,unit->m_coefsndbuf->data,unit->m_poles);
	//coef changes
		int nB = unit->filter->kernel_sizeB;
		int nA = unit->filter->kernel_sizeA;
		for(int i=0;i<nB;i++)
			unit->filter->KernelB[i] = G;//IN0(3+i);
		for(int i=0;i<nA;i++)
			unit->filter->KernelA[i] = -unit->m_coef[i + 1];//IN0(3+nB+i);
	}
	
	for (int i=0; i < inNumSamples; ++i)
	{
		out[i] = unit->filter->filter(in[i]);
	}
	//sintesisD(out,in,inNumSamples,unit->m_coef,unit->m_poles,unit->m_G);

}



/////////////////////////////////////////
PluginLoad(SonLPC)
{
	ft = inTable;
	DefineDtorUnit(SonLPC);
	DefineDtorUnit(SonLPCSynth);
}
