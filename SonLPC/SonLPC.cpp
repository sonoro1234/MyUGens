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
		//A[i]=0.54 - 0.46*cos(TWOPI*i/(N-1));
		A[i]=0.5 - 0.5*cos(TWOPI*i/(N));
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

	makewindow(A,N);
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
#define SCWrapClass(classname) extern "C" {\
	void classname##_Ctor(classname* unit);\
	void classname##_next(classname *unit, int inNumSamples);\
	void classname##_Dtor(classname *unit);\
}\
void classname##_Ctor(classname* unit){new(unit) classname(unit);}\
void classname##_Dtor(classname* unit){unit->~classname();}
static SndBuf * LPCGetBuffer(Unit * unit, uint32 bufnum, const char * ugenName, int inNumSamples)
{
	//Print("getting buffer %d\n", bufnum);
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
	
   Kaes[0]=coef[0]=coeftemp[0]=-1.;
   Kaes[1]=coef[1]=coeftemp[1]=K=R[0]!=0.? R[1]/R[0]:0.f;
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
		limite=sc_min(a,M);
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
		limite=sc_min(a,M);
   		for(int c=a-1,b=1;b<=limite;b++,c--){	
    	//	O[a]+=CO[b]*O[a-b];
		//	tempo+=CO[b]*O[a-b];
			tempo+=CO[b]*O[c];
   		}
		O[a]+=tempo;
	}
	   
}
struct SonLPCBase:public Unit
{
	SndBuf *m_fftsndbuf;
	int m_pos, m_fullbufsize, m_audiosize; // "fullbufsize" includes any zero-padding, "audiosize" does not.
	uint32 m_fftbufnum;
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
	m_fftbufnum = (uint32)ZIN0(0);
	m_fftsndbuf = LPCGetBuffer(unit, m_fftbufnum, "SonLPCBase", 1);
	if(!m_fftsndbuf)
		return;
	m_fullbufsize = m_fftsndbuf->samples;
	//space for poles and making hopsize a multiple of 64
	m_audiosize = m_fullbufsize - 64*4;//(MAX_POLES +1);
	//m_log2n_full  = LOG2CEIL(m_fullbufsize);
	m_pos = 0;
	ZOUT0(0) = ZIN0(0);
	return ;
}


//////////////////////////////////////////////////////////
struct SonLPC:public SonLPCBase
{
	SonLPC(Unit* unit);
	~SonLPC();
	float * m_inbuf;
	int m_inbuffpos;
	float *m_analwindow;
	scfft* m_scfft;
	scfft* m_iscfft;
	float m_coef[MAX_POLES + 2];
	float m_kaes[MAX_POLES + 2];
	float *m_buffer;
	float *m_error;
};

SCWrapClass(SonLPC);

SonLPC::SonLPC(Unit* unit):SonLPCBase(unit),m_inbuf(NULL),m_inbuffpos(0){

	int poles = (uint32)sc_min(MAX_POLES,sc_max(4, ZIN0(3)));
	m_analwindow = (float *)RTAlloc(unit->mWorld,m_audiosize * sizeof(float));
	makewindow(m_analwindow,m_audiosize);
	m_inbuf = (float *)RTAlloc(unit->mWorld,m_audiosize * sizeof(float));
	memset(m_inbuf, 0, m_audiosize * sizeof(float));
	m_buffer = (float *)RTAlloc(unit->mWorld,m_audiosize * sizeof(float));
	m_error = (float *)RTAlloc(unit->mWorld,m_audiosize * sizeof(float));
	
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
	//Print("anal hopsize %d\n",hopsize);
	if (INRATE(1) == calc_FullRate) {
		m_numSamples = unit->mWorld->mFullRate.mBufLength;
	} else {
		m_numSamples = 1;
	}
	
	SCWorld_Allocator alloc(ft, unit->mWorld);
	m_scfft = scfft_create(m_fullbufsize, m_fullbufsize,(SCFFT_WindowFunction) -1, m_fftsndbuf->data,m_fftsndbuf->data, kForward, alloc);
	m_iscfft = scfft_create(m_fullbufsize, m_fullbufsize, (SCFFT_WindowFunction) -1, m_fftsndbuf->data,m_fftsndbuf->data, kBackward, alloc);
	SETCALC(SonLPC_next);
	ZOUT0(0) = ZIN0(0);
	OUT(0)[1] = 1;
	OUT(0)[2] = poles;
}

SonLPC::~SonLPC(){
	//SCWorld_Allocator alloc(ft, unit->mWorld);
	//if(unit->m_scfft)
	//	scfft_destroy(unit->m_scfft, alloc);
	if(m_inbuf)
		RTFree(mWorld, m_inbuf);
	if(m_analwindow)
		RTFree(mWorld, m_analwindow);
	if(m_buffer)
		RTFree(mWorld, m_buffer);
	if(m_error)
		RTFree(mWorld, m_error);
}

/////////////////////////////////////////////////////////

void SonLPC_next(SonLPC *unit, int wronginNumSamples)
{
	//float *out = OUT(0);
	float *in = IN(1);
	int poles = (uint32)sc_min(MAX_POLES,sc_max(4, ZIN0(3)));
	int numSamples = unit->m_numSamples;

	/*
	// copy input incircular buffer
	// for pos to end
	float toend = sc_min(numSamples,m_audiosize - m_inbuffpos);
	memcpy(m_inbuff + m_inbuffpos, in, toend * sizeof(float));
	//for begin to pos
	float rest = numSamples - toend;
	memcpy(m_inbuff, in + toend, rest * sizeof(float));
	m_inbuffpos +=  toend;
	m_inbuffpos = m_inbufposs%m_audiosize;
	*/
	float *out = unit->m_inbuf + unit->m_pos + unit->m_shuntsize;
	// copy input
	memcpy(out, in, numSamples * sizeof(float));
	
	
	unit->m_pos += numSamples;
	
	if (unit->m_pos != unit->m_hopsize || !unit->m_fftsndbuf->data || unit->m_fftsndbuf->samples != unit->m_fullbufsize) {
		if(unit->m_pos == unit->m_hopsize)
			unit->m_pos = 0;
		ZOUT0(0) = -1.f;
		OUT(0)[1] = 1;
		OUT(0)[2] = poles;
	} else {

		unit->m_pos = 0;
		SndBuf * bbb = unit->m_fftsndbuf;
		LOCK_SNDBUF(bbb); 
		
		//copy to sndbuf
		//memcpy(unit->m_fftsndbuf->data, unit->m_inbuf, unit->m_fullbufsize * sizeof(float));
		for(int i = 0; i < unit->m_audiosize ; i++)
			unit->m_fftsndbuf->data[i] = unit->m_inbuf[i] * unit->m_analwindow[i];
			
		firstdiff(unit->m_fftsndbuf->data, unit->m_audiosize);
		
		for(int i = 0; i < unit->m_audiosize ; i++)
			unit->m_buffer[i] = unit->m_fftsndbuf->data[i];
			
		for(int i = unit->m_audiosize; i < unit->m_fullbufsize ; i++)
			unit->m_fftsndbuf->data[i] = 0.0f;
		
		
		//perform autocorrelation
		scfft_dofft(unit->m_scfft);
		
		int numbins = (unit->m_fftsndbuf->samples - 2) >> 1;
		float * data = unit->m_fftsndbuf->data;
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
	
		float reescalar = float(unit->m_fullbufsize)/float(unit->m_audiosize);

		for(int a=0;a<poles + 1;a++)
			R[a]=data[a]*reescalar;
			
		float G2 = lpcD(unit->m_coef,poles,unit->m_kaes,R);
		G2=sqrt(G2);
		
		calcerrorD(unit->m_buffer,unit->m_error,unit->m_audiosize,unit->m_coef,poles);
		for(int i = 0; i < unit->m_audiosize ; i++)
			unit->m_fftsndbuf->data[i] = unit->m_error[i];
		for(int i = unit->m_audiosize,j = 0; i < unit->m_fullbufsize ; i++,j++)
			unit->m_fftsndbuf->data[i] = unit->m_kaes[j];
		/////////////////////////////////////////////////////
		memmove(unit->m_inbuf, unit->m_inbuf + unit->m_hopsize, unit->m_shuntsize * sizeof(float));
		ZOUT0(0) = unit->m_fftbufnum;
		OUT(0)[1] = G2;
		OUT(0)[2] = poles;
	}
}
////////////////////////////////////////////////////////////////////////////////////////////
struct SonLPCSynth:public SonLPCBase
{
	SonLPCSynth(Unit* unit);
	~SonLPCSynth();
	float * m_olabuf;
	float * m_buffer;
	int m_pos;
	float * m_synwindow;
	bool m_windowinited;
};
SCWrapClass(SonLPCSynth);
SonLPCSynth::SonLPCSynth(Unit* unit):SonLPCBase(unit),m_pos(0),m_olabuf(NULL)
{
	m_synwindow = (float*)RTAlloc(unit->mWorld, m_audiosize * sizeof(float));
	m_windowinited = false;
	m_olabuf = (float*)RTAlloc(unit->mWorld, m_audiosize * sizeof(float));
	memset(m_olabuf, 0, m_audiosize * sizeof(float));
	m_buffer = (float*)RTAlloc(unit->mWorld, m_audiosize * sizeof(float));
	
	if (INRATE(1) == calc_FullRate) {
		m_numSamples = unit->mWorld->mFullRate.mBufLength;
	} else {
		m_numSamples = 1;
	}
	
	SETCALC(SonLPCSynth_next);
	//Print("done LPCSynth constructor\n");
}
SonLPCSynth::~SonLPCSynth(){
	//SCWorld_Allocator alloc(ft, unit->mWorld);
	//if(unit->m_scfft)
	//	scfft_destroy(unit->m_scfft, alloc);
	if(m_buffer)
		RTFree(mWorld, m_buffer);
	if(m_olabuf)
		RTFree(mWorld, m_olabuf);
	if(m_synwindow)
		RTFree(mWorld, m_synwindow);
}

void SonLPCSynth_next(SonLPCSynth *unit, int wronginNumSamples)
{
	
	float *out = OUT(0);
	float fbufnum = ZIN0(0);
	float G = IN(0)[1];
	int poles = IN(0)[2];
	int pos     = unit->m_pos;
	int audiosize = unit->m_audiosize;
	int numSamples = unit->m_numSamples;
	
	// Only run the IFFT if we're receiving a new block of input data - otherwise just output data already received
	if (fbufnum >= 0.f){
		//Print("poles %d G: %f buffnu %f\n",poles,G,fbufnum);
		float* error = unit->m_fftsndbuf->data;
		float* kaes = unit->m_fftsndbuf->data + audiosize;
		float coef[MAX_POLES + 2];
		calculalfa(coef, kaes, poles);
		G = 1;
		sintesisD(unit->m_buffer, error, audiosize, coef, poles, G);
		// Then shunt the "old" time-domain output down by one hop
		
		invfirstdiff(unit->m_buffer, audiosize);
		
		int hopsamps = pos;// + unit->mWorld->mFullRate.mBufLength;
		if(!unit->m_windowinited){
			makesynwindow(unit->m_synwindow, audiosize, hopsamps + unit->mWorld->mFullRate.mBufLength);
			unit->m_windowinited = true;
			Print("hopsamps is %d\n",hopsamps + unit->mWorld->mFullRate.mBufLength);
		}
		int shuntsamps = audiosize - hopsamps;
		if(hopsamps != audiosize)  // There's only copying to be done if the position isn't all the way to the end of the buffer
			memmove(unit->m_olabuf, unit->m_olabuf+hopsamps, shuntsamps * sizeof(float));

		// Then mix the "new" time-domain data in - adding at first, then just setting (copying) where the "old" is supposed to be zero.
			// NB we re-use the "pos" variable temporarily here for write rather than read
			for(pos = 0; pos < shuntsamps; ++pos){
				unit->m_olabuf[pos] += unit->m_buffer[pos] * unit->m_synwindow[pos];
			}

		//memcpy(unit->m_olabuf + shuntsamps, fftbuf + shuntsamps, (hopsamps) * sizeof(float));
		for(pos = shuntsamps; pos < audiosize; ++pos){
				unit->m_olabuf[pos] = unit->m_buffer[pos] * unit->m_synwindow[pos];
		}
		// Move the pointer back to zero, which is where playback will next begin
		pos = 0;

	} // End of has-the-chain-fired

	// Now we can output some stuff, as long as there is still data waiting to be output.
	// If there is NOT data waiting to be output, we output zero. (Either irregular/negative-overlap
	//     FFT firing, or FFT has given up, or at very start of execution.)
	if(pos >= audiosize)
		ClearUnitOutputs(unit, numSamples);
	else {
		memcpy(out, unit->m_olabuf + pos, numSamples * sizeof(float));
		pos += numSamples;
	}
	unit->m_pos = pos;
}
//////////////////////////////////////
struct SonLPCSynthIn:public SonLPCBase
{
	SonLPCSynthIn(Unit* unit);
	~SonLPCSynthIn();
	float * m_olabuf;
	float * m_buffer;
	int m_pos;
	float * m_synwindow;
	float *m_analwindow;
	bool m_windowinited;
	float * m_inbuffer;
	int m_inpos;
};
SCWrapClass(SonLPCSynthIn);
SonLPCSynthIn::SonLPCSynthIn(Unit* unit):SonLPCBase(unit),m_pos(0),m_olabuf(NULL),m_inpos(0)
{
	m_analwindow = (float *)RTAlloc(unit->mWorld,m_audiosize * sizeof(float));
	makewindow(m_analwindow,m_audiosize);
	m_synwindow = (float*)RTAlloc(unit->mWorld, m_audiosize * sizeof(float));
	m_windowinited = false;
	m_olabuf = (float*)RTAlloc(unit->mWorld, m_audiosize * sizeof(float));
	memset(m_olabuf, 0, m_audiosize * sizeof(float));
	m_buffer = (float*)RTAlloc(unit->mWorld, m_audiosize * sizeof(float));
	memset(m_buffer, 0, m_audiosize * sizeof(float));
	m_inbuffer = (float*)RTAlloc(unit->mWorld, m_audiosize * sizeof(float));
	memset(m_inbuffer, 0, m_audiosize * sizeof(float));
	
	if (INRATE(1) == calc_FullRate) {
		m_numSamples = unit->mWorld->mFullRate.mBufLength;
	} else {
		m_numSamples = 1;
	}
	
	SETCALC(SonLPCSynthIn_next);
	//Print("done LPCSynth constructor\n");
}
SonLPCSynthIn::~SonLPCSynthIn(){
	//SCWorld_Allocator alloc(ft, unit->mWorld);
	//if(unit->m_scfft)
	//	scfft_destroy(unit->m_scfft, alloc);
	if(m_buffer)
		RTFree(mWorld, m_buffer);
	if(m_inbuffer)
		RTFree(mWorld, m_inbuffer);
	if(m_olabuf)
		RTFree(mWorld, m_olabuf);
	if(m_synwindow)
		RTFree(mWorld, m_synwindow);
	if(m_analwindow)
		RTFree(mWorld, m_analwindow);
}

void SonLPCSynthIn_next(SonLPCSynthIn *unit, int wronginNumSamples)
{
	
	float *out = OUT(0);
	float fbufnum = ZIN0(0);
	float G = IN(0)[1];
	int poles = IN(0)[2];
	int pos     = unit->m_pos;
	int audiosize = unit->m_audiosize;
	int numSamples = unit->m_numSamples;
	////////get input
	float *in = IN(1);
	//memmove(unit->m_inbuffer + unit->m_numSamples, unit->m_inbuffer ,(audiosize - unit->m_numSamples)* sizeof(float));
	//memcpy(unit->m_inbuffer, in ,unit->m_numSamples* sizeof(float));
	memmove(unit->m_inbuffer , unit->m_inbuffer + unit->m_numSamples ,(audiosize - unit->m_numSamples)* sizeof(float));
	memcpy(unit->m_inbuffer + (audiosize - unit->m_numSamples), in ,unit->m_numSamples* sizeof(float));
//Print("buffnum %f, poles %d, G: %f \n",fbufnum,poles,G);	
	// Only run the IFFT if we're receiving a new block of input data - otherwise just output data already received
	if (fbufnum >= 0.f){

		float* error = unit->m_inbuffer; 
		
		//float* error = unit->m_fftsndbuf->data;
		float* kaes = unit->m_fftsndbuf->data + audiosize;
		float coef[MAX_POLES + 2];
		calculalfa(coef, kaes, poles);
		//G = 1;
		sintesisD(unit->m_buffer, error, audiosize, coef, poles, G);
		// Then shunt the "old" time-domain output down by one hop
		//memcpy(unit->m_buffer,error,audiosize*sizeof(float));
		//for(int i = 0; i < unit->m_audiosize ; i++)
		//	unit->m_buffer[i] = error[i];// * unit->m_analwindow[i];
		//memset(error, 0, unit->m_audiosize * sizeof(float));
		
		//invfirstdiff(unit->m_buffer, audiosize);
		
		int hopsamps = pos;
		if(!unit->m_windowinited){
			makesynwindow(unit->m_synwindow, audiosize, hopsamps + unit->mWorld->mFullRate.mBufLength);
			unit->m_windowinited = true;
			//Print("hopsamps is %d audiosize is %d\n",hopsamps + unit->mWorld->mFullRate.mBufLength, audiosize);
		}
		int shuntsamps = audiosize - hopsamps;
		if(hopsamps != audiosize)  // There's only copying to be done if the position isn't all the way to the end of the buffer
			memmove(unit->m_olabuf, unit->m_olabuf+hopsamps, shuntsamps * sizeof(float));

		// Then mix the "new" time-domain data in - adding at first, then just setting (copying) where the "old" is supposed to be zero.

			// NB we re-use the "pos" variable temporarily here for write rather than read
			for(pos = 0; pos < shuntsamps; ++pos){
				unit->m_olabuf[pos] += unit->m_buffer[pos] * unit->m_synwindow[pos];
			}
		//memcpy(unit->m_olabuf + shuntsamps, fftbuf + shuntsamps, (hopsamps) * sizeof(float));
		for(pos = shuntsamps; pos < audiosize; ++pos){
				unit->m_olabuf[pos] = unit->m_buffer[pos] * unit->m_synwindow[pos];
		}
		// Move the pointer back to zero, which is where playback will next begin
		pos = 0;

	} // End of has-the-chain-fired

	// Now we can output some stuff, as long as there is still data waiting to be output.
	// If there is NOT data waiting to be output, we output zero. (Either irregular/negative-overlap
	//     FFT firing, or FFT has given up, or at very start of execution.)
	if(pos >= audiosize)
		ClearUnitOutputs(unit, numSamples);
	else {
		memcpy(out, unit->m_olabuf + pos, numSamples * sizeof(float));
		pos += numSamples;
	}
	unit->m_pos = pos;
}
////////////////////////////////
struct SonLPCError:public SonLPCBase
{
	SonLPCError(Unit* unit);
	~SonLPCError();
	float * m_olabuf;
	float * m_buffer;
	int m_pos;
	float * m_synwindow;
	bool m_windowinited;
};
SCWrapClass(SonLPCError);
SonLPCError::SonLPCError(Unit* unit):SonLPCBase(unit),m_pos(0),m_olabuf(NULL)
{
	m_synwindow = (float*)RTAlloc(unit->mWorld, m_audiosize * sizeof(float));
	m_windowinited = false;
	m_olabuf = (float*)RTAlloc(unit->mWorld, m_audiosize * sizeof(float));
	memset(m_olabuf, 0, m_audiosize * sizeof(float));
	m_buffer = (float*)RTAlloc(unit->mWorld, m_audiosize * sizeof(float));
	
	if (INRATE(1) == calc_FullRate) {
		m_numSamples = unit->mWorld->mFullRate.mBufLength;
	} else {
		m_numSamples = 1;
	}
	
	SETCALC(SonLPCError_next);
}
SonLPCError::~SonLPCError(){
	//SCWorld_Allocator alloc(ft, unit->mWorld);
	//if(unit->m_scfft)
	//	scfft_destroy(unit->m_scfft, alloc);
	if(m_buffer)
		RTFree(mWorld, m_buffer);
	if(m_olabuf)
		RTFree(mWorld, m_olabuf);
	if(m_synwindow)
		RTFree(mWorld, m_synwindow);
}

void SonLPCError_next(SonLPCError *unit, int wronginNumSamples)
{
	
	float *out = OUT(0);
	float fbufnum = ZIN0(0);
	float G = IN(0)[1];
	float poles = IN(0)[2];
	int pos     = unit->m_pos;
	int audiosize = unit->m_audiosize;
	int numSamples = unit->m_numSamples;
	
	// Only run the IFFT if we're receiving a new block of input data - otherwise just output data already received
	if (fbufnum >= 0.f){

		float* error = unit->m_fftsndbuf->data;
		float* kaes = unit->m_fftsndbuf->data + audiosize;
		
		/*
		float coef[MAX_POLES + 2];
		calculalfa(coef, kaes, poles);
		G = 1;
		sintesisD(unit->m_buffer, error, audiosize, coef, poles, G);
		*/
		memcpy(unit->m_buffer, error, audiosize * sizeof(float));
		// Then shunt the "old" time-domain output down by one hop
		
		int hopsamps = pos;
		if(!unit->m_windowinited){
			makesynwindow(unit->m_synwindow, audiosize, hopsamps + unit->mWorld->mFullRate.mBufLength);
			unit->m_windowinited = true;
			Print("hopsamps is %d\n",hopsamps + unit->mWorld->mFullRate.mBufLength);
		}
		int shuntsamps = audiosize - hopsamps;
		if(hopsamps != audiosize)  // There's only copying to be done if the position isn't all the way to the end of the buffer
			memmove(unit->m_olabuf, unit->m_olabuf+hopsamps, shuntsamps * sizeof(float));

		// Then mix the "new" time-domain data in - adding at first, then just setting (copying) where the "old" is supposed to be zero.

			// NB we re-use the "pos" variable temporarily here for write rather than read
			for(pos = 0; pos < shuntsamps; ++pos){
				unit->m_olabuf[pos] += unit->m_buffer[pos] * unit->m_synwindow[pos];
			}

		//memcpy(unit->m_olabuf + shuntsamps, fftbuf + shuntsamps, (hopsamps) * sizeof(float));
		for(pos = shuntsamps; pos < audiosize; ++pos){
				unit->m_olabuf[pos] = unit->m_buffer[pos] * unit->m_synwindow[pos];
		}
		// Move the pointer back to zero, which is where playback will next begin
		pos = 0;

	} // End of has-the-chain-fired

	// Now we can output some stuff, as long as there is still data waiting to be output.
	// If there is NOT data waiting to be output, we output zero. (Either irregular/negative-overlap
	//     FFT firing, or FFT has given up, or at very start of execution.)
	if(pos >= audiosize)
		ClearUnitOutputs(unit, numSamples);
	else {
		memcpy(out, unit->m_olabuf + pos, numSamples * sizeof(float));
		pos += numSamples;
	}
	unit->m_pos = pos;
}
///////////////////////////////////////////////////
struct SonLPCMorph:public SonLPCSynth
{
	SonLPCMorph(Unit *unit);
	SndBuf * m_buf2;
	uint32 m_buf2num;
};
SCWrapClass(SonLPCMorph);
SonLPCMorph::SonLPCMorph(Unit* unit):SonLPCSynth(unit)
{
	m_buf2num = (uint32)ZIN0(1);
	m_buf2 = LPCGetBuffer(unit, m_buf2num, "SonLPCMorph", 1);
	if(!m_fftsndbuf)
		return;
	SETCALC(SonLPCMorph_next);
}
void SonLPCMorph_next(SonLPCMorph *unit, int wronginNumSamples)
{
	
	float *out = OUT(0);
	float fbufnum = ZIN0(0);
	float G1 = IN(0)[1];
	float poles = IN(0)[2];
	float fbufnum2 = ZIN0(1);
	float G2 = IN(1)[1];
	float poles2 = IN(1)[2];
	int pos     = unit->m_pos;
	
	int audiosize = unit->m_audiosize;
	int numSamples = unit->m_numSamples;
	
	// Only run the IFFT if we're receiving a new block of input data - otherwise just output data already received
	if (fbufnum >= 0.f){

		float* error = unit->m_buf2->data;
		float* kaes = unit->m_fftsndbuf->data + audiosize;
		float coef[MAX_POLES + 2];
		calculalfa(coef, kaes, poles);
		float G;
		if(G2!=0)
			G = G1/G2;
		else
			G = 1;
		sintesisD(unit->m_buffer, error, audiosize, coef, poles, G);
		// Then shunt the "old" time-domain output down by one hop
		
		invfirstdiff(unit->m_buffer, audiosize);
		
		int hopsamps = pos;
		if(!unit->m_windowinited){
			makesynwindow(unit->m_synwindow, audiosize, hopsamps + unit->mWorld->mFullRate.mBufLength);
			unit->m_windowinited = true;
			Print("hopsamps is %d\n",hopsamps + unit->mWorld->mFullRate.mBufLength);
		}
		int shuntsamps = audiosize - hopsamps;
		if(hopsamps != audiosize)  // There's only copying to be done if the position isn't all the way to the end of the buffer
			memmove(unit->m_olabuf, unit->m_olabuf+hopsamps, shuntsamps * sizeof(float));

		// Then mix the "new" time-domain data in - adding at first, then just setting (copying) where the "old" is supposed to be zero.

			// NB we re-use the "pos" variable temporarily here for write rather than read
			for(pos = 0; pos < shuntsamps; ++pos){
				unit->m_olabuf[pos] += unit->m_buffer[pos] * unit->m_synwindow[pos];
			}

		//memcpy(unit->m_olabuf + shuntsamps, fftbuf + shuntsamps, (hopsamps) * sizeof(float));
		for(pos = shuntsamps; pos < audiosize; ++pos){
				unit->m_olabuf[pos] = unit->m_buffer[pos] * unit->m_synwindow[pos];
		}
		// Move the pointer back to zero, which is where playback will next begin
		pos = 0;

	} // End of has-the-chain-fired

	// Now we can output some stuff, as long as there is still data waiting to be output.
	// If there is NOT data waiting to be output, we output zero. (Either irregular/negative-overlap
	//     FFT firing, or FFT has given up, or at very start of execution.)
	if(pos >= audiosize)
		ClearUnitOutputs(unit, numSamples);
	else {
		memcpy(out, unit->m_olabuf + pos, numSamples * sizeof(float));
		pos += numSamples;
	}
	unit->m_pos = pos;
}
////////////////////////////////////////////////////
struct SonLPCSynthInput : public SonLPCBase
{
	float m_coef[MAX_POLES + 2];
	float m_kaes[MAX_POLES + 2];
	uint32 m_poles;
	float m_G;
	SndBuf * m_coefsndbuf;
	uint32 m_coefbufnum;
	LTIv *filter;
	SonLPCSynthInput(SonLPCSynthInput* unit);
};


SCWrapClassT(SonLPCSynthInput);

SonLPCSynthInput::SonLPCSynthInput(SonLPCSynthInput *unit):SonLPCBase(unit)
{
	m_coefbufnum = (uint32)ZIN0(0);
	float *out = OUT(0);
	float G = IN(0)[1];
	float poles = IN(0)[2];
	
	m_G = G;
	m_poles = poles;
	int audiosize = unit->m_audiosize;
	m_coefsndbuf = LPCGetBuffer(unit, m_coefbufnum, "SonLPCSynth", 1);
	if(!m_coefsndbuf)
		return;
	
	float* error = m_coefsndbuf->data;
	float* kaes = m_coefsndbuf->data + audiosize;
	//m_maxpoles = sc_min(MAX_POLES,m_coefsndbuf->samples);
	Print("poles %g data %g\n",poles,(float)m_coefsndbuf->samples);
	calculalfa(m_coef,kaes,poles);
	
	int nB = 1;//(int)IN0(1);
	int nA = poles ;//(int)IN0(2);
	filter = new(unit) LTIv(unit,nB,nA);
	for(int i=0;i<nB;i++)
		this->filter->KernelB[i] = 0;//G;//IN0(3+i);
	for(int i=0;i<nA;i++)
		this->filter->KernelA[i] = 0;//-m_coef[i + 1];//IN0(3+nB+i);
	SETCALC(SonLPCSynthInput_next);
}
void SonLPCSynthInput_next(SonLPCSynthInput *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(1);
	
	float fbufnum = ZIN0(0);
	float G = IN(0)[1];
	float poles = IN(0)[2];
	int useG = ZIN0(2);
	int audiosize = unit->m_audiosize;
	if (fbufnum >= 0.f){
		if (useG){
			G = G;
		}else{
			G = 1.0;
		}
		unit->m_poles = poles;
		unit->m_G = G;
		float* error = unit->m_coefsndbuf->data;
		float* kaes = unit->m_coefsndbuf->data + audiosize;
		//Print("fbufnum %f %d\n",fbufnum,unit->m_poles);
		//for(int i=0;i<=unit->m_poles;i++){
		//	Print("%d es %f\n",i,unit->m_coefsndbuf->data[i]);
		//}
		calculalfa(unit->m_coef,kaes,unit->m_poles);
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
	DefineDtorUnit(SonLPCSynthIn);
	DefineDtorUnit(SonLPCError);
	DefineDtorUnit(SonLPCSynthInput);
	DefineDtorUnit(SonLPCMorph);
}
