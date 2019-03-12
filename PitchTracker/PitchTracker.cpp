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
#include "../DWGUGens/dwglib/DWG.cpp"
//#include "SC_PlugIn.h"
//#define MAXDELAY 1024
InterfaceTable *ft;


//////////////////////////////////////////////////////////
struct PitchTracker:public Unit
{
	SndBuf *m_coefsndbuf;
	int m_pos;// m_coefbufsize; // "fullbufsize" includes any zero-padding, "audiosize" does not.
	int m_coefbufnum;
	//uint32 m_maxpoles,m_poles;
	int m_hopsize, m_shuntsize; // These add up to m_audiosize

	int m_numSamples;
	World *mWorld;
	
	PitchTracker(Unit* unit);
	~PitchTracker();
	float * m_inbuf;
	//int m_inbufpos;
	int m_audiosize;
	
	scfft* m_scfft;
	scfft* m_iscfft;
	//float m_coef[MAX_POLES + 2];
	//float m_kaes[MAX_POLES + 2];
	float *m_fftbuf;
	float *m_nsdf;
	int m_fftbufsize;
	float m_freq,m_efreq,m_clarity;
	float m_threshold;
	int m_kind; //0 ACR, 1 NSDF, 2SDF
	int m_useefreq;//0 no, 1 yes
	bool m_dodump;
	float R0;
};

SCWrapClass(PitchTracker);

PitchTracker::PitchTracker(Unit* unit):m_inbuf(NULL),m_pos(0){

	mWorld = unit->mWorld;
	
	m_pos = 0;
	m_efreq = m_freq = ZIN0(3);
	m_clarity = 0.0;
	R0 = 0.0;
	m_dodump = false;
	m_audiosize = (int)ZIN0(1);
	m_fftbufsize = NEXTPOWEROFTWO(m_audiosize);
	m_fftbuf = (float *)RTAlloc(unit->mWorld,m_fftbufsize * sizeof(float));
	m_audiosize = m_fftbufsize /2 ;
	
	Print("m_audiosize %f , %d\n",ZIN0(1),m_audiosize);
	m_inbuf = (float *)RTAlloc(unit->mWorld,m_audiosize * sizeof(float));
	memset(m_inbuf, 0, m_audiosize * sizeof(float));
	
	m_nsdf = (float *)RTAlloc(unit->mWorld,m_audiosize * sizeof(float));
	memset(m_nsdf, 0, m_audiosize * sizeof(float));
	
	m_kind = (int)ZIN0(5);
	m_useefreq = (int)ZIN0(6);
	
	m_coefbufnum = (int)ZIN0(7);
	Print("PitchTracker m_audiosize %d fftsize %d m_coefbufnum %d\n",m_audiosize,m_fftbufsize,m_coefbufnum);
	m_coefsndbuf = NULL;
	if (m_coefbufnum > -1)
		m_coefsndbuf = ConvGetBuffer(unit, m_coefbufnum, "PitchTracker", 1);
	if(!m_coefsndbuf){
		m_coefbufnum = -1;
		Print("m_coefsndbuf %p\n",m_coefsndbuf);
	}else{
		Print("m_coefsndbuf->samples %d\n",m_coefsndbuf->samples);
		if(m_coefsndbuf->samples < m_audiosize){
			m_coefbufnum = -1;
			Print("buffer smaller than m_audiosize\n");
		}
	}
	
	int hopsize = (int)(sc_max(sc_min(ZIN0(2), 1.f), 0.f) * m_audiosize);
	if (hopsize < unit->mWorld->mFullRate.mBufLength) {
		Print("PitchTracker_Ctor: hopsize smaller than SC's block size (%i) - automatically corrected.\n", hopsize, unit->mWorld->mFullRate.mBufLength);
		hopsize = unit->mWorld->mFullRate.mBufLength;
	} else if (((int)(hopsize / unit->mWorld->mFullRate.mBufLength)) * unit->mWorld->mFullRate.mBufLength
				!= hopsize) {
		Print("PitchTracker_Ctor: hopsize (%i) not an exact multiple of SC's block size (%i) - automatically corrected.\n", hopsize, unit->mWorld->mFullRate.mBufLength);
		hopsize = ((int)(hopsize / unit->mWorld->mFullRate.mBufLength)) * unit->mWorld->mFullRate.mBufLength;
	}
	
	m_hopsize = hopsize;
	m_shuntsize = m_audiosize - m_hopsize;
	Print("anal hopsize %d\n",hopsize);
	if (INRATE(0) == calc_FullRate) {
		m_numSamples = unit->mWorld->mFullRate.mBufLength;
	} else {
		m_numSamples = 1;
	}
	
	SCWorld_Allocator alloc(ft, unit->mWorld);
	m_scfft = scfft_create(m_fftbufsize, m_audiosize,(SCFFT_WindowFunction) -1, m_inbuf,m_fftbuf, kForward, alloc);
	m_iscfft = scfft_create(m_fftbufsize, m_audiosize, (SCFFT_WindowFunction) -1, m_fftbuf, m_fftbuf, kBackward, alloc);
	//m_scfft = fftwf_plan_r2r_1d(size, unit->autocorrTime, unit->autocorrFFT, FFTW_R2HC, FFTW_ESTIMATE);
	//m_iscfft = fftwf_plan_r2r_1d(size, unit->autocorrFFT, unit->autocorrTime, FFTW_HC2R, FFTW_ESTIMATE);
	SETCALC(PitchTracker_next);
	
	ZOUT0(0) = m_freq;//ZIN0(0);
	ZOUT0(1) = m_clarity;
	ZOUT0(2) = R0;
	//OUT(0)[1] = m_coefbufnum;
	//OUT(0)[2] = 0;
	return ;
}

PitchTracker::~PitchTracker(){
	SCWorld_Allocator alloc(ft, mWorld);
	if(m_scfft)
		scfft_destroy(m_scfft, alloc);
	if(m_iscfft)
		scfft_destroy(m_iscfft, alloc);
	if(m_nsdf)
		RTFree(mWorld, m_nsdf);
	if(m_inbuf)
		RTFree(mWorld, m_inbuf);
	if(m_fftbuf)
		RTFree(mWorld, m_fftbuf);

}

/////////////////////////////////////////////////////////
float parabol_interpolation(float *data,int maxpos,float*maxval)
{
	float a = (data[maxpos+1]+data[maxpos-1])*0.5-data[maxpos];
	float b = (data[maxpos+1]-data[maxpos-1])*0.5;
	float maxposf = (float)maxpos;
	*maxval = data[maxpos];
	if(a >0.000001f){
		maxposf = maxposf-b/(2*a);
		*maxval = data[maxpos] + b*maxposf + a*maxposf*maxposf;
	}
	return maxposf;
}
void PeakFindFreqMax(PitchTracker *unit,float *data)
{
	float exppos = (float)FULLRATE/unit->m_efreq;
	int lpos = sc_max(exppos*0.75,1);
	int rpos = sc_min(exppos*1.45,unit->m_fftbufsize/2);
	float maxval = 0;
	int maxpos = 1;
	for(int i=lpos; i<rpos;i++){
		//if peak
		if ((data[i-1] < data[i]) && (data[i+1] < data[i]))
			if(data[i] > maxval){
				maxval = data[i];
				maxpos = i;
			}
	}
	//parabollic interp
	float maxposf = parabol_interpolation(data,maxpos,&maxval);
	unit->m_dodump = false;
	if(maxpos > 700 || maxposf < 0 || maxpos < 2){
		Print("PeakFindFreqMax maxpos %d lpos %d rpos %d\n",maxpos,lpos,rpos);
		unit->m_dodump = true;
	}
	unit->m_clarity = maxval;
	unit->m_freq = FULLRATE/maxposf;///SAMPLERATE;
}
void PeakFindFreqMin(PitchTracker *unit,float *data)
{
	float exppos = (float)FULLRATE/unit->m_efreq;
	int lpos = sc_max(exppos*0.75,1);
	int rpos = sc_min(exppos*1.45,unit->m_fftbufsize/2);
	float maxval = 100000;
	int maxpos = 1;
	for(int i=lpos; i<rpos;i++){
		//if peak
		if ((data[i-1] > data[i]) && (data[i+1] > data[i]))
			if(data[i] < maxval){
				maxval = data[i];
				maxpos = i;
			}
	}
	//parabollic interp
	float maxposf = parabol_interpolation(data,maxpos,&maxval);
	unit->m_clarity = maxval;
	unit->m_freq = FULLRATE/maxposf;///SAMPLERATE;
}
void PeakFind(PitchTracker *unit,float *data)
{
	if(unit->m_useefreq>0)
		return PeakFindFreqMax(unit,data);
	//get maxpeak
	float maxvalall = 0;
	float prevmaxvalall = 0;
	int maxpos = 1;
	int prevmaxpos = 1;
	bool zeropass = false;
	int inipos = 1;
	int lastpos = unit->m_audiosize*0.95;
	for(int i=1; i<unit->m_audiosize*0.95;i++){
		if(!zeropass){ 
			//if(data[i] > 0 && data[i-1] <=0) //zero crosing pos slope
			if ((data[i-1] > data[i]) && (data[i+1] >= data[i])){ //antipeak
				zeropass = true;
				inipos = i;
			}
		}else{
			//if peak
			if ((data[i-1] < data[i]) && (data[i+1] <= data[i]))
				if(data[i] > prevmaxvalall){
					prevmaxvalall = data[i];
					prevmaxpos = i;
				}
				//if negative antipeak accept prevmaxval
			if((data[i-1] > data[i]) && (data[i+1] >= data[i])){
					maxvalall = prevmaxvalall;
					maxpos = prevmaxpos;
					lastpos = i;
			}
		}
	}
	//search for peak,skip 0
	float maxval = 0;
	maxpos = 1;
	for(int i=inipos; i<lastpos;i++){
		//if peak
		if ((data[i-1] < data[i]) && (data[i+1] <= data[i]))
			if(data[i] > maxval){
				maxval = data[i];
				maxpos = i;
			}
		if(((data[i-1] > data[i]) && (data[i+1] >= data[i])) and maxval >= unit->m_threshold*maxvalall)
			break;

	}
	

	//parabollic interp
	float maxposf = parabol_interpolation(data,maxpos,&maxval);
	unit->m_clarity = maxval;
	unit->m_freq = FULLRATE/maxposf;///SAMPLERATE;
	if(unit->R0 < 1e-6){
		unit->m_clarity = 0;
		unit->m_freq = unit->m_efreq;///SAMPLERATE;
	}
	
	unit->m_dodump = false;
	if( maxpos ==1 || unit->m_freq> 1000){
		Print("maxpos %d %f inipos %d lastpos %d m_freq %f R0 %f maxvalall %f prevmaxvalall %f pos %d\n",maxpos,maxval,inipos,lastpos,unit->m_freq, unit->R0,maxvalall,prevmaxvalall,prevmaxpos);
		unit->m_dodump = true;
	}
}
float PeakFindBAK(PitchTracker *unit,float *data)
{
	//get maxpeak
	float maxval = 0;
	int maxpos = 1;
	bool zeropass = false;
	int inipos = 1;
	for(int i=1; i<unit->m_fftbufsize/2;i++){
		if(!zeropass){ 
			//if(data[i] > 0 && data[i-1] <=0) //zero crosing pos slope
			if ((data[i-1] > data[i]) && (data[i+1] > data[i])){ //antipeak
				zeropass = true;
				inipos = i;
			}
		}else{
			//if peak
			if ((data[i-1] < data[i]) && (data[i+1] < data[i]))
				if(data[i] > maxval){
					maxval = data[i];
					maxpos = i;
				}
			if(data[i] < 0 and maxval >= unit->m_threshold)
				break;
		}
	}

	//parabollic interp
	float a = (data[maxpos+1]+data[maxpos-1])*0.5-data[maxpos];
	float b = (data[maxpos+1]-data[maxpos-1])*0.5;
	float maxposf = (float)maxpos;
	if(a !=0.0f){
		maxposf = maxposf-b/(2*a);
	}
	
	unit->m_freq = FULLRATE/maxposf;///SAMPLERATE;
}

void PeakFindSDF(PitchTracker *unit,float *data)
{
	if(unit->m_useefreq>0)
		return PeakFindFreqMin(unit,data);
	
	//get minpeak
	float maxval = 100000; //huge
	int maxpos = 1;
	bool zeropass = false;
	int inipos = 1;
	for(int i=1; i<unit->m_fftbufsize/2;i++){
		if(!zeropass){ 
			//if(data[i] > 0 && data[i-1] <=0) //zero crosing pos slope
			if ((data[i-1] < data[i]) && (data[i+1] < data[i])){ //peak
				zeropass = true;
				inipos = i;
			}
		}else{
			//if antipeak
			if ((data[i-1] > data[i]) && (data[i+1] > data[i]))
				if(data[i] < maxval){
					maxval = data[i];
					maxpos = i;
				}
			if(((data[i-1] < data[i]) && (data[i+1] < data[i])) and maxval <= unit->m_threshold)
				break;
		}
	}

	//parabollic interp
	float a = (data[maxpos+1]+data[maxpos-1])*0.5-data[maxpos];
	float b = (data[maxpos+1]-data[maxpos-1])*0.5;
	float maxposf = (float)maxpos;
	if(a !=0.0f){
		maxposf = maxposf-b/(2*a);
	}
	
	unit->m_freq = FULLRATE/maxposf;///SAMPLERATE;
}

void PitchTracker_next(PitchTracker *unit, int wronginNumSamples)
{
	//float *out = OUT(0);
	float *in = IN(0);
	unit->m_efreq = ZIN0(3);
	unit->m_threshold = ZIN0(4);
	unit->m_kind = (int)ZIN0(5);
	unit->m_useefreq = (int)ZIN0(6);
	//int poles = (uint32)sc_min(unit->m_maxpoles,sc_max(4, ZIN0(3)));
	int numSamples = unit->m_numSamples;


	float *out = unit->m_inbuf + unit->m_pos + unit->m_shuntsize;
	// copy input
	memcpy(out, in, numSamples * sizeof(float));
	unit->m_pos += numSamples;
	
	if (unit->m_pos != unit->m_hopsize){// || !unit->m_fftsndbuf->data || unit->m_fftsndbuf->samples != unit->m_fullbufsize) {
		if(unit->m_pos == unit->m_hopsize)
			unit->m_pos = 0;
		//Print("SonLPC m_pos %u\n",unit->m_pos);
		ZOUT0(0) = unit->m_freq;
		ZOUT0(1) = unit->m_clarity;
		ZOUT0(2) = unit->R0;
	} else {

		unit->m_pos = 0; 
		/*
		//real autocorrelation
		float sum=0;
		for(int i=0; i < unit->m_audiosize;i++)
			sum += unit->m_inbuf[i]*unit->m_inbuf[i];
		*/
		//perform autocorrelation
		scfft_dofft(unit->m_scfft);
		
		int numbins = (unit->m_fftbufsize - 2) >> 1;
		float *data = unit->m_fftbuf;
		data[0] = data[0]*data[0];
		data[1] = data[1]*data[1];
		for(int i=1 ; i < numbins ; i++) {
			int re = i*2;
			int im = re +1;
			data[re] = data[re]*data[re] + data[im]*data[im];
			data[im] = 0.0f;
		}
		scfft_doifft(unit->m_iscfft);
		unit->R0 = data[0]/(float)unit->m_audiosize;
		//Print("R0 %f autoc[0] %f\n",unit->R0,sum);
		
		if(unit->m_kind==0){ //ACR
		
			//normalize
			float data0 = data[0];
			for(int i=0; i<unit->m_fftbufsize;i++){
				data[i] = data[i]/data0;
			}
			
			PeakFind(unit,data);
		
		}else if(unit->m_kind==1){ //NSDF
			//normalize
			// float data0 = data[0];
			// for(int i=0; i<unit->m_fftbufsize;i++){
				// data[i] = data[i]/data0;
			// }
			
			//nsdf
			double sumRightSq = data[0], sumLeftSq = data[0];
			float * input= unit->m_inbuf; 
			float * output= unit->m_nsdf; 
		
			for(int j=0; j<unit->m_audiosize; ++j) {
				float left= input[unit->m_audiosize-1-j];
				float right= input[j];
				
				sumLeftSq  -= left*left;///data[0];
				sumRightSq -= right*right;///data[0];
				//output[j] = 1.0 - (sumLeftSq + sumRightSq - 2*output[j]) / (sumLeftSq + sumRightSq);
				output[j] = 2.f*data[j] / static_cast<float>(sumLeftSq + sumRightSq);
				//output[j] = (static_cast<float>(sumLeftSq + sumRightSq) - 2*data[j]);///data[0];
			}
			
			data = output;
			
			PeakFind(unit,data);
		}else{//if(unit->m_kind==1){ //SDF
			//normalize
			// float data0 = data[0];
			// for(int i=0; i<unit->m_fftbufsize;i++){
				// data[i] = data[i]/data0;
			// }
			
			//nsdf
			double sumRightSq = data[0], sumLeftSq = data[0];
			float * input= unit->m_inbuf; 
			float * output= unit->m_nsdf; 
		
			for(int j=0; j<unit->m_audiosize; ++j) {
				float left= input[unit->m_audiosize-1-j];
				float right= input[j];
				
				sumLeftSq  -= left*left;///data[0];
				sumRightSq -= right*right;///data[0];
				//output[j] = 1.0 - (sumLeftSq + sumRightSq - 2*output[j]) / (sumLeftSq + sumRightSq);
				//output[j] = 2.f*data[j] / static_cast<float>(sumLeftSq + sumRightSq);
				output[j] = (static_cast<float>(sumLeftSq + sumRightSq) - 2*data[j]);///data[0];
			}
			
			data = output;
			
			PeakFindSDF(unit,data);
		}
		
		
		
		//if having sndbuf coppy data
		if(unit->m_coefbufnum > -1){
			SndBuf *buf = unit->m_coefsndbuf;
			if (buf->data && unit->m_dodump){
			int limit = sc_min(buf->samples,unit->m_audiosize);
			LOCK_SNDBUF(buf);
			for(int i=0; i<limit;i++){
				buf->data[i]=data[i];
			}
			}
		}

		/////////////////////////////////////////////////////
		memmove(unit->m_inbuf, unit->m_inbuf + unit->m_hopsize, unit->m_shuntsize * sizeof(float));
		//Print("SonLPC buffnum %u\n",unit->m_coefbufnum);
		ZOUT0(0) = unit->m_freq;//unit->m_coefbufnum;
		ZOUT0(1) = unit->m_clarity;
		ZOUT0(2) = unit->R0;
	}
}


/////////////////////////////////////////
PluginLoad(SonLPC)
{
	ft = inTable;
	DefineDtorUnit(PitchTracker);
}
