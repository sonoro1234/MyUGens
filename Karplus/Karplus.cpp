#include "SC_PlugIn.h"
//#include <algorithm>
static InterfaceTable *ft;
struct DelayUnit : public Unit
{
	float *m_dlybuf;

	float m_dsamp, m_fdelaylen;
	float m_delaytime, m_maxdelaytime;
	long m_iwrphase, m_idelaylen, m_mask;
	long m_numoutput;
};
struct FeedbackDelay : public DelayUnit
{
	float m_feedbk, m_decaytime;
};
struct CubicInterpolationUnit
{
	static const int minDelaySamples = 2;
};
struct CoefsFilterUnit
{
	float* m_coefsB;
	float* m_coefsA;
	float* m_valsX;
	float* m_valsY;
	long m_ncoefsB,m_ncoefsA,m_indexX,m_indexY,m_indexXmask,m_indexYmask;
};
struct AllPassInterp
{
	float lastY;
	float lastX;
};
struct Karplus : public FeedbackDelay, CubicInterpolationUnit,CoefsFilterUnit,AllPassInterp
{
	float m_lastsamp, m_prevtrig, m_coef;
	long m_inputsamps;
};

extern "C"
{
	void Karplus_Ctor(Karplus* unit);
	void Karplus_next_aa(Karplus *unit, int inNumSamples);
	void Karplus_next_aa_z(Karplus *unit, int inNumSamples);
	void Karplus_next_ka(Karplus *unit, int inNumSamples);
	void Karplus_next_ka_z(Karplus *unit, int inNumSamples);
	void Karplus_next_kk(Karplus *unit, int inNumSamples);
	void Karplus_next_kk_z(Karplus *unit, int inNumSamples);
}

///////////////////////////////
template <typename Unit>
static float CalcDelay(Unit *unit, float delaytime)
{
	float minDelay = Unit::minDelaySamples;
	float next_dsamp = delaytime * (float)SAMPLERATE;
	return sc_clip(next_dsamp, minDelay, unit->m_fdelaylen);
}
static float DoFilter(Karplus *unit,float val)
{
	float res = val * unit->m_coefsB[0];
	for(int i = 1; i< unit->m_ncoefsB; i++)
		res += unit->m_coefsB[i] * unit->m_valsX[(i + unit->m_indexX) & unit->m_indexXmask];
	unit->m_valsX[unit->m_indexX] = val;
	unit->m_indexX = (unit->m_indexX - 1) & unit->m_indexXmask;
	
	for(int i = 0; i< unit->m_ncoefsA; i++)
		res -= unit->m_coefsA[i] * unit->m_valsY[(i + unit->m_indexY) & unit->m_indexYmask];
	unit->m_indexY = (unit->m_indexY - 1) & unit->m_indexYmask;
	res = zapgremlins(res);
	unit->m_valsY[unit->m_indexY] = res;
	
	return res;	
}
//fractional delay with allpass
static float DoAllPass(float val,AllPassInterp *ap,float coef)
{
	//if(ap->lastY !=0.0f)
		//Print("ap.lastY %g",ap->lastY);
	//float mu = std::max(0.0f,val - 0.5f)*2 ;
	//float mu = std::max(0.0f,(float)fabs(val) - 0.1f)*1.0f/(1.0f-0.1f) ;
	//mu = std::max(0.0f,val) ;
	float mu =1.0f;
	if (fabs(val) > coef)
		mu = 0.0f;
	mu = (1-mu)/(1+mu);
	float res = mu*(val - ap->lastY) + ap->lastX;
	ap->lastX = val;
	ap->lastY = res;
	
	return res;
}
static bool DelayUnit_AllocDelayLine(DelayUnit *unit, const char * className)
{
	long delaybufsize = (long)ceil(unit->m_maxdelaytime * SAMPLERATE + 1.f);
	delaybufsize = delaybufsize + BUFLENGTH;
	delaybufsize = NEXTPOWEROFTWO(delaybufsize);  // round up to next power of two
	unit->m_fdelaylen = unit->m_idelaylen = delaybufsize;

	if (unit->m_dlybuf)
		RTFree(unit->mWorld, unit->m_dlybuf);
	unit->m_dlybuf = (float*)RTAlloc(unit->mWorld, delaybufsize * sizeof(float));
	
	//for(int i=0;i<delaybufsize;i++)
	//		unit->m_dlybuf[i]=0.f;
	memset(unit->m_dlybuf, 0, delaybufsize * sizeof(float));
	
	if (unit->m_dlybuf == NULL) {
		SETCALC(ft->fClearUnitOutputs);
		ClearUnitOutputs(unit, 1);

		if(unit->mWorld->mVerbosity > -2)
			Print("Failed to allocate memory for %s ugen.\n", className);
	}

	unit->m_mask = delaybufsize - 1;
	return (unit->m_dlybuf != NULL);
}
void Karplus_Dtor(Karplus *unit)
{
	RTFree(unit->mWorld, unit->m_dlybuf);
}

////////////////////////////////
static bool Karplus_setCoefs(Karplus* unit)
{
	unit->m_ncoefsA = IN0(6);
	unit->m_ncoefsB = IN0(7);
	long coefsbufsize = unit->m_ncoefsA + unit->m_ncoefsB;
	unit->m_coefsA = (float*)RTAlloc(unit->mWorld, (unit->m_ncoefsA + unit->m_ncoefsB) * sizeof(float));
	unit->m_coefsB = unit->m_coefsA + unit->m_ncoefsB;
	for(int i = 0;i < unit->m_ncoefsA; i++){
			unit->m_coefsA[i] = IN0(8 + i);
			//Print("coefA %d %g\n",i,unit->m_coefsA[i]);
	}
	for(int i = 0;i < unit->m_ncoefsB; i++){
			unit->m_coefsB[i] = IN0(8 + unit->m_ncoefsA + i);
			//Print("coefB %d %g\n",i,unit->m_coefsB[i]);
	}
	long Xbufsize = NEXTPOWEROFTWO(unit->m_ncoefsB); 
	long Ybufsize = NEXTPOWEROFTWO(unit->m_ncoefsA);
	unit->m_valsX = (float*)RTAlloc(unit->mWorld, (Xbufsize + Ybufsize) * sizeof(float));
	memset(unit->m_valsX,0,(Xbufsize + Ybufsize) * sizeof(float));
	unit->m_valsY = unit->m_valsX + Xbufsize;
	unit->m_indexX = unit->m_indexY = 0;
	unit->m_indexXmask = Xbufsize - 1;
	unit->m_indexYmask = Ybufsize - 1;
	//Print("Xbufsize %d m_indexXmask %d\n",Xbufsize,unit->m_indexXmask);
	//Print("Ybufsize %d m_indexYmask %d\n",Ybufsize,unit->m_indexYmask);
}
void Karplus_Ctor(Karplus *unit)
{
//	FeedbackDelay_Reset(unit);
	unit->m_maxdelaytime = IN0(2);
	unit->m_delaytime = IN0(3);
	unit->m_dlybuf = 0;
	bool allocationSucessful = DelayUnit_AllocDelayLine(unit, "Karplus");
	if (!allocationSucessful)
		return;
	Karplus_setCoefs(unit);
	
	unit->m_dsamp = CalcDelay(unit, unit->m_delaytime);

	unit->m_numoutput = 0;
	unit->m_iwrphase = 0;
	unit->m_feedbk = sc_CalcFeedback(unit->m_delaytime, unit->m_decaytime);
	
	unit->lastX=0.0f;
	unit->lastY=0.0f;
/*
	if (INRATE(1) == calc_FullRate) {
		if(INRATE(5) == calc_FullRate){
			SETCALC(Karplus_next_aa_z);
		} else {
			SETCALC(Karplus_next_ak_z); //ak
		}
	} else {
	    if(INRATE(5) == calc_FullRate){
			SETCALC(Karplus_next_ka_z); //ka
		} else {
			SETCALC(Karplus_next_kk_z); //kk
		}
	}
	*/
	SETCALC(Karplus_next_kk);
	OUT0(0) = unit->m_lastsamp = 0.f;
	unit->m_prevtrig = 0.f;
	unit->m_inputsamps = 0;
	unit->m_coef = IN0(5);
}
void Karplus_next_aa(Karplus *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float *trig = IN(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float *coef = IN(5);
	float lastsamp = unit->m_lastsamp;
	unsigned long inputsamps = unit->m_inputsamps;

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float thisin, curtrig;
	float prevtrig = unit->m_prevtrig;

	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
	for(int i = 0; i < inNumSamples; i++){
			curtrig = trig[i];
			if ((prevtrig <= 0.f) && (curtrig > 0.f)) {
			    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
			    }
			prevtrig = curtrig;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;
			out[i] = lastsamp = onepole;
			iwrphase++;
		};
	} else {

		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);

		float next_feedbk = sc_CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);

	    for(int i = 0; i < inNumSamples; i++){
			curtrig = trig[i];
			if ((prevtrig <= 0.f) && (curtrig > 0.f)) {
			    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
			    }
			prevtrig = curtrig;
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;
			out[i] = lastsamp = onepole;
			feedbk += feedbk_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
	}

	unit->m_prevtrig = prevtrig;
	unit->m_inputsamps = inputsamps;
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

}
void Karplus_next_aa_z(Karplus *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float *trig = IN(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float *coef = IN(5);
	float lastsamp = unit->m_lastsamp;

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float d0, d1, d2, d3;
	float thisin, curtrig;
	unsigned long inputsamps = unit->m_inputsamps;
	float prevtrig = unit->m_prevtrig;

	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		for(int i = 0; i < inNumSamples; i++){
			curtrig = trig[i];
			if ((prevtrig <= 0.f) && (curtrig > 0.f)) {
			    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
			    }
			prevtrig = curtrig;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin;
				out[i] = 0.f;
			} else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;
			out[i] = lastsamp = onepole;
			}
			iwrphase++;
		};
	} else {

		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);

		float next_feedbk = sc_CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);

		for(int i = 0; i < inNumSamples; i++) {
			curtrig = trig[i];
			if ((prevtrig <= 0.f) && (curtrig > 0.f)) {
			    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
			    }
			prevtrig = curtrig;
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin;
				out[i] = 0.f;
			    } else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;
			out[i] = lastsamp = onepole;
			}
			feedbk += feedbk_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
	}

	unit->m_inputsamps = inputsamps;
	unit->m_prevtrig = prevtrig;
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

	unit->m_numoutput += inNumSamples;
	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(Karplus_next_aa);
	}
}
void Karplus_next_ka(Karplus *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float trig = IN0(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float *coef = IN(5);
	float lastsamp = unit->m_lastsamp;
	unsigned long inputsamps = unit->m_inputsamps;

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float thisin;

	if ((unit->m_prevtrig <= 0.f) && (trig > 0.f)) {
	    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
	    }
	unit->m_prevtrig = trig;

	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
	for(int i = 0; i < inNumSamples; i++){
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;
			out[i] = lastsamp = onepole;
			iwrphase++;
		};
	} else {

		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);

		float next_feedbk = sc_CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);

	    for(int i = 0; i < inNumSamples; i++){
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;
			out[i] = lastsamp = onepole;
			feedbk += feedbk_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
	}

	unit->m_inputsamps = inputsamps;
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

}


void Karplus_next_ka_z(Karplus *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float trig = IN0(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float *coef = IN(5);
	float lastsamp = unit->m_lastsamp;

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float d0, d1, d2, d3;
	float thisin;
	unsigned long inputsamps = unit->m_inputsamps;

	if ((unit->m_prevtrig <= 0.f) && (trig > 0.f)) {
	    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
	    }

	unit->m_prevtrig = trig;

	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		for(int i = 0; i < inNumSamples; i++){
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin;
				out[i] = 0.f;
			} else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;
			out[i] = lastsamp = onepole;
			}
			iwrphase++;
		};
	} else {

		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);

		float next_feedbk = sc_CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);

		for(int i = 0; i < inNumSamples; i++) {
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin;
				out[i] = 0.f;
			    } else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float thiscoef = coef[i];
			float onepole = ((1. - fabs(thiscoef)) * value) + (thiscoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + feedbk * onepole;
			out[i] = lastsamp = onepole;
			}
			feedbk += feedbk_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
	}

	unit->m_inputsamps = inputsamps;
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

	unit->m_numoutput += inNumSamples;
	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(Karplus_next_ka);
	}
}

void Karplus_next_kk(Karplus *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float trig = IN0(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float coef = IN0(5);
	//float lastsamp = unit->m_lastsamp;
	unsigned long inputsamps = unit->m_inputsamps;

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float thisin;

	if ((unit->m_prevtrig <= 0.f) && (trig > 0.f)) {
	    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
	    }
	unit->m_prevtrig = trig;

	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime && coef == unit->m_coef) {
		long idsamp = (long)dsamp ;
		float frac = dsamp - idsamp;

	    for(int i = 0; i < inNumSamples; i++){
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			/*
			float onepole = ((1. - fabs(coef)) * value) + (coef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + (feedbk * onepole);
			out[i] = lastsamp = onepole; //value;
			*/
			//if(coef > 0.1){
				//value = atan(coef*value)/atan(coef);
				//value = value - value*value*value/3;
			//}
			/*
			if (value > 0  && value > coef)
				value = coef;
			else if(value < 0 and value < -coef)
				value = - coef;
			*/
			
			float fvalue = DoFilter(unit,value) * feedbk;
			//if(coef > 0.5)
				//fvalue = DoAllPass(fvalue,unit,coef);
			dlybuf[iwrphase & mask] = thisin + fvalue;
			out[i]  = fvalue; 
			
			iwrphase++;
		};
	} else {

		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);

		float next_feedbk = sc_CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);

		float curcoef = unit->m_coef;
		float coef_slope = CALCSLOPE(coef, curcoef);

	    for(int i = 0; i < inNumSamples; i++){
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			float d0 = dlybuf[irdphase0 & mask];
			float d1 = dlybuf[irdphase1 & mask];
			float d2 = dlybuf[irdphase2 & mask];
			float d3 = dlybuf[irdphase3 & mask];
			float value = cubicinterp(frac, d0, d1, d2, d3);
			/*
			float onepole = ((1. - fabs(curcoef)) * value) + (curcoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + (feedbk * onepole);
			out[i] = lastsamp = onepole; //value;
			*/
			//if(coef > 0.1){
				//value = atan(coef*value)/atan(coef);
				//value = value - pow(value,3)/3;
				//value = value - value*value*value/3;
			//}
			/*
			if (value > 0  && value > coef)
				value = coef;
			else if(value < 0 and value < -coef)
				value = - coef;
			*/
			
			float fvalue = DoFilter(unit,value) * feedbk;
			//if(coef > 0.5)
				//fvalue = DoAllPass(fvalue,unit,coef);
			dlybuf[iwrphase & mask] = thisin + fvalue;
			out[i]  = fvalue; 
			
			feedbk += feedbk_slope;
			curcoef += coef_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_coef = coef;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
	}

	unit->m_inputsamps = inputsamps;
	//unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

}


void Karplus_next_kk_z(Karplus *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	float trig = IN0(1);
	float delaytime = IN0(3);
	float decaytime = IN0(4);
	float coef = IN0(5);
	float lastsamp = unit->m_lastsamp;

	float *dlybuf = unit->m_dlybuf;
	long iwrphase = unit->m_iwrphase;
	float dsamp = unit->m_dsamp;
	float feedbk = unit->m_feedbk;
	long mask = unit->m_mask;
	float d0, d1, d2, d3;
	float thisin;
	unsigned long inputsamps = unit->m_inputsamps;

	if ((unit->m_prevtrig <= 0.f) && (trig > 0.f)) {
	    inputsamps = (long)(delaytime * unit->mRate->mSampleRate + .5f);
	    }
	unit->m_prevtrig = trig;

	if (delaytime == unit->m_delaytime && decaytime == unit->m_decaytime && coef == unit->m_coef) {
		long idsamp = (long)dsamp;
		float frac = dsamp - idsamp;
		for(int i = 0; i < inNumSamples; i++){

			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin;
				out[i] = 0.f;
			} else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
				float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(coef)) * value) + (coef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + (feedbk * onepole);
			out[i] = lastsamp = onepole; //value;
			}
			iwrphase++;
		};
	} else {

		float next_dsamp = CalcDelay(unit, delaytime);
		float dsamp_slope = CALCSLOPE(next_dsamp, dsamp);

		float next_feedbk = sc_CalcFeedback(delaytime, decaytime);
		float feedbk_slope = CALCSLOPE(next_feedbk, feedbk);

		float curcoef = unit->m_coef;
		float coef_slope = CALCSLOPE(coef, curcoef);

		for(int i = 0; i < inNumSamples; i++) {
			dsamp += dsamp_slope;
			long idsamp = (long)dsamp;
			float frac = dsamp - idsamp;
			long irdphase1 = iwrphase - idsamp;
			long irdphase2 = irdphase1 - 1;
			long irdphase3 = irdphase1 - 2;
			long irdphase0 = irdphase1 + 1;
			if (inputsamps > 0) {
			    thisin = in[i];
			    --inputsamps;
			    } else {
			    thisin = 0.f;
			    }
			if (irdphase0 < 0) {
				dlybuf[iwrphase & mask] = thisin;
				out[i] = 0.f;
			    } else {
				if (irdphase1 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
				} else if (irdphase2 < 0) {
					d1 = d2 = d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
				} else if (irdphase3 < 0) {
					d3 = 0.f;
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
				} else {
					d0 = dlybuf[irdphase0 & mask];
					d1 = dlybuf[irdphase1 & mask];
					d2 = dlybuf[irdphase2 & mask];
					d3 = dlybuf[irdphase3 & mask];
				}
			float value = cubicinterp(frac, d0, d1, d2, d3);
			float onepole = ((1. - fabs(curcoef)) * value) + (curcoef * lastsamp);
			dlybuf[iwrphase & mask] = thisin + (feedbk * onepole);
			out[i] = lastsamp = onepole; //value;
			}
			feedbk += feedbk_slope;
			curcoef += coef_slope;
			iwrphase++;
		};
		unit->m_feedbk = feedbk;
		unit->m_dsamp = dsamp;
		unit->m_delaytime = delaytime;
		unit->m_decaytime = decaytime;
		unit->m_coef = coef;
	}

	unit->m_inputsamps = inputsamps;
	unit->m_lastsamp = zapgremlins(lastsamp);
	unit->m_iwrphase = iwrphase;

	unit->m_numoutput += inNumSamples;
	if (unit->m_numoutput >= unit->m_idelaylen) {
		SETCALC(Karplus_next_kk);
	}
}

PluginLoad(Karplus){
	ft = inTable;
	#define DefineDelayUnit(name) \
	(*ft->fDefineUnit)(#name, sizeof(name), (UnitCtorFunc)&name##_Ctor, \
	(UnitDtorFunc)&DelayUnit_Dtor, 0);
	
	DefineDtorUnit(Karplus);

}