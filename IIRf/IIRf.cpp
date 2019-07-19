#include "SC_PlugIn.h"
#include <float.h>
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

static InterfaceTable *ft;
#include "DWG.cpp"

////////////////////////////////////////////////////
struct IIRf : public Unit
{
	LTIvT<double> *filter;
	IIRf(Unit* unit);
};

extern "C"
{
	void IIRf_Ctor(IIRf* unit);
	void IIRf_next(IIRf *unit, int inNumSamples);
	void IIRf_Dtor(IIRf *unit);
}

IIRf::IIRf(Unit *unit)
{
	int nB = (int)IN0(1);
	int nA = (int)IN0(2);
	filter = new(unit) LTIvT<double>(unit,nB,nA);
	for(int i=0;i<nB;i++)
		this->filter->KernelB[i] = IN0(3+i);
	for(int i=0;i<nA;i++)
		this->filter->KernelA[i] = IN0(3+nB+i);
	SETCALC(IIRf_next);
	(unit->mCalcFunc)(unit, 1);
}

void IIRf_Ctor(IIRf* unit)
{
	new(unit) IIRf(unit);
}
void IIRf_Dtor(IIRf* unit)
{
	unit->~IIRf();
}
void IIRf_next(IIRf *unit, int inNumSamples)
{
	float *out = OUT(0);
	float *in = IN(0);
	//coef changes
	int nB = unit->filter->kernel_sizeB;
	int nA = unit->filter->kernel_sizeA;
	for(int i=0;i<nB;i++)
		unit->filter->KernelB[i] = IN0(3+i);
	for(int i=0;i<nA;i++)
		unit->filter->KernelA[i] = IN0(3+nB+i);
		
	for (int i=0; i < inNumSamples; ++i)
	{
		out[i] = unit->filter->filter(in[i]);
	}

}

/////////////////////////////////////////
PluginLoad(IIRs)
{
	ft = inTable;
	DefineDtorUnit(IIRf);
}
