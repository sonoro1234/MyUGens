#define _GNU_SOURCE
#include <fenv.h>
#include "SC_PlugIn.h"
#include <float.h>
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif
#define TWOPI 2*M_PI
#define C_AIRE 331.45f	   //m/s  *sqrt(T/T0)
#define RO_AIRE 1.2929f //
#define MU_AIRE 1.846e-5f
#define PRANDTL .707f //0.737 //pagina 514 priston
#define HEAT_RATIO 1.402f   //Cp/Cv
#define CP_AIRE 1004.6f //julios/kg*K
#define TERM_COND_AIRE 2.624e-2f //elliot .02346
#define FRECMAX 22050.f//*4.f
#define PRESION_ATMOSFERICA 1.01325e5
//#include "../server/supernova/utilities/malloc_aligned.hpp"
//#include "simd_memory.hpp"
#include "simd_binary_arithmetic.hpp"
#include "simd_horizontal_functions.hpp"
//#include "simd_mix.hpp"
//#include <cstring>
//#include <stdio.h>

static InterfaceTable *ft;

typedef struct vec2dtag
{
	float x;
	float y;
}VECTOR2D;

///////////////////////////////////////////
//returns Sum(i=0,n-1) A[i]*B[i]
//Buf, A and B must be memory_aligned
float convolve_reversed(float *Buf,float *A,float *B,int n,int aligj)
{
	const unsigned int per_loop = nova::vec<float>::objects_per_cacheline;
	//int precount = ((intptr_t)A % (per_loop * sizeof(float)))/sizeof(float);

	int vsize = nova::vec<float>::size;
    int alig = ((intptr_t)A /sizeof(float))% vsize;
    int precount = (vsize - alig) % vsize;
    precount = sc_min(precount,n);
    int nmultiple = ((n - precount)/ per_loop)*per_loop;
	int nrest = (n - precount) - nmultiple;
	float sum = 0;


	for( int i = 0; i < precount; i++){
		sum += A[i]*B[i];
	}

    if(nmultiple > 0){
        nova::times_vec_simd(Buf,A + precount,B + precount,nmultiple);
        sum = nova::horizontal_sum_vec_simd(Buf,nmultiple);
    }

	if( nrest > 0){
		nova::times_vec(Buf,A + nmultiple + precount,B + nmultiple + precount, nrest);
		sum += nova::horizontal_sum_vec(Buf, nrest);
	}
	return sum;
}
float convolve_reversedS(float *Buf,float *A,float *B,int n)
{
	float sum = 0;
	for( int i = 0; i < n; i++){
		sum += A[i]*B[i];
	}
	return sum;
}
///////////////////////////
class CircularBuffer
{
	public:
	float * Buffer;
	void * BufferRAW;
	int buffer_size;
	int pointer;
	void Ctor(int size,Unit *unit){
		buffer_size = size; 
		pointer = 0;
		long vsize = sizeof(float)*nova::vec<float>::size;
		BufferRAW=RTAlloc(unit->mWorld,size * sizeof(float) + vsize);
		Buffer = (float *)((intptr_t)BufferRAW + vsize - ((intptr_t)BufferRAW & (vsize -1)));
		memset(Buffer, 0, size * sizeof(float));
	}
	void push(float a){
		pointer++;
		if(pointer >= buffer_size)
			pointer = 0;
		Buffer[pointer] = a;
	}
	void Dtor(Unit* unit){
		RTFree(unit->mWorld, BufferRAW);
	}
};

class SIMDConvolver
{
	float *Kernel[4];
	void *KernelRAW[4];
	float *Buffer;
	void *BufferRAW;
	int kernel_size;
	Unit * unit;
	public:
	void Ctor(int size,float* data,Unit* unit){
		kernel_size = size;
		unit = unit;
		long vsize = sizeof(float)*nova::vec<float>::size;
		for(int i=0; i<4 ; i++){
			KernelRAW[i]=RTAlloc(unit->mWorld,(kernel_size + i )* sizeof(float) + vsize);
			Kernel[i] = (float *)((intptr_t)KernelRAW[i] + vsize - ((intptr_t)KernelRAW[i] & (vsize -1)));
		}
		BufferRAW=RTAlloc(unit->mWorld,(kernel_size )* sizeof(float) + vsize);
		Buffer = (float *)((intptr_t)BufferRAW + vsize - ((intptr_t)BufferRAW & (vsize -1)));
		for(int g = 0; g < size; g++){
			Kernel[0][g] = data[size - 1 - g];
		}
		for(int i=1; i <4 ;i++){
			memcpy(Kernel[i] + i,Kernel[0],size * sizeof(float));
		}
	}
	void Dtor(Unit* unit){
		for(int i=0; i <4 ;i++){
			RTFree(unit->mWorld, KernelRAW[i]);
		}
		RTFree(unit->mWorld, BufferRAW);
	}
	float ConvolSIMD(float *Al,int pAl,int AlSize)
	{
		assert(AlSize >= kernel_size);
		float sum=0.;
		/*
		int pAl2;
		if (pAl > 0){
			pAl2 = pAl--;
		}else{
			pAl2 = AlSize - 1;
		}
		// return ((intptr_t)(ptr) & (intptr_t)(size * sizeof(float) - 1)) == 0;
		*/
		int pAl2 = (pAl - kernel_size + 1);
		if( pAl2 < 0 ){ pAl2 += AlSize;}
		pAl2 = pAl2%AlSize;
		int vsize = nova::vec<float>::size;
		int alig = pAl2 % vsize;
		sum += convolve_reversed(Buffer, Kernel[alig] + alig, Al + pAl2, AlSize - pAl2,alig);
		alig = (AlSize - pAl2) % vsize;
		alig = (vsize - alig) % vsize;
		sum += convolve_reversed(Buffer, Kernel[alig] + alig + AlSize - pAl2, Al , kernel_size - AlSize + pAl2,alig);

		return sum;
	}
	float ConvolSIMD(CircularBuffer buf)
	{
		int &AlSize = buf.buffer_size;
		assert(AlSize >= kernel_size);
		float sum=0.;

		int pAl2 = (buf.pointer - kernel_size + 1);
		pAl2 = pAl2%AlSize;
		if( pAl2 < 0 ){ pAl2 += AlSize;}
		
		int vsize = nova::vec<float>::size;
		int alig = pAl2 % vsize;
		int howmany = std::min(AlSize - pAl2,kernel_size);
		sum += convolve_reversed(Buffer, Kernel[alig] + alig, buf.Buffer + pAl2, howmany,alig);
		int howmany2 = kernel_size - howmany;
		alig = (howmany) % vsize;
		alig = (vsize - alig) % vsize;
		sum += convolve_reversed(Buffer, Kernel[alig] + alig + howmany, buf.Buffer , howmany2,alig);

		return sum;
	}
	float ConvolSIMD(CircularBuffer buf,int delay)
	{
		int &AlSize = buf.buffer_size;
		assert(AlSize >= kernel_size + delay);
		float sum=0.;

		int pAl2 = (buf.pointer - delay - kernel_size + 1);
		pAl2 = pAl2%AlSize;
		if( pAl2 < 0 ){ pAl2 += AlSize;}
		
		int vsize = nova::vec<float>::size;
		int alig = pAl2 % vsize;
		int howmany = std::min(AlSize - pAl2,kernel_size);
		sum += convolve_reversed(Buffer, Kernel[alig] + alig, buf.Buffer + pAl2, howmany,alig);
		int howmany2 = kernel_size - howmany;
		alig = (howmany) % vsize;
		alig = (vsize - alig) % vsize;
		sum += convolve_reversed(Buffer, Kernel[alig] + alig + howmany, buf.Buffer , howmany2,alig);

		return sum;
	}
};



///////////////////////////////////////////
 float Convol0BAK(float *Reflec,int TamReflec,float *Al,int TamAl,int pAl)
{
		float sum=0.;
		if (pAl<0)
			pAl+=TamAl;

		int pAl2=(pAl + TamReflec - 1 )%TamAl;
		for(int y=0;y<TamReflec;y++){
				sum+=Reflec[y]*Al[pAl2];
				//pAl2 = (pAl2++)%TamAl;
				pAl2++;
				if(pAl2 >=TamAl)
					pAl2 = 0;
		}
		return sum;
}

///////////////////////////////////////////
 float Convol0(float *Reflec[4],int TamReflec,float *Al,int pAl,float *Buffer)
{
		float sum=0.;

		int pAl2;
		if (pAl > 0){
			pAl2 = pAl--;
		}else{
			pAl2 = TamReflec - 1;
		}

		int y = 0;
		for(int x=pAl2; x<TamReflec; y++,x++){
			sum += Reflec[0][y]*Al[x];
		}
		for(int x=0; x<pAl2; y++,x++){
			sum += Reflec[0][y]*Al[x];
		}

		return sum;
}

///////////////////////////////////////////
 float Convol0SIMD(float *Reflec[4],int TamReflec,float *Al,int pAl,float *Buffer1)
{
		float sum=0.;

		int pAl2;
		if (pAl > 0){
			pAl2 = pAl--;
		}else{
			pAl2 = TamReflec - 1;
		}
		// return ((intptr_t)(ptr) & (intptr_t)(size * sizeof(float) - 1)) == 0;

		int vsize = nova::vec<float>::size;
		int alig = pAl2 % vsize;
		sum += convolve_reversed(Buffer1, Reflec[alig] + alig, Al + pAl2, TamReflec - pAl2,alig);
		alig = (TamReflec - pAl2) % vsize;
		alig = (vsize - alig) % vsize;
		sum += convolve_reversed(Buffer1, Reflec[alig] + alig + TamReflec - pAl2, Al , pAl2,alig);

		return sum;
}
//////////////////////////////////////////
//include local buffer test in one place
static SndBuf * ConvGetBuffer(Unit * unit, uint32 bufnum, const char * ugenName, int inNumSamples)
{
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
	SETCALC(*ClearUnitOutputs);
	ClearUnitOutputs(unit, inNumSamples);
	unit->mDone = true;
	return NULL;
}
////////////////////////////////////////////////////
struct AdachiAyers : public Unit
{
	CircularBuffer AlmacenU;
	CircularBuffer AlmacenP;
	SIMDConvolver convolver;
	SIMDConvolver convolver2;
	SIMDConvolver convolver3;
	//float Ref0;
	float Ref0f2;
	float Ref0f1;
	VECTOR2D X,Xold,Xoldold;
	float Uacoust;
	bool OK;
};

extern "C"
{
	void AdachiAyers_Ctor(AdachiAyers* unit);
	void AdachiAyers_next(AdachiAyers *unit, int inNumSamples);
	void AdachiAyers_Dtor(AdachiAyers *unit);
}

#define MAXDELAY 1500
//////////////////////////////////////////
void AdachiAyers_Ctor(AdachiAyers* unit)
{
//	_controlfp(_EM_INEXACT | _EM_UNDERFLOW | _EM_INVALID | _EM_ZERODIVIDE ,_MCW_EM);
//feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

	SndBuf *buf = ConvGetBuffer(unit, (int)ZIN0(3), "AdachiAyers", 1);
	SndBuf *buf2 = ConvGetBuffer(unit, (int)ZIN0(4), "AdachiAyers", 1);
	SndBuf *buf3 = ConvGetBuffer(unit, (int)ZIN0(5), "AdachiAyers", 1);
	
	if (buf == NULL || buf2 == NULL || buf3 == NULL){
		return ;
	}
	LOCK_SNDBUF_SHARED(buf);
	unit->convolver.Ctor(buf->frames - 1,buf->data + 1,(Unit *)unit);
	LOCK_SNDBUF_SHARED(buf2);
	unit->convolver2.Ctor(buf2->frames - 1,buf2->data + 1,(Unit *)unit);
	LOCK_SNDBUF_SHARED(buf3);
	unit->convolver3.Ctor(buf3->frames - 1,buf3->data + 1,(Unit *)unit);
	
	unit->AlmacenU.Ctor(std::max(buf->frames - 1,buf2->frames - 1) + MAXDELAY,unit);
	unit->AlmacenP.Ctor(buf3->frames - 1 + MAXDELAY,unit);
	Print("sizes reflec %d, %d, %d\n",buf->frames,buf2->frames,buf3->frames);
	//float Ref0 = TuboImped[0];
	//unit->Ref0= Ref0;
	//Print("Ref0 %f\n",Ref0);
	unit->Ref0f2=1./(1-buf3->data[0]);
	unit->Ref0f1=(buf->data[0]+buf2->data[0])/(1- buf3->data[0]);
	Print("Ref0 %f, %f\n",unit->Ref0f2,unit->Ref0f1);

	VECTOR2D Xequil={1.e-3,0.0e-3};

	unit->Xoldold.x = unit->Xold.x = unit->X.x = Xequil.x;
	unit->Xoldold.y = unit->Xold.y = unit->X.y = Xequil.y;
	unit->Uacoust=0;

	unit->OK = true;

    SETCALC(AdachiAyers_next);
    AdachiAyers_next(unit,1);
}

void AdachiAyers_Dtor(AdachiAyers* unit)
{
	unit->convolver.Dtor(unit);
	unit->convolver2.Dtor(unit);
	unit->convolver3.Dtor(unit);
	unit->AlmacenU.Dtor(unit);
	unit->AlmacenP.Dtor(unit);
}

void AdachiAyers_next(AdachiAyers *unit, int inNumSamples)
{


	float *out = OUT(0);
	float Flip = IN0(0);
	float P0 = IN0(1);
	float radio = IN0(2);
	float yequil = IN0(6);
	float gate = IN0(7);
	float delay =(int)IN0(8);

	if (gate == 0.0){
		Flip = 0;
		yequil = 0.002;
	}
	if(!unit->OK){
		for (int i=0; i < inNumSamples; ++i)
			out[i] = 0;
		return;
	}


	//Print("Flip %f, P0 %f, radio %f\n",Flip,P0,radio);
	float b=7.e-3;
	float d=2.e-3;
	VECTOR2D Xequil={1.e-3,yequil};//2.e-3};//0.01e-3};//.2e-3};////??
	VECTOR2D Xjoint={0.,1*4.e-3};
	float Qop =3.;//20.;//10.;//5.;//3.;
	float Qcl=.5,Q;
	float m=1.5/(TWOPI*TWOPI*Flip);
	float kop=1.5*Flip;
	float kcl=4.*kop,kefec;
	//float P0=5.e3;//6.e3;// 6.e3;//3.0e3; //.5-6. kPa
	//float &Slip = unit->Slip;
	float Slip;
    float A1 = radio*radio*M_PI;
	float invA1sq = 1/(A1*A1);
	float invA1 = 1/A1;
	//float &Ref0 = unit->Ref0; ;
	float Ref0f2 = unit->Ref0f2;;


	VECTOR2D &X = unit->X;
	VECTOR2D &Xold = unit->Xold;
	VECTOR2D &Xoldold = unit->Xoldold;

	float &Uacoust = unit->Uacoust;
	float Uacoustold;
	float Ug=0,Ulip=0;
	float P1=0.,Plip=0.;
	bool abierto;
	//float Z0=RO_AIRE*C_AIRE/A1;
	float factAreas;
    bool doprint = true;
	float T = SAMPLEDUR;
	float T2 = T*T;

	float Ref0f1divRho = unit->Ref0f1/RO_AIRE;

	for (int i=0; i < inNumSamples; ++i)
	{
		Uacoustold=Uacoust;
		////////nuevas areas
		Slip=sc_max(b*2.*X.y,0.);

		if(Slip>0.)
			abierto=true;
		else
			abierto=false;

		////convolucion
		float convol=0.;
		//convol=Convol0SIMD(Reflec,TamReflec,Almacen,puntAlm,unit->Buffer);//puntAlm+primero-1);
		//convol = unit->convolver.ConvolSIMD(unit->Almacen.Buffer,unit->Almacen.pointer,unit->Almacen.buffer_size);
		
		convol = unit->convolver.ConvolSIMD(unit->AlmacenU);
		convol += unit->convolver2.ConvolSIMD(unit->AlmacenU,delay);
		convol += unit->convolver3.ConvolSIMD(unit->AlmacenP,delay);
		float Ref0f2xconvol = Ref0f2*convol;
		//Ulip///////////////////////
		Ulip=(X.x-Xjoint.x)*(X.y-Xold.y)-((X.y-Xjoint.y)*(X.x-Xold.x));
		Ulip*=b/T;
			//Ulip=(Xjoint.y-X.y)*X.x-((Xjoint.y-Xold.y)*Xold.x);
			//Ulip*=b/T;
			//Ulip=0.;////////////////

		if(abierto){

			/////////ecuaciones para u
			//factAreas=(1./(A1*Slip)-(1./(A1*A1)));
			//factAreas=((A1-Slip)/(A1*A1*Slip));
			factAreas = invA1/Slip -invA1sq;
			//float Rv1=12*MU_AIRE*b*b*d/(Slip*Slip*Slip);	//viscosidad


			double A=.5/(Slip*Slip)-factAreas;
			double B=d/(Slip*T)+Ref0f1divRho ;//+ Rv1;
			double C=-d*Uacoustold/(Slip*T)+Ref0f1divRho*Ulip
				+(Ref0f2xconvol-P0)/RO_AIRE;

			double temporaiz=B*B-4.*A*C;
			if(temporaiz<0.){
				Uacoust = 0;
			}else{
				temporaiz=sqrt(temporaiz);

				double q;
				if(B>=0)
					q=-.5*(B+temporaiz);
				else
					q=-.5*(B-temporaiz);

				float solu1=q/A;
				//solu2=C/q;
				if(solu1>0){
					Uacoust=solu1;
				}else{
					solu1 = C/q;
					if(solu1 >0)
						Uacoust=solu1;
					else
						Uacoust = 0;
				}
				/*
				if((solu1>0 && solu2>0) || (solu1<0 && solu2<0)){
				//Print("iguales en signo %g %g\n",solu1,solu2);
				//if(fabs(solu1-Uacoustold)<fabs(solu2-Uacoustold))
				//	Uacoust=solu1;
				//else
				//	Uacoust=solu2;
					Uacoust = 0;
				//unit->OK = false;
				//return;
				}
				*/
			}

		}else{
			Uacoust=0.;
		}
		Ug = Uacoust + Ulip;
		///calculo P1 prepara convolucion
		P1 = Ref0f2xconvol + unit->Ref0f1*Ug;
		/*
		puntAlm++;
		if(puntAlm >= TamReflec)
			puntAlm = 0;
		Almacen[puntAlm] = P1 + Z0*Ug;
		*/
		//unit->Almacen.push(P1 + Z0*Ug);
		unit->AlmacenU.push(Ug);
		unit->AlmacenP.push(P1);
	 	//presion labios
		if(abierto){
			Plip=P1-Uacoust*Uacoust*factAreas*RO_AIRE;
			Q=Qop;
			kefec=kop;
		}else{
			Plip=P0;//sc_max(P1,P0);
			Q=Qcl;
			kefec=kcl;
		}

	/*
		float distYequil=X.y-Xequil.y;
		float spk=100.e+4;
		float FopNoLinearY=kop*spk*distYequil*distYequil*distYequil;
		float FclNoLinearY=0;
		if(!abierto){
			float spkcl=500.e4;
			FclNoLinearY=X.y*X.y*X.y*kcl*spkcl;
		}
		*/

		////////////////////////sistema muelles

		Xoldold.x=Xold.x;
		Xold.x=X.x;

		Xoldold.y=Xold.y;
		Xold.y=X.y;

		float ftrans = b * (P0 - P1);
		float facVelocidad = sqrt(m*kefec)/(Q*T);

		X.y = -T2/m *
		(facVelocidad*(Xold.y - Xoldold.y) + kefec*(Xold.y -Xequil.y) - 2*ftrans*(Xold.x -Xjoint.x)
			- 2* b*d*Plip) + 2*Xold.y - Xoldold.y;
		X.x = -T2/m *
		(facVelocidad*(Xold.x - Xoldold.x) + kefec*(Xold.x -Xequil.x) + 2*ftrans*(Xold.y -Xjoint.y)
			) + 2*Xold.x - Xoldold.x;

		//if((Xold.y -Xjoint.y) > 0 and doprint){
		//	Print("(Xold.y -Xjoint.y) mayor 0 \n");
		//	doprint = false;
		//}

		//Print("%f %f ->",Ulip,Uacoust);
		/////////output
		out[i] = Ug;//Plip;//Ug;

	}
	//assert(Uacoust == unit->Uacoust);
}

/////////////////////////////////////////
PluginLoad(AdachiAyers)
{
	ft = inTable;
	DefineDtorUnit(AdachiAyers);
}
