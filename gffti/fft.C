


// namespace goops
// {

#include <fftw3.h>
// 
// #include <cmath>

#include "plan.h"

#include <cat.h>

using namespace std;

using namespace cat;

//Generic Ctors
template <class RealType,class FourierType>
FFT<RealType,FourierType>::FFT():
direct_plan(),
inverse_plan()
{
}

//Generic Dtor
template <class RealType,class FourierType>
FFT<RealType,FourierType>::~FFT()
{
}

//Generic direct transform
template <class RealType,class FourierType>
void FFT<RealType,FourierType>::direct_transform(FourierType & fourierfield, const RealType & realfield)
{
	//	cout << field << endl;
	direct_plan.switch_data(fourierfield,realfield);
	direct_plan.execute();
	fourierfield/=realfield.size();
}

//Generic inverse transform
template <class RealType,class FourierType>
void FFT<RealType,FourierType>::inverse_transform(RealType & realfield,const FourierType & fourierfield)
{
	inverse_plan.switch_data(realfield,fourierfield);
	inverse_plan.execute();
}


///////////////////////////////

//Complex to Complex Ctor
template <>
template<int D>
FFT<cat::array<CS,D>,cat::array<CS,D> >::FFT():
direct_plan(FFTW_FORWARD),
inverse_plan(FFTW_BACKWARD)
{
}

template <>
template<int D>
FFT<cat::array<CS,D>,cat::array<CS,D> >::~FFT()
{
}

//1D Real to Real Ctor
//template <>
FFT<cat::array<RS,1>,cat::array<RS,1> >::FFT(const string & subtype__):
subtype(subtype__),
direct_plan(r2r_direct_kind()),
inverse_plan(r2r_inverse_kind())
{
}

FFT<cat::array<RS,1>,cat::array<RS,1> >::~FFT()
{
}

//1D Real to Real direct transform
//template <>
void FFT<cat::array<RS,1>,cat::array<RS,1> >::direct_transform(cat::array<RS,1> & fourierfield, const cat::array<RS,1> & realfield)
{
	if(subtype=="sin")
	{
		cat::array<RS,1> ff=cat::array<RS,1>(cat::tvector<int,1>(fourierfield.size()-2),fourierfield.data()+1);
		const cat::array<RS,1> rf=cat::array<RS,1>(cat::tvector<int,1>(realfield.size()-2),const_cast<RS *>(realfield.data())+1);
		direct_plan.switch_data(ff,rf);
		fourierfield(0)=0;
		fourierfield(fourierfield.size()-1)=0;
	}
	else
	{
		inverse_plan.switch_data(fourierfield,realfield);
		fourierfield(0)*=.5;
	}
	fourierfield/=(realfield.size()-1);
}

//1D Real to Real inverse transform
//template <>
void FFT<cat::array<RS,1>,cat::array<RS,1> >::inverse_transform(cat::array<RS,1> & realfield,const cat::array<RS,1> & fourierfield)
{
	if(subtype=="sin")
	{
		cat::array<RS,1> rf=cat::array<RS,1>(cat::tvector<int,1>(realfield.size()-2),realfield.data()+1);
		const cat::array<RS,1> ff=cat::array<RS,1>(cat::tvector<int,1>(fourierfield.size()-2),const_cast<RS *>(fourierfield.data())+1);
		inverse_plan.switch_data(rf,ff);
		realfield(0)=0;
		realfield(realfield.size()-1)=0;
	}
	else
	{
		inverse_plan.switch_data(realfield,fourierfield);
		RS av=realfield(0);
		realfield+=av;
	}
	realfield*=.5;
}

//template <>
fftw_r2r_kind FFT<cat::array<RS,1>,cat::array<RS,1> >::r2r_direct_kind()
{
	if (subtype=="sin")
		return FFTW_RODFT00;
	else
		return FFTW_REDFT00;
}

//template <>
fftw_r2r_kind FFT<cat::array<RS,1>,cat::array<RS,1> >::r2r_inverse_kind()
{
	if (subtype=="sin")
		return FFTW_RODFT00;
	else
		return FFTW_REDFT00;
}
