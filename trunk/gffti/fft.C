
#include <fftw3.h>

namespace goops
{

//Generic Ctor
template <class RealType,class FourierType>
FFT<RealType,FourierType>::FFT():
realdata(0),
fourierdata(0)
{
}

//Generic Dtor
template <class RealType,class FourierType>
FFT<RealType,FourierType>::~FFT()
{
	realdata=0;
	fourierdata=0;
}

//Generic direct transform
template <class RealType,class FourierType>
void FFT<RealType,FourierType>::direct_transform(FourierType & fourierfield,RealType & realfield)
{
	if(realdata!=realfield.data()||fourierdata!=fourierfield.data())
	{
		destroy_plans();
		realdata=realfield.data();
		fourierdata=fourierfield.data();
		size=fourierfield.size();
		create_plans();
	}
	fftw_execute(direct_plan);
	fourierfield/=size;
}

//Generic inverse transform
template <class RealType,class FourierType>
void FFT<RealType,FourierType>::inverse_transform(RealType & realfield,FourierType & fourierfield)
{
	if(realdata!=realfield.data()||fourierdata!=fourierfield.data())
	{
		destroy_plans();
		realdata=realfield.data();
		fourierdata=fourierfield.data();
		size=fourierfield.size();
		create_plans();
	}
	fftw_execute(inverse_plan);
	realfield/=size;
}

//Generic dstroy plan
template <class RealType,class FourierType>
	void FFT<RealType,FourierType>::destroy_plans()
{
	fftw_destroy_plan(direct_plan);
	fftw_destroy_plan(inverse_plan);
}



//Partial specialization




template <int D>
	void FFT<cat::array<RS,D>,cat::array<CS,D> >::create_plans()
{
	direct_plan=fftw_plan_dft_r2c_1d(&size,realdata,fourierdata,FFTW_ESTIMATE);
	inverse_plan=fftw_plan_dft_c2r_1d(&size,fourierdata,realdata,FFTW_ESTIMATE);
}




