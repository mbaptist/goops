
#include <fftw3.h>

namespace goops
{

//Generic Ctor
template <class RealType,class FourierType>
FFT_BASE<RealType,FourierType>::FFT_BASE()
{
	realdata=0;
	fourierdata=0;
}

//Generic Dtor
template <class RealType,class FourierType>
FFT_BASE<RealType,FourierType>::~FFT_BASE()
{
	realdata=0;
	fourierdata=0;
}

//Generic create plans
template <class RealType,class FourierType>
	void FFT_BASE<RealType,FourierType>::create_plans()
{
	do_create_plans();
	plans_created=1;
}

//Generic destroy plans
template <class RealType,class FourierType>
	void FFT_BASE<RealType,FourierType>::destroy_plans()
{
	if(plans_created)
	{
		fftw_destroy_plan(direct_plan);
		fftw_destroy_plan(inverse_plan);
	}
	plans_created=0;	
}

//Generic direct transform
template <class RealType,class FourierType>
void FFT_BASE<RealType,FourierType>::direct_transform(FourierType & fourierfield,RealType & realfield)
{
	switch_data(realfield,fourierfield);
	fftw_execute(direct_plan);
	fourierfield/=realsize;
}

//Generic inverse transform
template <class RealType,class FourierType>
void FFT_BASE<RealType,FourierType>::inverse_transform(RealType & realfield,FourierType & fourierfield)
{
	switch_data(realfield,fourierfield);
	fftw_execute(inverse_plan);
	realfield/=fouriersize;
}

//Generic inverse transform
template <class RealType,class FourierType>
	void FFT_BASE<RealType,FourierType>::switch_data(RealType & realfield,FourierType & fourierfield)
{
	if(realdata!=realfield.data()||fourierdata!=reinterpret_cast<FDT*>(fourierfield.data()))
	{
		destroy_plans();
		realdata=realfield.data();
		fourierdata=reinterpret_cast<FDT*>(fourierfield.data());
		realsize=realfield.size();
		fouriersize=fourierfield.size();
		create_plans();
	}
}

//Specialization


template <>
	template <int D>
	void FFT<cat::array<CS,D>,cat::array<CS,D> >::do_create_plans()
{
	direct_plan=fftw_plan_dft(D,&realsize,realdata,fourierdata,FFTW_FORWARD,FFTW_ESTIMATE);
	inverse_plan=fftw_plan_dft(D,&fouriersize,fourierdata,realdata,FFTW_FORWARD,FFTW_ESTIMATE);
}


template <>
	template <int D>
	void FFT<cat::array<RS,D>,cat::array<CS,D> >::do_create_plans()
{
	direct_plan=fftw_plan_dft_r2c(D,&realsize,realdata,fourierdata,FFTW_ESTIMATE);
	inverse_plan=fftw_plan_dft_c2r(D,&fouriersize,fourierdata,realdata,FFTW_ESTIMATE);
}




template <>
	template <int D>
	FFT<cat::array<RS,D>,cat::array<RS,D> >::FFT(const string & subtype):
	FFT_BASE<cat::array<RS,D>,cat::array<RS,D> >()
{
	if(subtype=="sin")
	{
		r2r_kind_direct=FFTW_RODFT00;
		r2r_kind_inverse=FFTW_RODFT00;
		//norm_factor_direct=(fouriersize+1)
	}
	else if(subtype=="cos")
	{
		r2r_kind_direct=FFTW_REDFT00;
		r2r_kind_inverse=FFTW_REDFT00;
	}
	else
	{
		cout << "Undefined subtype! Aborting..." << endl;
	}
}

template <>
	template <int D>
	void FFT<cat::array<RS,D>,cat::array<RS,D> >::do_create_plans()
{
	direct_plan=fftw_plan_r2r(D,&realsize,realdata-1,fourierdata,&r2r_kind_direct,FFTW_ESTIMATE);
	inverse_plan=fftw_plan_r2r(D,&fouriersize,fourierdata,realdata+1,&r2r_kind_inverse,FFTW_ESTIMATE);
}

// template <>
// 	template <int D>
// 	void FFT<cat::array<RS,D>,cat::array<RS,D> >::direct_transform(FourierType & fourierfield,RealType & realfield)
// {
// 	switch_data(realfield,fourierfield);
// 	fftw_execute(direct_plan);
// 	fourierfield/=realsize;
// }
// template <>
// 	template <int D>
// 	void FFT<cat::array<RS,D>,cat::array<RS,D> >::inverse_transform(RealType & realfield,FourierType & fourierfield)
// {
// 	switch_data(realfield,fourierfield);
// 	fftw_execute(direct_plan);
// 	fourierfield/=realsize;
// }
	

}
