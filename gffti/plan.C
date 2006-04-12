


// namespace goops
// {

// #include <fftw3.h>
// 
// #include <cmath>

using namespace std;

//Generic Ctor
template <class RealType,class FourierType>
FFT_BASE<RealType,FourierType>::FFT_BASE()
{
	realdata=0;
	fourierdata=0;
	plans_created=0;
}

//Generic Dtor
template <class RealType,class FourierType>
FFT_BASE<RealType,FourierType>::~FFT_BASE()
{
	destroy_plans();
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
void FFT_BASE<RealType,FourierType>::direct_transform(FourierType & fourierfield, const RealType & realfield)
{
	switch_data(*const_cast<RealType*>(&realfield),fourierfield);
	fftw_execute(direct_plan);
	fourierfield/=realsize;
}

//Generic inverse transform
template <class RealType,class FourierType>
void FFT_BASE<RealType,FourierType>::inverse_transform(RealType & realfield,const FourierType & fourierfield)
{
	switch_data(realfield,*const_cast<FourierType*>(&fourierfield));
	fftw_execute(inverse_plan);
	realfield/=fouriersize;
}

//Generic switch data
template <class RealType,class FourierType>
void FFT_BASE<RealType,FourierType>::switch_data(RealType & realfield,FourierType & fourierfield)
{
	RDT * realfield_ptr=reinterpret_cast<RDT*>(realfield.data());
	FDT * fourierfield_ptr=reinterpret_cast<FDT*>(fourierfield.data());
	if(realdata!=realfield_ptr||fourierdata!=fourierfield_ptr)
	{
		destroy_plans();
		realdata=realfield_ptr;
		fourierdata=fourierfield_ptr;
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
	if (r2r_kind_direct==FFTW_RODFT00)
	{
		realsize-=2;
		fouriersize-=2;
		direct_plan=fftw_plan_r2r(D,&realsize,realdata+1,fourierdata+1,&r2r_kind_direct,FFTW_ESTIMATE);
		inverse_plan=fftw_plan_r2r(D,&fouriersize,fourierdata+1,realdata+1,&r2r_kind_inverse,FFTW_ESTIMATE);
	}
	else
	{
	direct_plan=fftw_plan_r2r(D,&realsize,realdata,fourierdata,&r2r_kind_direct,FFTW_ESTIMATE);
	inverse_plan=fftw_plan_r2r(D,&fouriersize,fourierdata,realdata,&r2r_kind_inverse,FFTW_ESTIMATE);
	}
}

template <>
	template <int D>
	void FFT<cat::array<RS,D>,cat::array<RS,D> >::direct_transform(cat::array<RS,D> & fourierfield,const cat::array<RS,D> & realfield)
{
	switch_data(*const_cast<cat::array<RS,D>*>(&realfield),fourierfield);
	fftw_execute(direct_plan);
	if(r2r_kind_direct==FFTW_REDFT00)
		fourierfield/=(realsize-1);
	else
		fourierfield/=(realsize-1);
	fourierfield.data()[0]*=.5;
}
template <>
	template <int D>
	void FFT<cat::array<RS,D>,cat::array<RS,D> >::inverse_transform(cat::array<RS,D> & realfield,const cat::array<RS,D> & fourierfield)
{
	switch_data(realfield,*const_cast<cat::array<RS,D>*>(&fourierfield));
	fftw_execute(direct_plan);
	fourierfield*=.5;
	if (r2r_kind_direct==FFTW_RODFT00)
	{
	fourierfield.data()[0]=0;
	fourierfield.data()[fouriersize+1]=0;
	}
}
	

//}
