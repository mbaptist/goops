


// namespace goops
// {

// #include <fftw3.h>
// 
// #include <cmath>

using namespace std;

//Generic Ctor
template <class TypeOut,class TypeIn>
PlanBase<TypeOut,TypeIn>::PlanBase()
{
	realdata=0;
	fourierdata=0;
	plans_created=0;
}

//Generic Dtor
template <class TypeOut,class TypeIn>
PlanBase<TypeOut,TypeIn>::~PlanBase()
{
	destroy_plans();
	realdata=0;
	fourierdata=0;
}

//Generic create plan
template <class TypeOut,class TypeIn>
void PlanBase<TypeOut,TypeIn>::create_plans()
{
	do_create_plan();
	plans_created=1;
}

//Generic destroy plans
template <class TypeOut,class TypeIn>
void PlanBase<TypeOut,TypeIn>::destroy_plans()
{
	if(plans_created)
	{
		fftw_destroy_plan(direct_plan);
		fftw_destroy_plan(inverse_plan);
	}
	plans_created=0;	
}

template <class TypeOut,class TypeIn>
void PlanBase<TypeOut,TypeIn>::execute()
{
	fftw_execute(direct_plan);
}


//Generic switch data
template <class TypeOut,class TypeIn>
void PlanBase<TypeOut,TypeIn>::switch_data(TypeOut & outfield,const TypeIn & infield)
{
	TypeOutElement * outfield_ptr=reinterpret_cast<TypeOutElement *>(outfield.data());
	TypeInElement * infield_ptr=reinterpret_cast<TypeInElement *>(infield.data());
	if(datain!=infield_ptr||dataout!=outfield_ptr)
	{
		destroy_plan();
		datain=infield_ptr;
		dataout=outfield_ptr;
		realsize=infield.size();
		fouriersize=outfield.size();
		create_plan();
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
	void Plan<cat::array<RS,D>,cat::array<CS,D> >::do_create_plans()
{
	direct_plan=fftw_plan_dft_r2c(D,&realsize,realdata,fourierdata,FFTW_ESTIMATE);
	inverse_plan=fftw_plan_dft_c2r(D,&fouriersize,fourierdata,realdata,FFTW_ESTIMATE);
}




template <>
	template <int D>
	Plan<cat::array<RS,D>,cat::array<RS,D> >::Plan(const string & subtype):
	PlanBase<cat::array<RS,D>,cat::array<RS,D> >()
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
	void Plan<cat::array<RS,D>,cat::array<RS,D> >::do_create_plans()
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
	void Plan<cat::array<RS,D>,cat::array<RS,D> >::direct_transform(cat::array<RS,D> & fourierfield,const cat::array<RS,D> & realfield)
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
	void Plan<cat::array<RS,D>,cat::array<RS,D> >::inverse_transform(cat::array<RS,D> & realfield,const cat::array<RS,D> & fourierfield)
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
