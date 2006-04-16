


// namespace goops
// {

// #include <fftw3.h>
// 
// #include <cmath>

using namespace std;


//Class PlanBase implementation

//Ctor
template <class TypeOut,class TypeIn>
PlanBase<TypeOut,TypeIn>::PlanBase():
dataout(0),
datain(0),
plan_exists(0)
{
}

//Generic Dtor
template <class TypeOut,class TypeIn>
PlanBase<TypeOut,TypeIn>::~PlanBase()
{
	destroy_plan();
	realdata=0;
	fourierdata=0;
}

//Generic create plan
template <class TypeOut,class TypeIn>
void PlanBase<TypeOut,TypeIn>::create_plan()
{
	do_create_plan();
	plans_created=1;
}

//Generic destroy plans
template <class TypeOut,class TypeIn>
void PlanBase<TypeOut,TypeIn>::destroy_plan()
{
	if(plan_exists)
	{
		fftw_destroy_plan(direct_plan);
		fftw_destroy_plan(inverse_plan);
	}
	plan_exists=0;
}

template <class TypeOut,class TypeIn>
void PlanBase<TypeOut,TypeIn>::execute()
{
	fftw_execute(direct_plan);
}

//Generic switch data
template <class TypeOut,class TypeIn>
void PlanBase<TypeOut,TypeIn>::switch_data(TypeOut & fieldout,const TypeIn & fieldin)
{
	fftTypeOut * fieldout_ptr=reinterpret_cast<fftTypeOut *>(fieldout.data());
	fftTypeIn * fieldin_ptr=reinterpret_cast<fftTypeIn *>(fieldin.data());
	if(dataout!=fieldout_ptr||datain!=fieldin_ptr)
	{
		destroy_plan();
		dataout=fieldout_ptr;
		dataoutshape=fieldout.shape().data();
		dataoutsize=evalsize(dataoutshape,);
		datain=fieldin_ptr;
		datainshape=fieldin.shape().data();
		datainsize=evalsize(datainshape,);
		create_plan();
	}
}


//Class Plan implementation

//Specialisation

//Complex to Complex transforms
//Ctor
template <>
template <int D>
void Plan<cat::array<CS,D>,cat::array<CS,D> >::Plan(string & direction):
PlanBase()
{
	if (direction=="direct"||direction=="forward")
		fftdirection=FFTW_FORWARD;
	else if (direction=="inverse"||direction=="backward")
		fftdirection=FFTW_BACKWARD;
}
//Plan creation
template <>
template <int D>
void Plan<cat::array<CS,D>,cat::array<CS,D> >::do_create_plan()
{
	plan=fftw_plan_dft(D,datainshape,datain,dataout,fftdirection,FFTW_ESTIMATE);
}

//Real to Complex transforms
template <>
template <int D>
void Plan<cat::array<CS,D>,cat::array<RS,D> >::do_create_plan()
{
	fftw_plan_dft_r2c(D,datainshape,datain,dataout,FFTW_ESTIMATE);
}
//Complex to Real transforms
template <>
	template <int D>
	void Plan<cat::array<RS,D>,cat::array<CS,D> >::do_create_plans()
{
	direct_plan=fftw_plan_dft_r2c(D,&realsize,realdata,fourierdata,FFTW_ESTIMATE);
	plan=fftw_plan_dft_c2r(D,datainshape,datain,dataout,FFTW_ESTIMATE);
}



//1D Real to Real
template <>
Plan<cat::array<RS,1>,cat::array<RS,1> >::Plan(const string & subtype__,const string & direction__):
PlanBase<cat::array<RS,1>,cat::array<RS,1> >(),
subtype(subtype__),
direction(direction__)
{
	if(subtype=="sin")
	{
		if(direction=="direct"||direction=="forward")
			r2r_kind=FFTW_RODFT00;
		else if(direction=="inverse"||direction=="backward")
			r2r_kind=FFTW_RODFT00;
	}
	else if(subtype=="cos")
	{
		if(direction=="direct"||direction=="forward")
			r2r_kind=FFTW_RODFT00;
		else if(direction=="inverse"||direction=="backward")
			r2r_kind=FFTW_RODFT00;
	}
	else
	{
		cout << "Undefined subtype! Aborting..." << endl;
	}
}

template <>
void Plan<cat::array<RS,1>,cat::array<RS,1> >::do_create_plan()
{
		direct_plan=fftw_plan_r2r_1d(r2r_size,fftdatain,fftdataout,r2r_kind,FFTW_ESTIMATE);
}

//1D Real to Real switch data
template <>
void Plan<cat::array<RS,1>,cat::array<RS,1> >::switch_data(cat::array<RS,1> & fieldout,const cat::array<RS,1> & fieldin)
{
	fftTypeOut * fieldout_ptr=reinterpret_cast<fftTypeOut *>(fieldout.data());
	fftTypeIn * fieldin_ptr=reinterpret_cast<fftTypeIn *>(fieldin.data());
	if(dataout!=fieldout_ptr||datain!=fieldin_ptr)
	{
		destroy_plan();
		dataout=fieldout_ptr;
		dataoutshape=fieldout.shape().data();
		dataoutsize=evalsize(dataoutshape,);
		datain=fieldin_ptr;
		datainshape=fieldin.shape().data();
		datainsize=evalsize(datainshape,);
		if(subtype=="sin")
		{
			r2r_datain=datain+1;
			r2r_dataout=dataout+1;
			r2r_size=datainsize-2;
			if (direction=="direct"||direction=="forward")
			{
				normfactor=datainsize-1;
				normfactor_zero=.5;
			}
			else if(direction=="inverse"||direction=="backward")
			{
				normfactor=.5;
				normfactor_zero=1.;
			}
		}
		else if(subtype=="cos")
		{
			r2r_datain=datain;
			r2r_dataout=dataout;
			r2r_size=datainsize;
			if (direction=="direct"||direction=="forward")
			{
				normfactor=1./(datainsize-1);
				normfactor_zero=1.;
			}
			else if(direction=="inverse"||direction=="backward")
			{
				normfactor=.5;
				normfactor_zero=1.;
			}
		}
		create_plan();
	}

	//1D Real to Real normalise
	template <>
		void Plan<cat::array<RS,1>,cat::array<RS,1> >::normalise()
	{
		for(int i=0;i<outdatasize;++i)
			outdata[i]*=normfactor;			
		outdata[0]*=normfactor_zero;
		if (direction=="inverse"||direction=="backward")
		{
			outdata[0]=0.;
			outdata[outdatasize-1]=0.;
		}
	}

//}
