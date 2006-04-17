
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
	datain=0;
	dataout=0;
}

//Generic create plan
template <class TypeOut,class TypeIn>
void PlanBase<TypeOut,TypeIn>::create_plan()
{
	do_create_plan();
	plan_exists=1;
}

//Generic destroy plan
template <class TypeOut,class TypeIn>
void PlanBase<TypeOut,TypeIn>::destroy_plan()
{
	if(plan_exists)
		fftw_destroy_plan(plan);
	plan_exists=0;
}

template <class TypeOut,class TypeIn>
void PlanBase<TypeOut,TypeIn>::execute()
{
	fftw_execute(plan);
}

//Generic switch data
template <class TypeOut,class TypeIn>
void PlanBase<TypeOut,TypeIn>::switch_data(TypeOut & fieldout,TypeIn & fieldin)
{
	fftTypeOut * fieldout_ptr=reinterpret_cast<fftTypeOut *>(fieldout.data());
	fftTypeIn * fieldin_ptr=reinterpret_cast<fftTypeIn *>(fieldin.data());
	if(dataout!=fieldout_ptr||datain!=fieldin_ptr)
	{
		destroy_plan();
		dataout=fieldout_ptr;
		dataoutshape=fieldout.shape().data();
		dataoutsize=evalsize(dataoutshape,fftTypes<TypeOut>::Rank);
		datain=fieldin_ptr;
		datainshape=fieldin.shape().data();
		datainsize=evalsize(datainshape,fftTypes<TypeOut>::Rank);
		create_plan();
	}
}


template <class TypeOut,class TypeIn>
int PlanBase<TypeOut,TypeIn>::evalsize(int * shape,int rank)
{
	int size_=1;
	for( int i=0;i<rank;++i)
		size_*=shape[rank];
	return size_;
}

//Class Plan implementation

//Specialization

//Complex to Complex transforms
//Ctor
template <>
template <int D>
Plan<cat::array<CS,D>,const cat::array<CS,D> >::Plan(const string & direction):
PlanBase<cat::array<CS,D>,const cat::array<CS,D> >()
{
	if (direction=="direct"||direction=="forward")
		fftdirection=FFTW_FORWARD;
	else if (direction=="inverse"||direction=="backward")
		fftdirection=FFTW_BACKWARD;
}
//Plan creation
template <>
template <int D>
void Plan<cat::array<CS,D>,const cat::array<CS,D> >::do_create_plan()
{
	plan=fftw_plan_dft(D,datainshape,datain,dataout,fftdirection,FFTW_ESTIMATE);
}

//Real to Complex transforms
template <>
template <int D>
void Plan<cat::array<CS,D>,const cat::array<RS,D> >::do_create_plan()
{
	fftw_plan_dft_r2c(D,datainshape,datain,dataout,FFTW_ESTIMATE);
}
//Complex to Real transforms
template <>
template <int D>
void Plan<cat::array<RS,D>,const cat::array<CS,D> >::do_create_plan()
{
	plan=fftw_plan_dft_c2r(D,datainshape,datain,dataout,FFTW_ESTIMATE);
}


//1D Real to Real
//template <>
Plan<cat::array<RS,1>,const cat::array<RS,1> >::Plan(const string & subtype__,const string & direction__):
PlanBase<cat::array<RS,1>,const cat::array<RS,1> >(),
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

//template <>
void Plan<cat::array<RS,1>,const cat::array<RS,1> >::do_create_plan()
{
	plan=fftw_plan_r2r_1d(r2r_size,r2r_datain,r2r_dataout,r2r_kind,FFTW_ESTIMATE);
}

//1D Real to Real switch data
//template <>
void Plan<cat::array<RS,1>,const cat::array<RS,1> >::switch_data(cat::array<RS,1> & fieldout,const cat::array<RS,1> & fieldin)
{
	fftTypeOut * fieldout_ptr=reinterpret_cast<fftTypeOut *>(const_cast<fftTypeOut *>(fieldout.data()));
	fftTypeIn * fieldin_ptr=reinterpret_cast<fftTypeIn *>(const_cast<fftTypeIn *>(fieldin.data()));
	if(dataout!=fieldout_ptr||datain!=fieldin_ptr)
	{
		destroy_plan();
		dataout=fieldout_ptr;
		dataoutshape=fieldout.shape().data();
		dataoutsize=evalsize(dataoutshape,1);
		datain=fieldin_ptr;
		datainshape=const_cast<int *>(fieldin.shape().data());
		datainsize=evalsize(datainshape,1);
		if(subtype=="sin")
		{
			r2r_datain=datain+1;
			r2r_dataout=dataout+1;
			r2r_size=datainsize-2;
			if (direction=="direct"||direction=="forward")
			{
				transparcel=0	;			
				normfactor=1./(datainsize-1);
				normfactor_zero=1.;
			}
			else if(direction=="inverse"||direction=="backward")
			{
				transparcel=0	;			
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
				normfactor_zero=.5;
			}
			else if(direction=="inverse"||direction=="backward")
			{			
				normfactor=.5;
				normfactor_zero=1.;
			}
		}
		create_plan();
	}
}

	//1D Real to Real normalise
//	template <>
	void Plan<cat::array<RS,1>,const cat::array<RS,1> >::normalise()
	{
		for(int i=0;i<dataoutsize;++i)
		{
			dataout[i]+=transparcel;			
			dataout[i]*=normfactor;
		}
		dataout[0]*=normfactor_zero;
		if (subtype=="sin")
		{
			dataout[0]=0.;
			dataout[dataoutsize-1]=0.;
		}
	}

//}
