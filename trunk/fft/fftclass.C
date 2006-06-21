


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
	//direct_plan.switch_data(fourierfield,realfield);
	//direct_plan.execute();
	direct_plan.execute(fourierfield,realfield);
	fourierfield/=realfield.size();
}

//Generic inverse transform
template <class RealType,class FourierType>
void FFT<RealType,FourierType>::inverse_transform(RealType & realfield,const FourierType & fourierfield)
{
	//inverse_plan.switch_data(realfield,fourierfield);
	//inverse_plan.execute();
	inverse_plan.execute(realfield,fourierfield);
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
template <>
template <int D>
FFT<cat::array<RS,D>,cat::array<RS,D> >::FFT(const string & subtype__):
subtype(subtype__),
direct_plan(r2r_direct_kind(subtype__)),
inverse_plan(r2r_inverse_kind(subtype__))
{
	if (subtype=="sin")
	{
		direct_transform_imp=&FFT<cat::array<RS,D>,cat::array<RS,D> >::sin_direct_transform;
		inverse_transform_imp=&FFT<cat::array<RS,D>,cat::array<RS,D> >::sin_inverse_transform;
	}
	else
	{
		direct_transform_imp=&FFT<cat::array<RS,D>,cat::array<RS,D> >::cos_direct_transform;
		inverse_transform_imp=&FFT<cat::array<RS,D>,cat::array<RS,D> >::cos_inverse_transform;
	}
}

template <>
template <int D>
FFT<cat::array<RS,D>,cat::array<RS,D> >::~FFT()
{
}

//1D Real to Real direct transform
template <>
template <int D>
void FFT<cat::array<RS,D>,cat::array<RS,D> >::direct_transform(cat::array<RS,D> & fourierfield, const cat::array<RS,D> & realfield)
{
	(this->*direct_transform_imp)(fourierfield,realfield);
}

//1D Real to Real inverse transform
template <>
template <int D>
void FFT<cat::array<RS,D>,cat::array<RS,D> >::inverse_transform(cat::array<RS,D> & realfield,const cat::array<RS,D> & fourierfield)
{
	(this->*inverse_transform_imp)(realfield,fourierfield);
}

//1D Real to Real direct sine transform 
template <>
template <int D>
void FFT<cat::array<RS,D>,cat::array<RS,D> >::sin_direct_transform(cat::array<RS,D> & fourierfield, const cat::array<RS,D> & realfield)
{
		//		cat::array<RS,D> ff=cat::array<RS,D>(cat::tvector<int,D>(fourierfield.size()-2),fourierfield.data()+1);
		//  const cat::array<RS,D> rf=cat::array<RS,D>(cat::tvector<int,D>(realfield.size()-2),const_cast<RS *>(realfield.data())+1);
		
		size_t len=fourierfield.length()-2;
		cat::array<RS,D> ff=cat::array<RS,D>(cat::tvector<int,D>(fourierfield.size()-2),fourierfield.ordering(),fourierfield.stride(),len,fourierfield.data()+1);
		
// 		array<T,D>::array(tvector<int,D> & shape__,
// 		                  tvector<int,D> & ordering__,
// 		                  tvector<int,D> & stride__,
// 		                  size_t & length__,
// 		                  T * data__):
		
		len=realfield.length()-2;
		const cat::tvector<int,D> sha(realfield.size()-2);
		const cat::tvector<int,D> ord(realfield.ordering());
		cat::tvector<int,D> str(realfield.stride());
		const RS * p = realfield.data()+1;
		const cat::array<RS,D> rf(sha,ord,str,len,p);
		
	//direct_plan.switch_data(ff,rf);
	//direct_plan.execute();
	direct_plan.execute(ff,rf);
	
	fourierfield(0)=0;
	fourierfield*=1./(realfield.size()-1);
	fourierfield(fourierfield.size()-1)=0;
}

//1D Real to Real inverse sine transform
template <>
template <int D>
void FFT<cat::array<RS,D>,cat::array<RS,D> >::sin_inverse_transform(cat::array<RS,D> & realfield,const cat::array<RS,D> & fourierfield)
{
	//cout << "SINE INVERSE TRANSFORM" << endl;
		//cat::array<RS,D> rf=cat::array<RS,D>(cat::tvector<int,D>(realfield.size()-2),realfield.data()+1);
		size_t len=realfield.length()-2;
			cat::array<RS,D> rf(cat::tvector<int,D>(realfield.size()-2),
			                    realfield.ordering(),
			                    realfield.stride(),
			                    len,
			                    realfield.data()+1);
		
// 		array<T,D>::array(tvector<int,D> & shape__,
// 		                  tvector<int,D> & ordering__,
// 		                  tvector<int,D> & stride__,
// 		                  size_t & length__,
// 		                  T * data__):

		//const cat::array<RS,D> ff=cat::array<RS,D>(cat::tvector<int,D>(fourierfield.size()-2),const_cast<RS *>(fourierfield.data())+1);
		
		len=fourierfield.length()-2;
		const cat::tvector<int,D> sha(fourierfield.size()-2);
		const cat::tvector<int,D> ord(fourierfield.ordering());
		cat::tvector<int,D> str(fourierfield.stride());
		const RS * p = fourierfield.data()+1;
		const cat::array<RS,D> ff(sha,ord,str,len,p);
		
		//cout << rf << "\n\n" << ff << endl;
	//inverse_plan.switch_data(rf,ff);
	//inverse_plan.execute();
	inverse_plan.execute(rf,ff);
	
		//cout << rf << "\n\n" << ff << endl;
		//exit(0);
		
		//cout << realfield << endl;
	realfield(0)=0;
	realfield*=.5;
	realfield(realfield.size()-1)=0;
}
//1D Real to Real direct co-sine transform 
template <>
template <int D>
void FFT<cat::array<RS,D>,cat::array<RS,D> >::cos_direct_transform(cat::array<RS,D> & fourierfield, const cat::array<RS,D> & realfield)
{
	//cout << "COS direct transform" << endl;
	//direct_plan.switch_data(fourierfield,realfield);
	//direct_plan.execute();
	direct_plan.execute(fourierfield,realfield);
	fourierfield(0)*=.5;
	fourierfield*=1./(realfield.size()-1);
}

//1D Real to Real inverse co-sine transform
template <>
template <int D>
void FFT<cat::array<RS,D>,cat::array<RS,D> >::cos_inverse_transform(cat::array<RS,D> & realfield,const cat::array<RS,D> & fourierfield)
{
	// cout << "COS inverse transform" << endl;
// 		inverse_plan.switch_data(realfield,fourierfield);
// 		inverse_plan.execute();
	inverse_plan.execute(realfield,fourierfield);
	
	RS f0=fourierfield(0);
	RS fnm1=fourierfield(fourierfield.size()-1);
	typename cat::array<RS,D>::iterator iter(realfield);
	int sign=1;
	for(iter=realfield.begin();iter!=realfield.end();++iter)
	{
		(*iter)+=f0+sign*fnm1;
		(*iter)*=.5;
		sign*=-1;
	}
// 	realfield+=f0;
// 	realfield*=.5;
	
}

template <>
template <int D>
fftw_r2r_kind FFT<cat::array<RS,D>,cat::array<RS,D> >::r2r_direct_kind(const string & subtype__)
{
	if (subtype__=="sin")
		return FFTW_RODFT00;
	else
		return FFTW_REDFT00;
}

template <>
template <int D>
fftw_r2r_kind FFT<cat::array<RS,D>,cat::array<RS,D> >::r2r_inverse_kind(const string & subtype__)
{
	//cout << subtype << endl;
	if (subtype__=="sin")
		return FFTW_RODFT00;
	else
		return FFTW_REDFT00;
}
