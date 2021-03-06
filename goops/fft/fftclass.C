/*

Copyright 2004,2005,2006 Manuel Baptista

This file is part of GOOPS

GOOPS is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

GOOPS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/




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
    FFT<cat::Array<CS,D>,cat::Array<CS,D> >::FFT():
    direct_plan(FFTW_FORWARD),
    inverse_plan(FFTW_BACKWARD)
    {
}

template <>
    template<int D>
    FFT<cat::Array<CS,D>,cat::Array<CS,D> >::~FFT()
    {
}

//1D Real to Real Ctor
template <>
    template <int D>
    FFT<cat::Array<RS,D>,cat::Array<RS,D> >::FFT(const string & subtype__):
    subtype(subtype__),
    direct_plan(r2r_direct_kind(subtype__)),
    inverse_plan(r2r_inverse_kind(subtype__))
    {
    if (subtype=="sin")
    {
	direct_transform_imp=&FFT<cat::Array<RS,D>,cat::Array<RS,D> >::sin_direct_transform;
	inverse_transform_imp=&FFT<cat::Array<RS,D>,cat::Array<RS,D> >::sin_inverse_transform;
    }
    else
    {
	direct_transform_imp=&FFT<cat::Array<RS,D>,cat::Array<RS,D> >::cos_direct_transform;
	inverse_transform_imp=&FFT<cat::Array<RS,D>,cat::Array<RS,D> >::cos_inverse_transform;
    }
}

template <>
    template <int D>
    FFT<cat::Array<RS,D>,cat::Array<RS,D> >::~FFT()
    {
}

//1D Real to Real direct transform
template <>
    template <int D>
    void FFT<cat::Array<RS,D>,cat::Array<RS,D> >::direct_transform(cat::Array<RS,D> & fourierfield, const cat::Array<RS,D> & realfield)
    {
    (this->*direct_transform_imp)(fourierfield,realfield);
}

//1D Real to Real inverse transform
template <>
    template <int D>
    void FFT<cat::Array<RS,D>,cat::Array<RS,D> >::inverse_transform(cat::Array<RS,D> & realfield,const cat::Array<RS,D> & fourierfield)
    {
    (this->*inverse_transform_imp)(realfield,fourierfield);
}

//1D Real to Real direct sine transform 
template <>
    template <int D>
    void FFT<cat::Array<RS,D>,cat::Array<RS,D> >::sin_direct_transform(cat::Array<RS,D> & fourierfield, const cat::Array<RS,D> & realfield)
    {
    //		cat::Array<RS,D> ff=cat::Array<RS,D>(cat::tvector<int,D>(fourierfield.size()-2),fourierfield.data()+1);
    //  const cat::Array<RS,D> rf=cat::Array<RS,D>(cat::tvector<int,D>(realfield.size()-2),const_cast<RS *>(realfield.data())+1);
    
    size_t len=fourierfield.length()-2;
    cat::Array<RS,D> ff=cat::Array<RS,D>(cat::Tvector<int,D>(fourierfield.size()-2),fourierfield.ordering(),fourierfield.stride(),len,fourierfield.data()+1);
    
    // 		Array<T,D>::Array(tvector<int,D> & shape__,
    // 		                  tvector<int,D> & ordering__,
    // 		                  tvector<int,D> & stride__,
    // 		                  size_t & length__,
    // 		                  T * data__):
    
    len=realfield.length()-2;
    const cat::Tvector<int,D> sha(realfield.size()-2);
    const cat::Tvector<int,D> ord(realfield.ordering());
    cat::Tvector<int,D> str(realfield.stride());
    const RS * p = realfield.data()+1;
    const cat::Array<RS,D> rf(sha,ord,str,len,p);
    
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
    void FFT<cat::Array<RS,D>,cat::Array<RS,D> >::sin_inverse_transform(cat::Array<RS,D> & realfield,const cat::Array<RS,D> & fourierfield)
    {
    //cout << "SINE INVERSE TRANSFORM" << endl;
    //cat::Array<RS,D> rf=cat::Array<RS,D>(cat::tvector<int,D>(realfield.size()-2),realfield.data()+1);
    size_t len=realfield.length()-2;
    cat::Array<RS,D> rf(cat::Tvector<int,D>(realfield.size()-2),
			realfield.ordering(),
			realfield.stride(),
			len,
			realfield.data()+1);
    
    // 		Array<T,D>::Array(tvector<int,D> & shape__,
    // 		                  tvector<int,D> & ordering__,
    // 		                  tvector<int,D> & stride__,
    // 		                  size_t & length__,
    // 		                  T * data__):
    
    //const cat::Array<RS,D> ff=cat::Array<RS,D>(cat::tvector<int,D>(fourierfield.size()-2),const_cast<RS *>(fourierfield.data())+1);
    
    len=fourierfield.length()-2;
    const cat::Tvector<int,D> sha(fourierfield.size()-2);
    const cat::Tvector<int,D> ord(fourierfield.ordering());
    cat::Tvector<int,D> str(fourierfield.stride());
    const RS * p = fourierfield.data()+1;
    const cat::Array<RS,D> ff(sha,ord,str,len,p);
    
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
    void FFT<cat::Array<RS,D>,cat::Array<RS,D> >::cos_direct_transform(cat::Array<RS,D> & fourierfield, const cat::Array<RS,D> & realfield)
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
    void FFT<cat::Array<RS,D>,cat::Array<RS,D> >::cos_inverse_transform(cat::Array<RS,D> & realfield,const cat::Array<RS,D> & fourierfield)
    {
    // cout << "COS inverse transform" << endl;
    // 		inverse_plan.switch_data(realfield,fourierfield);
    // 		inverse_plan.execute();
    
    //cat::Array<RS,D> ff(fourierfield);
    //ff(ff.size()-1)=0;
    
    inverse_plan.execute(realfield,fourierfield);
    
    RS f0=fourierfield(0);
    RS fnm1=fourierfield(fourierfield.size()-1);
    // 	typename cat::Array<RS,D>::iterator iter(realfield);
    // 	int sign=1;
    // 	for(iter=realfield.begin();iter!=realfield.end();++iter)
    // 	{
    // 		(*iter)+=f0+sign*fnm1;
    // 		(*iter)*=.5;
    // 		sign*=-1;
    // 	}
    realfield+=f0;
    realfield*=.5;
    
}

template <>
    template <int D>
    fftw_r2r_kind FFT<cat::Array<RS,D>,cat::Array<RS,D> >::r2r_direct_kind(const string & subtype__)
    {
    if (subtype__=="sin")
	return FFTW_RODFT00;
    else
	return FFTW_REDFT00;
}

template <>
    template <int D>
    fftw_r2r_kind FFT<cat::Array<RS,D>,cat::Array<RS,D> >::r2r_inverse_kind(const string & subtype__)
    {
    //cout << subtype << endl;
    if (subtype__=="sin")
	return FFTW_RODFT00;
    else
	return FFTW_REDFT00;
}
