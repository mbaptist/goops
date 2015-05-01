// -*- C++ -*-
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



/////////////////////////////////////////////////
// Real TO COMPLEX TRANSFORM FOR VECTOR FIELDS //
/////////////////////////////////////////////////

#ifndef V_RFFT_H
#define V_RFFT_H

#include <complex>
using std::complex;

#include <fftw3.h>

#include "fftclass.h"

#include <cat.h>
using namespace cat;

#include "../goops_types.h"

//CLASS V_RFFT
//DECLARATION
template <int D,int N>
    class v_rfft
    {
    //Types for scalar fields
    typedef cat::Array<Real,D> RT;
    typedef cat::Array<complex<Real>,D> CT; 
    //Types for vector fields
    typedef cat::Array<cat::Tvector<Real,N>,D> VRT;
    typedef cat::Array<cat::Tvector<complex<Real>,N>,D> VCT;
    public:
    //Constructor
    v_rfft(const cat::Tvector<int,D> & size_):s_rfft_obj(){};
    //Destructor
    ~v_rfft(){};
    //Transformation functions for vectors
    void direct_transform(VCT& u_hat,const VRT& u);
    void inverse_transform(VRT& u,const VCT& u_hat);  
    
    private:
    FFT<cat::Array<double,D>,cat::Array<complex<double>,D> > s_rfft_obj;
    
};

//IMPLEMENTATION
//direct transform (vector)
template <int D,int N>
    void v_rfft<D,N>::direct_transform(VCT& u_hat,const VRT& u)
    {
    for(int comp=0;comp<N;++comp)
    {
	RT s_u(u.shape());
	CT s_u_hat(u_hat.shape());
	
	for(int i=0;i<s_u.size();++i)
	    s_u.data()[i]=(u.data()[i])[comp];
	
	
	s_rfft_obj.direct_transform(s_u_hat,s_u);
	
	for(int i=0;i<s_u_hat.size();++i)
	    (u_hat.data()[i])[comp]=s_u_hat.data()[i];  
	
    }
}

//inverse transform (vector)
template <int D,int N>
    void v_rfft<D,N>::inverse_transform(VRT& u,const VCT& u_hat)
    {
    for(int comp=0;comp<N;++comp)
    {
	RT s_u(u.shape());
	CT s_u_hat(u_hat.shape());
	
	for(int i=0;i<s_u_hat.size();++i)
	    s_u_hat.data()[i]=(u_hat.data()[i])[comp];
	
	
	s_rfft_obj.inverse_transform(s_u,s_u_hat);
	
	for(int i=0;i<s_u.size();++i)
	    (u.data()[i])[comp]=s_u.data()[i];  
	
    }
}






#endif
