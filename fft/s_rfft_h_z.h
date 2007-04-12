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



/////////////////////////////////////////////////////////////////////////
// Real TO COMPLEX TRANSFORM FOR 3D SCALAR FIELDS WITH SIN/COS ALONG Z //
/////////////////////////////////////////////////////////////////////////

#ifndef S_RFFT_H_Z_H
#define S_RFFT_H_Z_H



// namespace goops
// {

#include "fftclass.h"
#include "../goops_types.h"
#include <cat.h>

#include <string>

using namespace cat;

//CLASS S_RFFT_H_Z
//DECLARATION
class s_rfft_h_z
{
    private:
    //Types for scalar fields
    typedef cat::Array<Real,1> RT;
    typedef cat::Array<complex<Real>,1> CT; 
    typedef cat::Array<Real,2> RT2;
    typedef cat::Array<complex<Real>,2> CT2; 
    typedef cat::Array<Real,3> RT3;
    typedef cat::Array<complex<Real>,3> CT3; 
    public:
    //Constructor
    s_rfft_h_z(std::string subtype);
    //Destructor
    ~s_rfft_h_z();
    void direct_transform(CT3& u_hat,const RT3& u);
    void inverse_transform(RT3& u,const CT3& u_hat);
    private:
    //cat::tvector <int,3> size;
    FFT<cat::Array<double,2>, cat::Array<complex<double>,2> > s_rfft_obj;
    FFT<cat::Array<double,1>, cat::Array<double,1> > fftz_obj;
    s_rfft_h_z(const s_rfft_h_z &);
    s_rfft_h_z();
    
};

//}
#endif
