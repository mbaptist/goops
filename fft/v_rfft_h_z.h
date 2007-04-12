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


#ifndef V_RFFT_H_Z_CCS_H
#define V_RFFT_H_Z_CCS_H

//namespace goops
//{


#include <complex>
using std::complex;

#include "../goops_types.h"

#include <cat.h>

#include "fftclass.h"
#include "s_rfft_h_z.h"


//CLASS V_RFFT
//DECLARATION
class v_rfft_h_z
{
    //Types for scalar fields
    typedef cat::Array<Real,3> RT;
    typedef cat::Array<complex<Real>,3> CT; 
    //Types for vector fields
    typedef cat::Array<cat::Tvector<Real,3>,3> VRT;
    typedef cat::Array<cat::Tvector<complex<Real>,3>,3> VCT;
    public:
    //Constructor
    v_rfft_h_z(const string & subtype_x, const string & subtype_y,const string & subtype_z):
	s_rfft_x_obj(subtype_x),
	s_rfft_y_obj(subtype_y),
	s_rfft_z_obj(subtype_z)
	{
    };
    //Destructor
    ~v_rfft_h_z()
    {
    };
    //Transformation functions for vectors
    void direct_transform(VCT& u_hat,const VRT& u);
    void inverse_transform(VRT& u,const VCT& u_hat);  
    private:
    s_rfft_h_z s_rfft_x_obj;
    s_rfft_h_z s_rfft_y_obj;
    s_rfft_h_z s_rfft_z_obj;
};


//}

#endif
