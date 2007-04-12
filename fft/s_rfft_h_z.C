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
// namespace goops
// {
#include "s_rfft_h_z.h"
#include "../goops_types.h"
#include "fftclass.h"

#include <cat.h>
using namespace cat;

#include <string>

using namespace std;



//CLASS S_RFFT_H_Z
//IMPLEMENTATION

//Constructor
s_rfft_h_z::s_rfft_h_z(std::string subtype):
    s_rfft_obj(),
    fftz_obj(subtype)
    {
}
//Destructor
s_rfft_h_z::~s_rfft_h_z()
{
}

void s_rfft_h_z::direct_transform(CT3& u_hat,const RT3& u)
{
    cat::Tvector<int,3> size(u.shape());
    
    //Copies u to a working array
    RT3 work(u);
    //Transform along z
    for(int i=0;i<size[0];++i)
	for(int j=0;j<size[1];++j)
	{
	RT uz(size[2]);
	RT uz_hat(size[2]);
	for (int k=0;k<size[2];++k)
	    uz(k)=work(i,j,k);
	fftz_obj.direct_transform(uz_hat,uz);
	for (int k=0;k<size[2];++k)
	    work(i,j,k)=uz_hat(k);
    }
    
    // cout << "h2 " << sum(work*work) << endl;
    
    //Transform along x and y
    for (int k=0;k<size[2];++k)
    {
	RT2 uxy(size[0],size[1]);
	CT2 uxy_hat(size[0],size[1]/2+1);
	for(int i=0;i<size[0];++i)
	    for(int j=0;j<size[1];++j)
		uxy(i,j)=work(i,j,k);
	s_rfft_obj.direct_transform(uxy_hat,uxy);
	for(int i=0;i<size[0];++i)
	    for(int j=0;j<size[1]/2+1;++j)
		u_hat(i,j,k)=uxy_hat(i,j);
    }
}

void s_rfft_h_z::inverse_transform(RT3& u,const CT3& u_hat)
{  
    cat::Tvector<int,3> size(u.shape());
    
    //Creates a working array from u
    RT3 work(u.shape());
    
    
    
    //Transform along x and y
    for (int k=0;k<size[2];++k)
    {
	RT2 uxy(size[0],size[1]);
	CT2 uxy_hat(size[0],size[1]/2+1);
	for(int i=0;i<size[0];++i)
	    for(int j=0;j<size[1]/2+1;++j)
		uxy_hat(i,j)=u_hat(i,j,k);
	s_rfft_obj.inverse_transform(uxy,uxy_hat);
	for(int i=0;i<size[0];++i)
	    for(int j=0;j<size[1];++j)
		work(i,j,k)=uxy(i,j);
    }
    
    //Transform along z
    for(int i=0;i<size[0];++i)
	for(int j=0;j<size[1];++j)
	{
	RT uz(size[2]);
	RT uz_hat(size[2]);
	for (int k=0;k<size[2];++k)
	    uz_hat(k)=work(i,j,k);
	fftz_obj.inverse_transform(uz,uz_hat);
	for (int k=0;k<size[2];++k)
	    work(i,j,k)=uz(k);
    }
    //Copy work to u
    u=work;
}

//}

