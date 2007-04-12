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


#include "spectral_fourier_base.h"

#include "globals.h"

#include "../goops_types.h"
#include "../fft/fft.h"

#include <cat.h>

#include <iomanip>

#include <complex>
#include <fstream>

using namespace std;
using namespace cat;

//Construtor
SpectralFourierBase::SpectralFourierBase(const int & n1__,const int & n2__,const int & n3__,
					 const int & n1_hat__,const int & n2_hat__,const int & n3_hat__,
					 const Real & l1__,const Real & l2__,const Real & l3__,
					 const Real & ar1__,const Real & ar2__,const Real & ar3__):
    n1(n1__),
    n2(n2__),
    n3(n3__),
    n1_hat(n1_hat__),
    n2_hat(n2_hat__),
    n3_hat(n3_hat__),
    l1(l1__),
    l2(l2__),
    l3(l3__),
    ar1(ar1__),
    ar2(ar2__),
    ar3(ar3__),
    wv(n1_hat__,n2_hat__,n3_hat__),
    wv2(n1_hat__,n2_hat__,n3_hat__),
    wnmax(0.),
    wnstep(1.),
    dealiasing_limit(0.),
    dealiasing_mask(n1_hat__,n2_hat__,n3_hat__)
    {
}

//Desctructor
SpectralFourierBase::~SpectralFourierBase()
{
}

//Dealiasing
//performs dealiasing for second order non-linearities
//scalar fields
void SpectralFourierBase::dealias(CSF& field) const
{
    field*=dealiasing_mask;
}
//vector fields
void SpectralFourierBase::dealias(CVF& field) const
{
    field*=dealiasing_mask;
}

//Laplacian of scalar field
CSF SpectralFourierBase::lap_hat(const CSF & field)
{
    return CSF(-wv2*field);
}
//Laplacian of vector field
CVF SpectralFourierBase::lap_hat(const CVF & field)
{
    return CVF(-wv2*field);
}

//Solve lap(f)=g in fourier space - scalars
CSF SpectralFourierBase::poisson_hat(const CSF & field)
{
    CSF out(field);
    wv2(0,0,0)=1.;
    out/=(-wv2);
    wv2(0,0,0)=1e-30;
    out(0,0,0)=0.;
    return out;
}
//Solve lap(f)=g in fourier space - vectors
CVF SpectralFourierBase::poisson_hat(const CVF & field)
{
    CVF out(field);
    wv2(0,0,0)=1.;
    out/=(-wv2);
    wv2(0,0,0)=1e-30;
    out(0,0,0)=0.;
    return out;
}

//Printing non-vanishing harmonics in a field
//scalar fields in fourier space
void SpectralFourierBase::pnvh_hat(const CSF & field)
{
    CSF::const_iterator field_iter(field);
    for(field_iter=field.begin();field_iter!=field.end();++field_iter)
	if(abs(*field_iter)>1e-10)
	    cout << field_iter.indices() << " " << *field_iter << endl;
}
//vector fields in fourier space
void SpectralFourierBase::pnvh_hat(const CVF & field)
{
    CVF::const_iterator field_iter(field);
    for (field_iter=field.begin();field_iter!=field.end();++field_iter)
	if(norm(*field_iter)>1e-10)
	    cout << field_iter.indices() << " " << *field_iter << endl;
}

Real SpectralFourierBase::energy(const RSF & field)
{
    //Performs the integral in <field,field> using the trapezoidal rule
    Real tmp=0;
    for(int i=0;i<n1;++i)
	for(int j=0;j<n2;++j)
	    for(int k=0;k<n3;++k)
	    {
	    if(k==0||k==n3-1)
		tmp+=field(i,j,k)*field(i,j,k);
	    else
		tmp+=2.*field(i,j,k)*field(i,j,k);
	}
    return  .5*tmp*(l1*l2*l3)/(n1*n2*(n3-1))/2.;
}


Real SpectralFourierBase::energy(const RVF & field)
{
    //Performs the integral in <field,field> using the trapezoidal rule
    Real tmp=0;
    for(int i=0;i<n1;++i)
	for(int j=0;j<n2;++j)
	    for(int k=0;k<n3;++k)
	    {
	    if(k==0||k==n3-1)
		tmp+=norm_sq(field(i,j,k));
	    else
		tmp+=2.*norm_sq(field(i,j,k));
	}
    return  .5*tmp*(l1*l2*l3)/(n1*n2*(n3-1))/2.;
}

