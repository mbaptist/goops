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



#ifndef SPECTRAL_FOURIER_BASE_H
#define SPECTRAL_FOURIER_BASE_H

#include "../goops_types.h"
#include "../fft/fft.h"

#include <cat.h>

//Base Class to Fourier Spectral Methods
class SpectralFourierBase 
{
    //Members
    protected:
    //sizes
    int n1,n2,n3;//number of points of the grid in real space
    int n1_hat,n2_hat,n3_hat;//number of points of the grid in fourierspace
    Real l1,l2,l3;//physical dimensions
    Real ar1,ar2,ar3;//aspect ratios
    RVF wv;//wavevectors
    RSF wv2;//Square of the norm of wavevectors
    Real wnmax;
    int nwn;
    Real wnstep;
    Real dealiasing_limit;
    cat::Array<bool,3> dealiasing_mask;
    
    //Ctor
    SpectralFourierBase(const int & n1__,const int & n2__,const int & n3__,
			const int & n1_hat__,const int & n2_hat__,const int & n3_hat__,
			const Real & l1__,const Real & l2__,const Real & l3__,
			const Real & ar1__,const Real & ar2__,const Real & ar3__);
    //Dtor
    virtual ~SpectralFourierBase();
    //Forbidden Ctors
    private:
    SpectralFourierBase();
    SpectralFourierBase(const SpectralFourierBase &);	
    //Public methods
    public:
    
    //dialiasing (for second order non-linearities)
    void dealias(CSF& field) const;//scalar field
    void dealias(CVF& field) const;//vector field
    
    //Laplacian
    CSF lap_hat(const CSF & field);//vector field
    CVF lap_hat(const CVF & field);//vector field
    
    //Solve lap(f)=g in fourier space
    CSF poisson_hat(const CSF & field);//scalar field
    CVF poisson_hat(const CVF & field);//vector field
    
    //print non-vanishing harmonics in fourier space
    void pnvh_hat(const CSF & field);//scalar field
    void pnvh_hat(const CVF & field);//vector field
    
    //function that evaluates the energy in real space
    Real energy(const RSF & field);
    Real energy(const RVF & field);
};

#endif







