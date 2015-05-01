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

#ifndef SPECTRAL_FOURIER_LAYER_H
#define SPECTRAL_FOURIER_LAYER_H

#include "../goops_types.h"
#include "../fft/fft.h"

#include "spectral_fourier_base.h"

#include <cat.h>

class SpectralFourierLayer : public SpectralFourierBase
{
    //Members
    
    private:
    using SpectralFourierBase::n1;
    using SpectralFourierBase::n2;
    using SpectralFourierBase::n3;
    using SpectralFourierBase::n1_hat;
    using SpectralFourierBase::n2_hat;
    using SpectralFourierBase::n3_hat;
    using SpectralFourierBase::l1;
    using SpectralFourierBase::l2;
    using SpectralFourierBase::l3;
    using SpectralFourierBase::ar1;
    using SpectralFourierBase::ar2;
    using SpectralFourierBase::ar3;
    public:
    using SpectralFourierBase::wv;
    using SpectralFourierBase::wv2;
    using SpectralFourierBase::wnmax;
    using SpectralFourierBase::nwn;
    using SpectralFourierBase::wnstep;
    
    public:
    //Object to perform ffts
    v_rfft_h_z fft_ccs;
    v_rfft_h_z fft_ssc;
    s_rfft_h_z sfft_s;
    s_rfft_h_z sfft_c;
    
    //Constructor and destructor
    public:
    SpectralFourierLayer(const int & n1__,
			 const int & n2__,
			 const int & n3__,
			 const Real & l1__,
			 const Real & l2__,
			 const Real & l3__);
    ~SpectralFourierLayer();
    
    private:
    //Forbidden Ctors
    SpectralFourierLayer();
    SpectralFourierLayer(const SpectralFourierLayer &);
    
    private:
    void initialise();
    
    //Public methods
    public:
    
    //deivatives in horizontal directions - scalars
    CSF d_dx_index_hat(const CSF & field,const int index);
    //deivatives in horizontal directions - vectors
    CVF d_dx_index_hat(const CVF & field, const int index);
    
    //Gradient of scalar field
    CVF grad_hat(const CSF & field,const bool kind);
    //Divergence of vector field
    CSF div_hat(const CVF & field,const bool kind);
    //Curl of vector field
    CVF curl_hat(const CVF & field,const bool kind);
    
    //Extract gradient part
    //3D vector field
    CVF remove_gradient(CVF & field,const bool kind);
    
    //L2 scalar product
    Real scalar_prod(const CSF & x,const CSF & y,const bool kind) const;//scalar fields
    Real scalar_prod(const CSF & x,const CSF & y) const;// To be moved to linops
    Real scalar_prod(const CVF & x,const CVF & y,const bool kind) const;//vector fields
    Real scalar_prod(const CVF & x,const CVF & y) const;// To be moved to linops
    
    //Energy spectrum
    cat::Array<Real,1> eval_energ_spec(const CSF & field,const bool & kind);//scalar fields
    cat::Array<Real,1> eval_energ_spec(const CVF & field,const bool & kind);//vector fields
    
    //print non-vanishing harmonics in real space
    using SpectralFourierBase::pnvh_hat;
    void pnvh(const RSF & field,const bool & kind);//scalar fields
    void pnvh(const RVF & field,const bool & kind);//vector fields
    
};

#endif







