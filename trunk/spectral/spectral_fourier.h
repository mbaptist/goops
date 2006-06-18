// -*- C++ -*-

#ifndef SPECTRAL_FOURIER_H
#define SPECTRAL_FOURIER_H

#include "spectral_fourier_base.h"

#include "../goops_types.h"
#include "../fft/fft.h"

#include <cat.h>

class SpectralFourier : public SpectralFourierBase
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
  v_rfft<3,3> fft;
  FFT<RSF,CSF> sfft;
  
        
  //Constructor and destructor
public:
  SpectralFourier(const int & n1__,
		   const int & n2__,
		   const int & n3__,
		   const Real & l1__,
		   const Real & l2__,
		   const Real & l3__);
  ~SpectralFourier();
	
private:
	//Forbidden Ctors
	SpectralFourier();
	SpectralFourier(const SpectralFourier &);
	
private:
	void initialise();
	
      //Public methods
  public:
  
  //deivatives in respect to coordinate x_index - scalars
	CSF d_dx_index_hat(const CSF & field,const int index);
  //deivatives in respect to coordinate x_index - vectors
	CVF d_dx_index_hat(const CVF & field, const int index);

  //Gradient of scalar field
  CVF grad_hat(const CSF & field);
  //Divergence of vector field
  CSF div_hat(const CVF & field);
  //Curl of vector field
  CVF curl_hat(const CVF & field);

  //remove gradient part after of a vfield
  CVF remove_gradient(CVF & field);
	
  //L2 scalar product
	Real scalar_prod(const CSF & x,const CSF & y) const;//scalar fields
	Real scalar_prod(const CVF & x,const CVF & y) const;//vector fields

    //Energy spectrum
	cat::array<Real,1> eval_energ_spec(const CSF & field);//scalar fields
	cat::array<Real,1> eval_energ_spec(const CVF & field);//vector fields

  //print non-vanishing harmonics
	using SpectralFourierBase::pnvh_hat;
	void pnvh(const RVF & field);//scalar fields
	void pnvh(const RSF & field);//vector fields
  
};

#endif
