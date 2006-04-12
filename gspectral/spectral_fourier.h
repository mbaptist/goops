// -*- C++ -*-

#ifndef SPECTRAL_FOURIER_H
#define SPECTRAL_FOURIER_H

#include "../goops_types.h"
#include "../gffti/gffti.h"

#include <cat.h>

class spectral_fourier
{
  //Members
private:
  //sizes
  int n1,n2,n3;
  Real l1,l2,l3;
  
public:
  //wavevectors
  RVF wv;
    
  //Square of norm of wavevectors 
  RSF wv2;
  
private:
  Real dealiasing_limit;
  cat::array<bool,3> dealiasing_mask;

public:
  //Object to perform ffts
  v_rfft<3,3> fft;
  FFT<RSF,CSF> sfft;
  
        
  //Constructor and destructor
public:
  spectral_fourier(const int & n1__,
		   const int & n2__,
		   const int & n3__,
		   const Real & l1__,
		   const Real & l2__,
		   const Real & l3__);
  ~spectral_fourier();
      
      //Public methods
  public:
  
  //perform dialiasing on fields
  void dealias(CSF& field) const;//scalar
  void dealias(CVF& field) const;//vector

  
  //Solve lap(f)=g in fourier space
  CSF poisson_hat(const CSF & field);
  CVF poisson_hat(const CVF & field);



 //Derivatives in fourier space
 

  //Differential operators in fourier space
  //deivatives in horizontal directions - scalars
  CSF d_dx_index_hat(const CSF & field,const int index);
  //deivatives in horizontal directions - vectors
  CVF d_dx_index_hat(const CVF & field, const int index);

  //Gradient of scalar field
  CVF grad_hat(const CSF & field);
  //Divergence of vector field
  CSF div_hat(const CVF & field);
  //Curl of vector field
  CVF curl_hat(const CVF & field);
  //Laplacian of scalar field
  CSF lap_hat(const CSF & field);
  //Laplacian of vector field
  CVF lap_hat(const CVF & field);


  //remove gradient part after of a vfield
  CVF remove_gradient(CVF & field);


  //L2 scalar product
  Real scalar_prod(const CSF & x,
		     const CSF & y) const;
  Real scalar_prod(const CVF & x,
		     const CVF & y) const;

    //Energy spectrum
  cat::array<Real,1> eval_energ_spec(const CSF & field);
  cat::array<Real,1> eval_energ_spec(const CVF & field);
  cat::array<Real,1> eval_energ_spec(const CSF & field,const int & npoints);
  cat::array<Real,1> eval_energ_spec(const CVF & field,const int & npoints);


  //print non-vanishing harmonics
  void pnvh(const CVF & field);
  void pnvh(const CSF & field);
  
  void pnvh(const RVF & field);
  void pnvh(const RSF & field);
  
};

#endif







