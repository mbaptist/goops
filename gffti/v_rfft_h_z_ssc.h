// -*- C++ -*-

/////////////////////////////////////////////////
// Real TO COMPLEX TRANSFORM FOR VECTOR FIELDS //
/////////////////////////////////////////////////


#ifndef V_RFFT_H_Z_SSC_H
#define V_RFFT_H_Z_SSC_H

// namespace goops
// {

#include <complex>
using std::complex;

#include "../goops_types.h"

#include <cat.h>

#include "fft.h"
#include "s_rfft_h_z.h"


//CLASS V_RFFT
//DECLARATION
class v_rfft_h_z_ssc
{
 //Types for scalar fields
  typedef cat::array<Real,3> RT;
  typedef cat::array<complex<Real>,3> CT; 
 //Types for vector fields
  typedef cat::array<cat::tvector<Real,3>,3> VRT;
  typedef cat::array<cat::tvector<complex<Real>,3>,3> VCT;

 public:
  //Constructor
  v_rfft_h_z_ssc(cat::tvector<int,3> size_):
    s_rfft_x_obj("sin"),
		s_rfft_y_obj("sin"),
    s_rfft_z_obj("cos"){};
  //Destructor
  ~v_rfft_h_z_ssc(){};
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
