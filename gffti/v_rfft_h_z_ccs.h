// -*- C++ -*-

/////////////////////////////////////////////////
// Real TO COMPLEX TRANSFORM FOR VECTOR FIELDS //
/////////////////////////////////////////////////


#ifndef V_RFFT_H_Z_CCS_H
#define V_RFFT_H_Z_CCS_H

#include <complex>
using std::complex;

#include "../goops_types.h"

#include <cat.h>

#include "s_cosfft_1d.h"
#include "s_sinfft_1d.h"
#include "s_rfft_h_z.h"


//CLASS V_RFFT
//DECLARATION
class v_rfft_h_z_ccs
{
  //Types for scalar fields
  typedef cat::array<Real,3> RT;
  typedef cat::array<complex<Real>,3> CT; 
 //Types for vector fields
  typedef cat::array<cat::tvector<Real,3>,3> VRT;
  typedef cat::array<cat::tvector<complex<Real>,3>,3> VCT;
 public:
  //Constructor
  v_rfft_h_z_ccs(cat::tvector<int,3> size_):
    s_rfft_cos_obj(size_),
    s_rfft_sin_obj(size_){};
  //Destructor
  ~v_rfft_h_z_ccs(){};
  //Transformation functions for vectors
  void direct_transform(VCT& u_hat,const VRT& u);
  void inverse_transform(VRT& u,const VCT& u_hat);  
 private:
  s_rfft_h_z<s_cosfft_1d> s_rfft_cos_obj;
  s_rfft_h_z<s_sinfft_1d> s_rfft_sin_obj;
};



#endif
