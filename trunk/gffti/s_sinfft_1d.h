// -*- C++ -*-

/////////////////////////////////////////
// SINE TRANSFORM FOR 1D SCALAR FIELDS //
/////////////////////////////////////////

#ifndef S_SINFFT_1D_H
#define S_SINFFT_1D_H

#include <fftw3.h>

#include <cat.h>

#include "../goops_types.h"

//CLASS S_SINFFT_1D
//DECLARATION
class s_sinfft_1d
{
  //Type for real scalar fields
  typedef cat::array<Real,1> RT; 
 public:
  //Constructor
  s_sinfft_1d(int size_);
  //Destructor
  ~s_sinfft_1d();
  int size;//Size of the vector to be transformed, along each dimension
  //plan variables
  fftw_plan direct_plan;//direct transforms
  fftw_plan inverse_plan;//inverse transforms
  //pointers to local data
  double *real_data;//Pointer to real data
  //Transformation functions
  void direct_transform(RT& u_hat,const RT& u);
  void inverse_transform(RT& u,const RT& u_hat);
};


#endif
