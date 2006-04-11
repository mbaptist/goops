 
/////////////////////////////////////////
// SINE TRANSFORM FOR 1D SCALAR FIELDS //
/////////////////////////////////////////

#include "s_sinfft_1d.h"

#include <cat.h>

using namespace cat;

//IMPLEMENTATION
//constructor 
s_sinfft_1d::s_sinfft_1d(int size_):size(size_)
{
  //allocate temporary space for inline transform (as a real variable)
  //real_data=static_cast<Real *>(fftw_malloc(sizeof(Real) * (size-2)));
  real_data=new double[size-2];
  //initialize temporary space with zeros
  for(int i=0;i<size-2;++i)
    real_data[i]=0;
  //create plan for direct transforms
  direct_plan=fftw_plan_r2r_1d
    (size-2,real_data,real_data,FFTW_RODFT00,FFTW_ESTIMATE);
  //create plan for inverse transforms
  inverse_plan=fftw_plan_r2r_1d
    (size-2,real_data,real_data,FFTW_RODFT00,FFTW_ESTIMATE);
}
//Destructor
s_sinfft_1d::~s_sinfft_1d()
{  
  fftw_destroy_plan(direct_plan);
  fftw_destroy_plan(inverse_plan);
  //fftw_free(real_data);
  delete[] real_data;
}

//direct transform (scalar)
void s_sinfft_1d::direct_transform(RT& u_hat,const RT& u)
{
  //read real data from input
  for(int i=1;i<size-1;++i)
    real_data[i-1]=u(i);

  //perform the transform
  fftw_execute(direct_plan);

  //Normalise
  for (int i=0;i<size-2;++i)
    real_data[i] /= ((size-1));

  //write real data in output
  u_hat(0)=0;
  for(int i=1;i<size-1;++i)
    u_hat(i)=real_data[i-1];
  u_hat(size-1)=0;
}

//inverse transform (scalar)
void s_sinfft_1d::inverse_transform(RT& u,const RT& u_hat)
{
  //read real data from input
  for(int i=1;i<size-1;++i)
    real_data[i-1]=u_hat(i);

  //perform the transform
  fftw_execute(inverse_plan);

  for (int i=0;i<size-2;++i)
    real_data[i] /= 2.;

  //write real data in output
  u(0)=0;
  for(int i=1;i<size-1;++i)
    u(i)=real_data[i-1];
  u(size-1)=0;
}
