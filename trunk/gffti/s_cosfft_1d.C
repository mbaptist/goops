 
///////////////////////////////////////////
// COSINE TRANSFORM FOR 1D SCALAR FIELDS //
///////////////////////////////////////////

#include "s_cosfft_1d.h"

#include <fftw3.h>

#include <cat.h>
using namespace cat;

//IMPLEMENTATION
//constructor 

#define PLAN fftw_plan_r2r_1d

s_cosfft_1d::s_cosfft_1d(int size_):size(size_)
{
  //allocate temporary space for inline transform (as a real variable)
  //real_data=static_cast<Real *>(fftw_malloc(sizeof(Real) * (size)));
  real_data=new double[size];
  //initialize temporary space with zeros
  for(int i=0;i<size;++i)
    real_data[i]=0;
  //create plan for direct transforms
  direct_plan =fftw_plan_r2r_1d
    (size,real_data,real_data,FFTW_REDFT00,FFTW_ESTIMATE);
  //create plan for inverse transforms
  inverse_plan=fftw_plan_r2r_1d
    (size,real_data,real_data,FFTW_REDFT00,FFTW_ESTIMATE);
}
//Destructor
s_cosfft_1d::~s_cosfft_1d()
{  
  fftw_destroy_plan(direct_plan);
  fftw_destroy_plan(inverse_plan);
  //fftw_free(real_data);
  delete[] real_data;
}

//direct transform (scalar)
void s_cosfft_1d::direct_transform(RT& u_hat,const RT& u)
{
  //read real data from input
  for(int i=0;i<size;++i)
    real_data[i]=u(i);
 

  //perform the transform
  fftw_execute(direct_plan);

  // //Normalise
   for (int i=0;i<size;++i)
    real_data[i] /= (size-1);

  //Factor 2
  real_data[0] *= .5;



  //write real data in output
  for(int i=0;i<size;++i)
    u_hat(i)=real_data[i];
}

//inverse transform (scalar)
void s_cosfft_1d::inverse_transform(RT& u,const RT& u_hat)
{
  //read real data from input
  for(int i=0;i<size;++i)
    real_data[i]=u_hat(i);

  real_data[size-1]=0;

  //perform the transform
  fftw_execute(inverse_plan);
  
  for (int i=0;i<size;++i)
    {
      real_data[i] += u_hat(0); 
      real_data[i] *= .5;
    }

  //write real data in output
  for(int i=0;i<size;++i)
    u(i)=real_data[i];

}
