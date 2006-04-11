// -*- C++ -*-

////////////////////////////////////////////////////
// COMPLEX TO COMPLEX TRANSFORM FOR SCALAR FIELDS //
////////////////////////////////////////////////////

#ifndef S_FFT_H
#define S_FFT_H

#include <complex>
using std::complex;

#include <fftw3.h>

#include <cat.h>
using namespace cat;

#include "../goops_types.h"

//CLASS S_RFFT
//DECLARATION
template <int D>
class s_fft
{
  //Types for scalar fields
  typedef cat::array<complex<Real>,D> CT; 
 public:
  //Constructor
  s_fft(cat::tvector<int,D> size_);
  //Destructor
  ~s_fft();
  int *size;//Size of the vector to be transformed, along each dimension
  int total_size;
  //plan variables
  fftw_plan direct_plan;//direct transforms
  fftw_plan inverse_plan;//inverse transforms
  //pointers to local data
  fftw_complex *data;//Pointer to complex data
  //Transformation functions
  void direct_transform(CT& u_hat,const CT& u);
  void inverse_transform(CT& u,const CT& u_hat);
};

//IMPLEMENTATION
//constructor
template <int D> 
s_fft<D>::s_fft(cat::tvector<int,D> size_)
{
  //initialize size
  size=new int[D];
  for (int i=0;i<D;++i)
    size[i]=size_[i];
  //evaluate total size
  total_size=1;
  for (int i=0;i<D;++i)
    total_size *= size[i];
  //allocate temporary space for inline transform
  data=static_cast<fftw_complex *>
    (fftw_malloc(sizeof(fftw_complex) * total_size));
  //initialize temporary space with zeros
  for(int i=0;i<total_size;++i)
    {
      data[i][0]=0;
      data[i][1]=0;
    }
 //create plan for direct transforms
  direct_plan=fftw_plan_dft(D,size,data,data,FFTW_FORWARD,FFTW_ESTIMATE);
  //create plan for inverse transforms
  inverse_plan=fftw_plan_dft(D,size,data,data,FFTW_BACKWARD,FFTW_ESTIMATE); 
}
//Destructor
template <int D> 
s_fft<D>::~s_fft()
{  
  fftw_destroy_plan(direct_plan);
  fftw_destroy_plan(inverse_plan);
  data=0;
  fftw_free(data);
  delete[] size;
}

//direct transform (scalar)
template <int D>
void s_fft<D>::direct_transform(CT& u_hat,const CT& u)
{
  //read real data from input field via pointer
  for(int i=0;i<u.size();++i)
    {
      data[i][0]=(u.data()[i]).real();
      data[i][1]=(u.data()[i]).imag();
    }

  //perform the transform
  fftw_execute(direct_plan);

  //Normalize output
  for (int i=0;i<total_size;i++)
    {
      data[i][0] /= total_size;
      data[i][1] /= total_size;
    }

  //write complex data in output field via pointer
  for(int i=0;i<u_hat.size();++i)
    u_hat.data()[i]=complex<Real>(data[i][0],data[i][1]);
}

//inverse transform (scalar)
template <int D>
void s_fft<D>::inverse_transform(CT& u,const CT& u_hat)
{
  //read real data from input field via pointer
  for(int i=0;i<u_hat.size();++i)
    {
      data[i][0]=(u_hat.data()[i]).real();
      data[i][1]=(u_hat.data()[i]).imag();
    }

  //perform the transform
  fftw_execute(inverse_plan);

  //write complex data in output field via pointer
  for(int i=0;i<u.size();++i)
    u.data()[i]=complex<Real>(data[i][0],data[i][1]);
}

#endif
