// -*- C++ -*-

/////////////////////////////////////////////////
// Real TO COMPLEX TRANSFORM FOR VECTOR FIELDS //
/////////////////////////////////////////////////

#ifndef V_RFFT_H
#define V_RFFT_H

#include <complex>
using std::complex;

#include <fftw3.h>

#include "s_rfft.h"

#include <cat.h>
using namespace cat;

#include "../goops_types.h"

//CLASS V_RFFT
//DECLARATION
template <int D,int N>
class v_rfft
{
  //Types for scalar fields
  typedef cat::array<Real,D> RT;
  typedef cat::array<complex<Real>,D> CT; 
  //Types for vector fields
  typedef cat::array<cat::tvector<Real,N>,D> VRT;
  typedef cat::array<cat::tvector<complex<Real>,N>,D> VCT;
 public:
  //Constructor
  v_rfft(const cat::tvector<int,D> & size_):s_rfft_obj(size_){};
  //Destructor
  ~v_rfft(){};
  //Transformation functions for vectors
  void direct_transform(VCT& u_hat,const VRT& u);
  void inverse_transform(VRT& u,const VCT& u_hat);  

 private:
  s_rfft<D> s_rfft_obj;

};

//IMPLEMENTATION
//direct transform (vector)
template <int D,int N>
void v_rfft<D,N>::direct_transform(VCT& u_hat,const VRT& u)
{
  for(int comp=0;comp<N;++comp)
    {
      RT s_u(u.shape());
      CT s_u_hat(u_hat.shape());

      for(int i=0;i<s_u.size();++i)
	s_u.data()[i]=(u.data()[i])[comp];
      

      s_rfft_obj.direct_transform(s_u_hat,s_u);

      for(int i=0;i<s_u_hat.size();++i)
	(u_hat.data()[i])[comp]=s_u_hat.data()[i];  

    }
}

//inverse transform (vector)
template <int D,int N>
void v_rfft<D,N>::inverse_transform(VRT& u,const VCT& u_hat)
{
  for(int comp=0;comp<N;++comp)
    {
      RT s_u(u.shape());
      CT s_u_hat(u_hat.shape());

      for(int i=0;i<s_u_hat.size();++i)
	s_u_hat.data()[i]=(u_hat.data()[i])[comp];
      

      s_rfft_obj.inverse_transform(s_u,s_u_hat);

      for(int i=0;i<s_u.size();++i)
	(u.data()[i])[comp]=s_u.data()[i];  

    }
}






#endif
