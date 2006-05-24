// -*- C++ -*-

////////////////////////////////////////////////////////////////////////////
// COMPLEX TO COMPLEX TRANSFORM FOR 3D SCALAR FIELDS WITH SIN/COS ALONG Z //
////////////////////////////////////////////////////////////////////////////

#ifndef S_FFT_H_Z_H
#define S_FFT_H_Z_H

#include "s_fft.h"

#include <cat.h>
using namespace cat;

//CLASS S_FFT_H_Z
//DECLARATION
template <class FFTZ>
class s_fft_h_z
{
  //Types for scalar fields
  typedef cat::array<Real,1> RT;
  typedef cat::array<complex<Real>,1> CT; 
  typedef cat::array<Real,2> RT2;
  typedef cat::array<complex<Real>,2> CT2; 
  typedef cat::array<Real,3> RT3;
  typedef cat::array<complex<Real>,3> CT3; 
public:
  //Constructor
  s_fft_h_z(cat::tvector<int,3> size_):
    size(size_),
    s_fft_obj(cat::tvector<int,2>(size[0],size[1])),
    fftz_obj(size[2]){};
  //Destructor
  ~s_fft_h_z(){};
  void direct_transform(CT3& u_hat,const CT3& u);
  void inverse_transform(CT3& u,const CT3& u_hat);
private:
  cat::tvector <int,3> size;
  s_fft<2> s_fft_obj;
  FFTZ fftz_obj;
};

template <class FFTZ>
void s_fft_h_z<FFTZ>::direct_transform(CT3& u_hat,const CT3& u)
{
  //Copies u to a working array
  CT3 work(u);
  //Transform along z
  for(int i=0;i<size[0];++i)
    for(int j=0;j<size[1];++j)
      {
	cat::tvector<int,1> ss;
	RT uz_re(size[2]);
	RT uz_hat_re(size[2]);
	RT uz_im(size[2]);
	RT uz_hat_im(size[2]);
	for (int k=0;k<size[2];++k)
	  {
	    uz_re(k)=work(i,j,k).real();
	    uz_im(k)=work(i,j,k).imag();
	  }
	fftz_obj.direct_transform(uz_hat_re,uz_re);
	fftz_obj.direct_transform(uz_hat_im,uz_im);
	for (int k=0;k<size[2];++k)
	  work(i,j,k)=complex<Real>(uz_hat_re(k),uz_hat_im(k));
      }

  //Transform along x and y
  for (int k=0;k<size[2];++k)
    {
      CT2 uxy(size[0],size[1]);
      CT2 uxy_hat(size[0],size[1]);
      for(int i=0;i<size[0];++i)
	for(int j=0;j<size[1];++j)
	  uxy(i,j)=work(i,j,k);
      s_fft_obj.direct_transform(uxy_hat,uxy);
      for(int i=0;i<size[0];++i)
	for(int j=0;j<size[1];++j)
	  u_hat(i,j,k)=uxy_hat(i,j);
    }

}
	

template <class FFTZ>
void s_fft_h_z<FFTZ>::inverse_transform(CT3& u,const CT3& u_hat)
{  

  //Creates a working array from u
  CT3 work(u);
  
  //Transform along x and y
  for (int k=0;k<size[2];++k)
    {
      CT2 uxy(size[0],size[1]);
      CT2 uxy_hat(size[0],size[1]);
      for(int i=0;i<size[0];++i)
	for(int j=0;j<size[1];++j)
	  uxy_hat(i,j)=u_hat(i,j,k);
      s_fft_obj.inverse_transform(uxy,uxy_hat);
      for(int i=0;i<size[0];++i)
	for(int j=0;j<size[1];++j)
	  work(i,j,k)=uxy(i,j);	  
    }
  //Transform along z
  for(int i=0;i<size[0];++i)
    for(int j=0;j<size[1];++j)
      {
	RT uz_re(size[2]);
	RT uz_hat_re(size[2]);
	RT uz_im(size[2]);
	RT uz_hat_im(size[2]);
	for (int k=0;k<size[2];++k)
	  {
	    uz_hat_re(k)=work(i,j,k).real();
	    uz_hat_im(k)=work(i,j,k).imag();
	  }
	fftz_obj.inverse_transform(uz_re,uz_hat_re);
	fftz_obj.inverse_transform(uz_im,uz_hat_im);
	for (int k=0;k<size[2];++k)
	  work(i,j,k)=complex<Real>(uz_re(k),uz_im(k));
      }
  //Copy work to u
  u=work;
}


#endif
