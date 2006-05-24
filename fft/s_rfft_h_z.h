// -*- C++ -*-

/////////////////////////////////////////////////////////////////////////
// Real TO COMPLEX TRANSFORM FOR 3D SCALAR FIELDS WITH SIN/COS ALONG Z //
/////////////////////////////////////////////////////////////////////////

#ifndef S_RFFT_H_Z_H
#define S_RFFT_H_Z_H

#include "s_rfft.h"

#include <cat.h>
using namespace cat;

//CLASS S_RFFT_H_Z
//DECLARATION
template <class FFTZ>
class s_rfft_h_z
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
  s_rfft_h_z(cat::tvector<int,3> size_):
    size(size_),
    s_rfft_obj(cat::tvector<int,2>(size[0],size[1])),
    fftz_obj(size[2]){};
  //Destructor
  ~s_rfft_h_z(){};
  void direct_transform(CT3& u_hat,const RT3& u);
  void inverse_transform(RT3& u,const CT3& u_hat);
private:
  cat::tvector <int,3> size;
  s_rfft<2> s_rfft_obj;
  FFTZ fftz_obj;
};

template <class FFTZ>
void s_rfft_h_z<FFTZ>::direct_transform(CT3& u_hat,const RT3& u)
{
  //Copies u to a working array
  RT3 work(u);
  //Transform along z
  for(int i=0;i<size[0];++i)
    for(int j=0;j<size[1];++j)
      {
	RT uz(size[2]);
	RT uz_hat(size[2]);
	for (int k=0;k<size[2];++k)
	  uz(k)=work(i,j,k);
	fftz_obj.direct_transform(uz_hat,uz);
	for (int k=0;k<size[2];++k)
	  work(i,j,k)=uz_hat(k);
      }

  // cout << "h2 " << sum(work*work) << endl;

  //Transform along x and y
  for (int k=0;k<size[2];++k)
    {
      RT2 uxy(size[0],size[1]);
      CT2 uxy_hat(size[0],size[1]/2+1);
      for(int i=0;i<size[0];++i)
	for(int j=0;j<size[1];++j)
	  uxy(i,j)=work(i,j,k);
      s_rfft_obj.direct_transform(uxy_hat,uxy);
      for(int i=0;i<size[0];++i)
	for(int j=0;j<size[1]/2+1;++j)
	  u_hat(i,j,k)=uxy_hat(i,j);
    }

#if 0 
  cout << "h1" << endl;
    for(int k1=0;k1<size[0];++k1)
    for(int k2=0;k2<size[1]/2+1;++k2)    
      for(int k3=0;k3<size[2];++k3)
	if ( abs( u_hat(k1,k2,k3)) > 1e-12 )
	  cout << u_hat(k1,k2,k3) << endl;
    cout << endl;
#endif

}
	

template <class FFTZ>
void s_rfft_h_z<FFTZ>::inverse_transform(RT3& u,const CT3& u_hat)
{  

  //Creates a working array from u
  RT3 work(u.shape());

  //Transform along x and y
  for (int k=0;k<size[2];++k)
    {
      RT2 uxy(size[0],size[1]);
      CT2 uxy_hat(size[0],size[1]/2+1);
      for(int i=0;i<size[0];++i)
	for(int j=0;j<size[1]/2+1;++j)
	  uxy_hat(i,j)=u_hat(i,j,k);
      s_rfft_obj.inverse_transform(uxy,uxy_hat);
      for(int i=0;i<size[0];++i)
	for(int j=0;j<size[1];++j)
	  work(i,j,k)=uxy(i,j);	  
    }

  //Transform along z
  for(int i=0;i<size[0];++i)
    for(int j=0;j<size[1];++j)
      {
	RT uz(size[2]);
	RT uz_hat(size[2]);
	for (int k=0;k<size[2];++k)
	  uz_hat(k)=work(i,j,k);
	fftz_obj.inverse_transform(uz,uz_hat);
	for (int k=0;k<size[2];++k)
	  work(i,j,k)=uz(k);
      }
  //Copy work to u
  u=work;
}


#endif
