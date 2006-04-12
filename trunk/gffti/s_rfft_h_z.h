// -*- C++ -*-

/////////////////////////////////////////////////////////////////////////
// Real TO COMPLEX TRANSFORM FOR 3D SCALAR FIELDS WITH SIN/COS ALONG Z //
/////////////////////////////////////////////////////////////////////////

#ifndef S_RFFT_H_Z_H
#define S_RFFT_H_Z_H



// namespace goops
// {

#include "fft.h"
#include "../goops_types.h"
#include <cat.h>

#include <string>

using namespace cat;

//CLASS S_RFFT_H_Z
//DECLARATION
class s_rfft_h_z
{
private:
//Types for scalar fields
  typedef cat::array<Real,1> RT;
  typedef cat::array<complex<Real>,1> CT; 
  typedef cat::array<Real,2> RT2;
  typedef cat::array<complex<Real>,2> CT2; 
  typedef cat::array<Real,3> RT3;
  typedef cat::array<complex<Real>,3> CT3; 
public:
  //Constructor
	s_rfft_h_z(std::string subtype);
  //Destructor
  ~s_rfft_h_z();
  void direct_transform(CT3& u_hat,const RT3& u);
  void inverse_transform(RT3& u,const CT3& u_hat);
private:
	//cat::tvector <int,3> size;
	FFT<cat::array<double,2>, cat::array<complex<double>,2> > s_rfft_obj;
	FFT<cat::array<double,1>, cat::array<double,1> > fftz_obj;
};

//}
#endif
