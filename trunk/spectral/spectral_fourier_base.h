// -*- C++ -*-

#ifndef SPECTRAL_FOURIER_BASE_H
#define SPECTRAL_FOURIER_BASE_H

#include "../goops_types.h"
#include "../fft/fft.h"

#include <cat.h>

//Base Class to Fourier Spectral Methods
class SpectralFourierBase 
{
  //Members
protected:
  //sizes
	int n1,n2,n3;//number of points of the grid in real space
	int n1_hat,n2_hat,n3_hat;//number of points of the grid in fourierspace
	Real l1,l2,l3;//physical dimensions
	Real ar1,ar2,ar3;//aspect ratios
	RVF wv;//wavevectors
	RSF wv2;//Square of the norm of wavevectors
	Real wnmax;
	int nwn;
	Real wnstep;
  Real dealiasing_limit;
  cat::array<bool,3> dealiasing_mask;

	//Ctor
  SpectralFourierBase(const int & n1__,const int & n2__,const int & n3__,
                      const int & n1_hat__,const int & n2_hat__,const int & n3_hat__,
                      const Real & l1__,const Real & l2__,const Real & l3__,
                      const Real & ar1__,const Real & ar2__,const Real & ar3__);
	//Dtor
  virtual ~SpectralFourierBase();
	//Forbidden Ctors
private:
	SpectralFourierBase();
	SpectralFourierBase(const SpectralFourierBase &);	
      //Public methods
public:
  
  //dialiasing (for second order non-linearities)
	void dealias(CSF& field) const;//scalar field
	void dealias(CVF& field) const;//vector field

  //Laplacian
	CSF lap_hat(const CSF & field);//vector field
	CVF lap_hat(const CVF & field);//vector field

  //Solve lap(f)=g in fourier space
	CSF poisson_hat(const CSF & field);//scalar field
	CVF poisson_hat(const CVF & field);//vector field
/*	 template <class T>
		T poisson_hat(const ArrayExpression<T> & field);*/
		
  //print non-vanishing harmonics in fourier space
	void pnvh_hat(const CSF & field);//scalar field
	void pnvh_hat(const CVF & field);//vector field

	//function that evaluates the energy in real space
	Real energy(const RSF & field);
	Real energy(const RVF & field);
};

#endif







