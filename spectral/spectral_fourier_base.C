
#include "spectral_fourier_base.h"

#include "globals.h"

#include "../goops_types.h"
#include "../fft/fft.h"

#include <cat.h>

#include <iomanip>

#include <complex>
#include <fstream>

using namespace std;
using namespace cat;

//Construtor
SpectralFourierBase::SpectralFourierBase(const int & n1__,const int & n2__,const int & n3__,
                                         const int & n1_hat__,const int & n2_hat__,const int & n3_hat__,
                                         const Real & l1__,const Real & l2__,const Real & l3__,
                                         const Real & ar1__,const Real & ar2__,const Real & ar3__):
n1(n1__),
n2(n2__),
n3(n3__),
n1_hat(n1_hat__),
n2_hat(n2_hat__),
n3_hat(n3_hat__),
l1(l1__),
l2(l2__),
l3(l3__),
ar1(ar1__),
ar2(ar2__),
ar3(ar3__),
wv(n1_hat__,n2_hat__,n3_hat__),
wv2(n1_hat__,n2_hat__,n3_hat__),
dealiasing_limit(0),
dealiasing_mask(n1_hat__,n2_hat__,n3_hat__)
{
}

  //Desctructor
SpectralFourierBase::~SpectralFourierBase()
{
}

  //Dealiasing
  //performs dealiasing for second order non-linearities
  //scalar fields
void SpectralFourierBase::dealias(CSF& field) const
{
	field*=dealiasing_mask;
}
  //vector fields
void SpectralFourierBase::dealias(CVF& field) const
{
	field*=dealiasing_mask;
}

//Laplacian of scalar field
CSF SpectralFourierBase::lap_hat(const CSF & field)
{
	return CSF(-wv2*field);
}
//Laplacian of vector field
CVF SpectralFourierBase::lap_hat(const CVF & field)
{
	return CVF(-wv2*field);
}

  //Solve lap(f)=g in fourier space - scalars
CSF SpectralFourierBase::poisson_hat(const CSF & field)
{
	CSF out(field.shape());
	out=field;
	wv2(0,0,0)=1;
	out/=(-wv2);
	wv2(0,0,0)=1e-30;
	return out;
}
  //Solve lap(f)=g in fourier space - vectors
CVF SpectralFourierBase::poisson_hat(const CVF & field)
{
	CVF out(field.shape());
	out=field;
	wv2(0,0,0)=1;
	out/=(-wv2);
	wv2(0,0,0)=1e-30;
	return out;
}

//Printing non-vanishing harmonics in a field
//scalar fields in fourier space
void SpectralFourierBase::pnvh_hat(const CSF & field)
{
	CSF::const_iterator field_iter(field);
	for(field_iter=field.begin();field_iter!=field.end();++field_iter)
		if(abs(*field_iter)>1e-10)
			cout << field_iter.indices() << " " << *field_iter << endl;
}
  //vector fields in fourier space
void SpectralFourierBase::pnvh_hat(const CVF & field)
{
	CVF::const_iterator field_iter(field);
	for (field_iter=field.begin();field_iter!=field.end();++field_iter)
		if(norm(*field_iter)>1e-10)
			cout << field_iter.indices() << " " << *field_iter << endl;
}
