
#include "spectral_fourier.h"

#include "globals.h"

#include "../goops_types.h"
#include "../fft/fft.h"

#include <cat.h>

using namespace std;
using namespace cat;

//Construtor
SpectralFourier::SpectralFourier(const int & n1__,const int & n2__,const int & n3__,
                                 const Real & l1__,const Real & l2__,const Real & l3__):
SpectralFourierBase(n1__,n2__,n3__,
                    n1__,n2__,n3__/2+1,
                    l1__,l2__,l3__,
                    2*M_PI/l1__,2*M_PI/l2__,M_PI/l3__),
fft(cat::tvector<int,3>(n1,n2,n3)),
sfft()
{
	initialise();
}

void SpectralFourier::initialise()
{
  //evaluate wave vectors
	for (int k=0;k<n3/2+1;++k)
	{
		for (int j=0;j<n2/2+1;++j)
		{
			for (int i=0;i<n1/2+1;++i)
				wv(i,j,k)=cat::tvector<Real,3>(ar1*i,ar2*j,ar3*k);
			for (int i=n1/2+1;i<n1;++i)
				wv(i,j,k)=cat::tvector<Real,3>(ar1*(i-n1),ar2*j,ar3*k);
		}
		for (int j=n2/2+1;j<n2;++j)
		{
			for (int i=0;i<n1/2+1;++i)
				wv(i,j,k)=cat::tvector<Real,3>(ar1*i,ar2*(j-n2),ar3*k);
			for (int i=n1/2+1;i<n1;++i)
				wv(i,j,k)=cat::tvector<Real,3>(ar1*(i-n1),ar2*(j-n2),ar3*k);
		}
	}
	
	wv(0,0,0)=1e-30;
	
      //evaluate the square of norm of wavevectors
	wv2=norm_sq(wv);
	
	wnmax=max(wv2);
	int mx12=n1>n2?n1/2+1:n2/2+1;
	int mx=mx12>n3?mx12:n3;
	int mn12=n1<n2?n1/2+1:n2/2+1;
	int mn=mn12<n3?mn12:n3;
	nwn=mn;
	wnstep=wnmax/(nwn-1);
	
	dealiasing_limit=4./9.*max(wv2);
	
  //initialise dealiasing mask
	cat::array<bool,3>::iterator dm_iter(dealiasing_mask);
	RSF::const_iterator wv2_iter(wv2);
	for(dm_iter=dealiasing_mask.begin(),wv2_iter=wv2.begin();
	    dm_iter!=dealiasing_mask.end(),wv2_iter!=wv2.end();
	    ++dm_iter,++wv2_iter)
	{
		if((*wv2_iter)>dealiasing_limit)
			(*dm_iter)=0;
		else
			(*dm_iter)=1;
	}
}


  //Desctructor
SpectralFourier::~SpectralFourier()
{
}


 //deivatives in respect to coordinate x_index - scalars
CSF SpectralFourier::d_dx_index_hat(const CSF & field,const int index)
{
	CSF out(field.shape());
	out=field;
	out*=wv[index];
	out*=I;
	return out;
}
  //deivatives in respect to coordinate x_index - vectors
CVF SpectralFourier::d_dx_index_hat(const CVF & field, const int index)
{
	CVF out;
	out[0]=this->d_dx_index_hat(field[0],index);
	out[1]=this->d_dx_index_hat(field[1],index);
	out[2]=this->d_dx_index_hat(field[2],index);
	return out;
}


  //Gradient of scalar field
CVF SpectralFourier::grad_hat(const CSF & field)
{
	CVF out(field.shape());
	out=cat::tvector<Complex,3>(I,I,I);
	out*=wv;
	out*=field;
	return out;
}
  //Divergence of vector field
CSF SpectralFourier::div_hat(const CVF & field)
{
	CVF outc(field.shape());
	outc=cat::tvector<Complex,3>(I,I,I);
	outc*=wv;
	return CSF(dot_product(outc,field));
}

  //Curl of vector field
CVF SpectralFourier::curl_hat(const CVF & field)
{
	CVF outc(field.shape());
    //Note that the z derivative acts on the 1st 2 components
	outc=cat::tvector<Complex,3>(I,I,I);
	outc*=wv;
	return CVF(cross_product(outc,field));
}



  //Extract gradient part
  //3D vector field
CVF SpectralFourier::remove_gradient(CVF & field)
{
	CVF out(field.shape());
	CSF sout(field.shape());
	//sout=div_hat(field);
	wv2(0,0,0)=1;
	sout=div_hat(field)/(-wv2);
	wv2(0,0,0)=1e-30;
	out=grad_hat(sout);
	field-=out;
	return out;
}



  //L2 scalar product
  //scalars
Real SpectralFourier::scalar_prod(const CSF & x,
                                  const CSF & y) const
{
	Real out=sum(real(x*conj(y)));
	out*=2.;
	for(int i=0;i<n1;++i)
		for(int j=0;j<n2;++j)
			out-=(x(i,j,0)*conj(y(i,j,0))).real();
	return out;
}
  //vectors
Real SpectralFourier::scalar_prod(const CVF & x,
                                  const CVF & y) const
{
	return Real(scalar_prod(x[0],y[0])+
	            scalar_prod(x[1],y[1])+
	            scalar_prod(x[2],y[2]));
}


//Evaluate energy spectrum
//dividing the sphere in npoints shells 
//scalar fields
cat::array<Real,1> SpectralFourier::eval_energ_spec(const CSF & field)
{
	cat::array<Real,1> out(nwn-1);
	out=0;
	CSF::const_iterator field_iterator(field);
	RSF::iterator wv2_iterator(wv2);
	for(field_iterator=field.begin(),
	    wv2_iterator=wv2.begin();
	    field_iterator!=field.end(),
	    wv2_iterator!=wv2.end();
	    ++field_iterator,
	    ++wv2_iterator)
		out(static_cast<int>(sqrt(*wv2_iterator)/wnstep))+=
		((*field_iterator)*conj(*field_iterator)).real();
	return out;
}
//vector fields
cat::array<Real,1> SpectralFourier::eval_energ_spec(const CVF & field)
{
	cat::array<Real,1> out(nwn-1);
	out=0;
	CVF::const_iterator field_iterator(field);
	RSF::iterator wv2_iterator(wv2);
	for(field_iterator=field.begin(),
	    wv2_iterator=wv2.begin();
	    field_iterator!=field.end(),
	    wv2_iterator!=wv2.end();
	    ++field_iterator,
	    ++wv2_iterator)
		out(static_cast<int>(sqrt(*wv2_iterator)/wnstep))+=
		norm_sq(*field_iterator);
	return out;
}


//Printing non-vanishing harmonics in a field

void SpectralFourier::pnvh(const RSF & field)
{
	CSF field_hat(field.shape());
	sfft.direct_transform(field_hat,field);
	SpectralFourierBase::pnvh_hat(field_hat);
}

void SpectralFourier::pnvh(const RVF & field)
{
	CVF field_hat(field.shape());
	fft.direct_transform(field_hat,field);
	SpectralFourierBase::pnvh_hat(field_hat);
}
