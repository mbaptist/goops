
#include "spectral_fourier_layer.h"

#include "globals.h"

#include "../goops_types.h"
#include "../fft/fft.h"

#include<cmath>
#include<complex>

#include <cat.h>

using namespace std;
using namespace cat;

//Construtor
SpectralFourierLayer::SpectralFourierLayer(const int & n1__,const int & n2__,const int & n3__,
                                           const Real & l1__,const Real & l2__,const Real & l3__):
SpectralFourierBase(n1__,n2__,n3__,
                    n1__,n2__/2+1,n3__,
                    l1__,l2__,l3__,
                    2*M_PI/l1__,2*M_PI/l2__,M_PI/l3__),
fft_ccs("cos","cos","sin"),
fft_ssc("sin","sin","cos"),
sfft_s("sin"),
sfft_c("cos")
{
	initialise();
}

void SpectralFourierLayer::initialise()
{
	
//evaluate wave vectors
	for (int k=0;k<n3_hat;++k)
		for (int j=0;j<n2_hat;++j)
		{
			for (int i=0;i<n1_hat/2+1;++i)
				wv(i,j,k)=cat::tvector<Real,3>(ar1*i,ar2*j,ar3*k);
			for (int i=n1_hat/2+1;i<n1_hat;++i)
				wv(i,j,k)=cat::tvector<Real,3>(ar1*(i-n1_hat),ar2*j,ar3*k);
		}
	wv(0,0,0)=1e-30;
	
  //evaluate the square of norm of wavevectors
	wv2=norm_sq(wv);

	wnmax=sqrt(max(wv2));
	int mx12=n1>n2?(n1/2+1):(n2/2+1);
	int mx=mx12>n3?mx12:n3;
	int mn12=n1<n2?(n1/2+1):(n2/2+1);
	int mn=mn12<n3?mn12:n3;
	nwn=mx;
	wnstep=wnmax/(nwn-1);
	
	//Set dealiasing limit
	Real wvmx0=max(wv[0]);
	Real wvmx1=max(wv[1]);
	Real wvmx2=max(wv[2]);
  Real min_wvmx01=(wvmx0<wvmx1?wvmx0:wvmx1);
  Real min_wvmx=	(min_wvmx01<wvmx2?min_wvmx01:wvmx2);
	dealiasing_limit=4./9.*pow(min_wvmx,2);
	
	
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
SpectralFourierLayer::~SpectralFourierLayer()
{
}

//Derivatives in respect to coordinate x_index - scalars
CSF SpectralFourierLayer::d_dx_index_hat(const CSF & field,const int index)
{
	assert(index==0||index==1);
	return CSF(I*wv[index]*field);
}
  //deivatives in respect to coordinate x_index - vectors
CVF SpectralFourierLayer::d_dx_index_hat(const CVF & field, const int index)
{
	assert(index==0||index==1);
	CVF out(field.shape());
	out[0]=this->d_dx_index_hat(field[0],index);
	out[1]=this->d_dx_index_hat(field[1],index);
	out[2]=this->d_dx_index_hat(field[2],index);
	return out;
}

//Gradient of scalar field
CVF SpectralFourierLayer::grad_hat(const CSF & field,const bool kind)
{
return CVF(cat::tvector<Complex,3>(I,I,(kind?-1:1))*wv*field);
}

//Divergence of vector field
CSF SpectralFourierLayer::div_hat(const CVF & field,const bool kind)
{
return CSF(dot_product(cat::tvector<Complex,3>(I,I,(kind?-1:1))*wv,field));
}

//Curl of vector field
CVF SpectralFourierLayer::curl_hat(const CVF & field,const bool kind)
{
  //Note that the z derivative acts on the 1st 2 components
return CVF(cross_product(cat::tvector<Complex,3>(I,I,(kind?1:-1))*wv,field));
}

//Extract gradient part 
//3D vector field
CVF SpectralFourierLayer::remove_gradient(CVF & field,const bool kind)
{
	CVF out(field.shape());
	wv2(0,0,0)=1;
	out=grad_hat(-div_hat(field,kind)/wv2,!kind);
	wv2(0,0,0)=1e-30;
	field-=out;
	return out;
}

//L2 scalar product
//scalars
Real SpectralFourierLayer::scalar_prod(const CSF & x,const CSF & y,const bool kind) const
{
#if 1
	Real zerokz=0;
	Real otherkz=0;
	
	if (kind==1)
	{
		for(int i=1;i<n1/2+1;++i)
		{
			zerokz+=real(x(i,0,0)*conj(y(i,0,0)));
		}
		for(int i=0;i<n1;++i)
			for(int j=1;j<n2/2+1;++j)
				zerokz+=real(x(i,j,0)*conj(y(i,j,0)));
		zerokz*=2;
		zerokz+=real(x(0,0,0)*conj(y(0,0,0)));
	}
	
	for(int k=1;k<n3;++k)
		otherkz+=real(x(0,0,k)*conj(y(0,0,k)));
	otherkz*=.5;
	for(int i=1;i<n1/2+1;++i)
		for(int k=1;k<n3;++k)
			otherkz+=real(x(i,0,k)*conj(y(i,0,k)));
	for(int i=0;i<n1;++i)
		for(int j=1;j<n2/2+1;++j)
			for(int k=1;k<n3;++k)
				otherkz+=real(x(i,j,k)*conj(y(i,j,k)));
	return Real(zerokz+otherkz);
#endif
#if 0
	Real zerokz=0;
	Real otherkz=0;

	if (kind==1)//if cosine
	{
		for(int i=0;i<n1;++i)
			for(int j=1;j<n2/2+1;++j)
					zerokz+=real(x(i,j,0)*conj(y(i,j,0)));
		zerokz*=2;
		for(int i=0;i<n1;++i)
				zerokz+=real(x(i,0,0)*conj(y(i,0,0)));
	}
	//sin and cosine
	for(int i=0;i<n1;++i)
		for(int k=1;k<n3;++k)
			otherkz+=real(x(i,0,k)*conj(y(i,0,k)));
	otherkz*=.5;
	for(int i=0;i<n1;++i)
		for(int j=1;j<n2/2+1;++j)
			for(int k=1;k<n3;++k)
				otherkz+=real(x(i,j,k)*conj(y(i,j,k)));
	return Real(zerokz+otherkz);
	
#endif
#if 0
	Real out;
	out=0;
	CSF::const_iterator x_iter(x);
	CSF::const_iterator y_iter(y);
	RSF::const_iterator wv2_iter(wv2);
	for(x_iter=x.begin(),y_iter=y.begin(),wv2_iter=wv2.begin();
	    x_iter!=x.end();
	    ++x_iter,++y_iter,++wv2_iter)
	{
		if((x_iter.indices())[2]!=0)
		{
			if((x_iter.indices())[1]==0)
				out+=.5*((*x_iter)*conj(*y_iter)).real();
			else
				out+=((*x_iter)*conj(*y_iter)).real();
		}
		else if(kind==1)
		{
				if((x_iter.indices())[1]==0)
					out+=((*x_iter)*conj(*y_iter)).real();
				else
					out+=2.*((*x_iter)*conj(*y_iter)).real();
			}
		}
	
	return out;
#endif
	
}
//scalars kind s
Real SpectralFourierLayer::scalar_prod(const CSF & x,const CSF & y) const
{
	return Real(scalar_prod(x,y,0));
}

//vectors
Real SpectralFourierLayer::scalar_prod(const CVF & x,const CVF & y,const bool kind) const
{
	return Real(scalar_prod(x[0],y[0],!kind)+scalar_prod(x[1],y[1],!kind)+scalar_prod(x[2],y[2],kind));
}
//vectors kind ccs
Real SpectralFourierLayer::scalar_prod(const CVF & x,const CVF & y) const
{
	return Real(scalar_prod(x[0],y[0],1)+scalar_prod(x[1],y[1],1)+scalar_prod(x[2],y[2],0));
}

//Evaluate energy spectrum
//dividing the sphere in nwn-1 shells
//scalar fields
cat::array<Real,1> SpectralFourierLayer::eval_energ_spec(const CSF & field,const bool & kind)
{
	cat::array<Real,1> out(nwn);
	out=0;	
	CSF::const_iterator field_iterator(field);
	RSF::iterator wv2_iterator(wv2);
	for(field_iterator=field.begin(),
	    wv2_iterator=wv2.begin();
	    field_iterator!=field.end(),
	    wv2_iterator!=wv2.end();
	    ++field_iterator,
	    ++wv2_iterator)
	{
		double mf=1.;
		if((field_iterator.indices())[2]!=0)
			mf*=1.;
		else if(kind==1)
			mf*=2.;
		else if(kind==0)
			mf*=0.;
		int index=static_cast<int>(sqrt(*wv2_iterator)/wnstep);
		//cout << index << endl;
		if((field_iterator.indices())[1]==0)
			out(index)+=mf*.5*((*field_iterator)*conj(*field_iterator)).real();
		else
			out(index)+=mf*((*field_iterator)*conj(*field_iterator)).real();
	}
	out*=.5/wnstep;
	return out;
}
//vector fields
cat::array<Real,1> SpectralFourierLayer::eval_energ_spec(const CVF & field,const bool & kind)
{
	return cat::array<Real,1>(eval_energ_spec(field[0],!kind)+eval_energ_spec(field[1],!kind)+eval_energ_spec(field[2],kind));
}


//Printing non-vanishing harmonics in a field

void SpectralFourierLayer::pnvh(const RSF & field,const bool & kind)
{
	CSF field_hat(field.shape());
	if(kind==0)
		sfft_s.direct_transform(field_hat,field);
	else
		sfft_c.direct_transform(field_hat,field);
	SpectralFourierBase::pnvh_hat(field_hat);
}

void SpectralFourierLayer::pnvh(const RVF & field,const bool & kind)
{
	CVF field_hat(field.shape());
	if (kind==0)
		fft_ccs.direct_transform(field_hat,field);
	else
		fft_ssc.direct_transform(field_hat,field);
	SpectralFourierBase::pnvh_hat(field_hat);
}

