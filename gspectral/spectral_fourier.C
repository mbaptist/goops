
#include "spectral_fourier.h"

#include "globals.h"

#include "../goops_types.h"
#include "../gffti/gffti.h"

#include <cat.h>

#include <iomanip>

#include <complex>
#include <fstream>

using namespace std;
using namespace cat;

//Construtor
spectral_fourier::spectral_fourier(const int & n1__,
				   const int & n2__,
				   const int & n3__,
				   const Real & l1__,
				   const Real & l2__,
				   const Real & l3__):
  n1(n1__),
  n2(n2__),
  n3(n3__),
  l1(l1__),
  l2(l2__),
  l3(l3__),
  wv(n1,n2,n3/2+1),
  wv2(n1,n2,n3/2+1),
  dealiasing_limit(0),
  dealiasing_mask(n1,n2,n3/2+1),
  fft(cat::tvector<int,3>(n1,n2,n3)),
  sfft()
{

  //define aspect ratios
  Real ar1=2*M_PI/l1;
  Real ar2=2*M_PI/l2;
  Real ar3=2*M_PI/l3;
  
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

      dealiasing_limit=4./9.*max(wv2);

      //initialise dealiasing mask
      cat::array_iterator<bool,3> dealiasing_mask_iterator(dealiasing_mask);
      cat::array_const_iterator<Real,3> wv2_iterator(wv2);
      for(dealiasing_mask_iterator=dealiasing_mask.begin(),
	    wv2_iterator=wv2.begin();
	  dealiasing_mask_iterator!=dealiasing_mask.end(),
	    wv2_iterator!=wv2.end();
	  ++dealiasing_mask_iterator,
	    ++wv2_iterator)
	if((*wv2_iterator)>dealiasing_limit)
	  *dealiasing_mask_iterator=0;
	else
	  *dealiasing_mask_iterator=1;

}


  //Desctructor
  spectral_fourier::~spectral_fourier()
    {
    }


  //Dealiasing
  //performs dealiasing for second order non-linearities
  //scalar fields
  void spectral_fourier::dealias(CSF& field) const
  {
    field*=dealiasing_mask;
  }
  //vector fields
  void spectral_fourier::dealias(CVF& field) const
  {
    field*=dealiasing_mask;
  }


  //Solve lap(f)=g in fourier space
  CSF spectral_fourier::poisson_hat(const CSF & field)
  {
    CSF out(field.shape());
    out=field;
    wv2(0,0,0)=1;
    out/=(-wv2);
    wv2(0,0,0)=1e-30;
    return out;
  }

  //Solve lap(f)=g in fourier space
  CVF spectral_fourier::poisson_hat(const CVF & field)
  {
    CVF out(field.shape());
    out=field;
    wv2(0,0,0)=1;
    out/=(-wv2);
    wv2(0,0,0)=1e-30;
    return out;
  }


  //Derivatives in fourier space 


  //Differential operators in fourier space
  //deivatives in horizontal directions - scalars
  CSF spectral_fourier::d_dx_index_hat(const CSF & field,const int index)
  {
    CSF out(field.shape());
    out=field;
    out*=wv[index];
    out*=I;
    return out;
  }
  //deivatives in horizontal directions - vectors
  CVF spectral_fourier::d_dx_index_hat(const CVF & field, const int index)
  {
    CVF out;
    out[0]=d_dx_index_hat(field[0],index);
    out[1]=d_dx_index_hat(field[1],index);
    out[2]=d_dx_index_hat(field[2],index);
    return out;
  }



  //Gradient of scalar field
  CVF spectral_fourier::grad_hat(const CSF & field)
  {
    CVF out(field.shape());
    out=cat::tvector<Complex,3>(I,I,I);
    out*=wv;
    out*=field;
    return out;
  }
  //Divergence of vector field
  CSF spectral_fourier::div_hat(const CVF & field)
  {
    CVF outc(field.shape());
    outc=cat::tvector<Complex,3>(I,I,I);
    outc*=wv;
    return CSF(dot_product(outc,field));
  }

  //Curl of vector field
  CVF spectral_fourier::curl_hat(const CVF & field)
  {
    CVF outc(field.shape());
    //Note that the z derivative acts on the 1st 2 components
    outc=cat::tvector<Complex,3>(I,I,I);
    outc*=wv;
    return CVF(cross_product(outc,field));
  }


  //Laplacian of scalar field
  CSF spectral_fourier::lap_hat(const CSF & field)
  {
    CSF out(field.shape());
    cat::array_iterator<Complex,3> out_iterator(out);
    cat::array_const_iterator<Complex,3> field_iterator(field);
    cat::array_const_iterator<Real,3> wv2_iterator(wv2);
    for(out_iterator=out.begin(),
	  field_iterator=field.begin(),
	  wv2_iterator=wv2.begin();
	out_iterator!=out.end(),
	  field_iterator!=field.end(),
	  wv2_iterator!=wv2.end();
	++out_iterator,
	  ++field_iterator,
	  ++wv2_iterator)
      {
	(*out_iterator)=(*field_iterator);
	(*out_iterator)*=(-(*wv2_iterator));
      }
    return out;
  }
  //Laplacian of vector field
  CVF spectral_fourier::lap_hat(const CVF & field)
  {
    CVF out(field.shape());
    cat::array_iterator<cat::tvector<Complex,3>,3> out_iterator(out);
    cat::array_const_iterator<cat::tvector<Complex,3>,3> field_iterator(field);
    cat::array_const_iterator<Real,3> wv2_iterator(wv2);
    for(out_iterator=out.begin(),
	  field_iterator=field.begin(),
	  wv2_iterator=wv2.begin();
	out_iterator!=out.end(),
	  field_iterator!=field.end(),
	  wv2_iterator!=wv2.end();
	++out_iterator,
	  ++field_iterator,
	  ++wv2_iterator)
      {
	(*out_iterator)=(*field_iterator);
	(*out_iterator)*=(-(*wv2_iterator));
      }
    return out;
  }

  //Extract gradient part 
  //3D vector field
  CVF spectral_fourier::remove_gradient(CVF & field)
  {
    CVF out(field.shape());
    CSF sout(field.shape());
    sout=div_hat(field);
    wv2(0,0,0)=1;
    sout/=(-wv2);
    wv2(0,0,0)=1e-30;
    out=grad_hat(sout);
    field-=out;
    return out;
  }



  //L2 scalar product
  //scalars
  Real spectral_fourier::scalar_prod(const CSF & x,
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
  Real spectral_fourier::scalar_prod(const CVF & x,
				       const CVF & y) const
  {
    return Real(scalar_prod(x[0],y[0])+
		  scalar_prod(x[1],y[1])+
		  scalar_prod(x[2],y[2]));
  }


//Evaluate energy spectrum
//scalar fields
cat::array<Real,1> spectral_fourier::eval_energ_spec(const CSF & field)
{
  int mx12=n1>n2?n1:n2;
  int mx=mx12>n3?mx12:n3;
  return eval_energ_spec(field,mx);
}
//vector fields
cat::array<Real,1> spectral_fourier::eval_energ_spec(const CVF & field)
{
  int mx12=n1>n2?n1:n2;
  int mx=mx12>n3?mx12:n3;
  return eval_energ_spec(field,mx);
}
//Evaluate energy spectrum
//dividing the sphere in npoints shells 
//scalar fields
cat::array<Real,1> spectral_fourier::eval_energ_spec(const CSF & field,const int & npoints)
{
  cat::array<Real,1> out(npoints);
  out=0;
  Real wv2step=max(wv2)/(npoints-1);
  cat::array_const_iterator<Complex,3> field_iterator(field);
  cat::array_iterator<Real,3> wv2_iterator(wv2);
  for(field_iterator=field.begin(),
	wv2_iterator=wv2.begin();
      field_iterator!=field.end(),
	wv2_iterator!=wv2.end();
      ++field_iterator,
	++wv2_iterator)
    out(static_cast<int>((*wv2_iterator)/wv2step))+=
      ((*field_iterator)*conj(*field_iterator)).real();
  return out;
}
//vector fields
cat::array<Real,1> spectral_fourier::eval_energ_spec(const CVF & field,const int & npoints)
{
  cat::array<Real,1> out(npoints);
  out=0;
  Real wv2step=max(wv2)/(npoints-1);
  cat::array_const_iterator<cat::tvector<Complex,3>,3> field_iterator(field);
  cat::array_iterator<Real,3> wv2_iterator(wv2);
  for(field_iterator=field.begin(),
	wv2_iterator=wv2.begin();
      field_iterator!=field.end(),
	wv2_iterator!=wv2.end();
      ++field_iterator,
	++wv2_iterator)
    out(static_cast<int>((*wv2_iterator)/wv2step))+=
      norm_sq(*field_iterator);
  return out;
}


//Printing non-vanishing harmonics in a field
//scalar
void spectral_fourier::pnvh(const CSF & field)
{
    for(int i=0;i<field.shape()[0];i++)
      for(int j=0;j<field.shape()[1];j++)
	for(int k=0;k<field.shape()[2];k++)
	  if(abs(field(i,j,k))>1e-10)
	    cout << i << " " << j << " " << k << " " << setprecision(20) << field(i,j,k) << endl;
  }
  //vector
  void spectral_fourier::pnvh(const CVF & field)
  {
    for(int i=0;i<field.shape()[0];i++)
      for(int j=0;j<field.shape()[1];j++)
	for(int k=0;k<field.shape()[2];k++)
	  if(norm(field(i,j,k))>1e-10)
	    cout << i << " " << j << " " << k << " " << setprecision(20) << field(i,j,k) << endl;
  }


  void spectral_fourier::pnvh(const RSF & field)
  {
    CSF field_hat(field.shape());
    sfft.direct_transform(field_hat,field);
    pnvh(field_hat);
  }

  void spectral_fourier::pnvh(const RVF & field)
  {
    CVF field_hat(field.shape());
    fft.direct_transform(field_hat,field);
    pnvh(field_hat);  
  }
