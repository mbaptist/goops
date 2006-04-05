
#include "spectral_fourier_layer.h"

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
spectral_fourier_layer::spectral_fourier_layer(const int & n1__,
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
  wv(n1,n2/2+1,n3),
  wv2(n1,n2/2+1,n3),
  dealiasing_limit(0),
  dealiasing_mask(n1,n2/2+1,n3),
  fft_ccs(cat::tvector<int,3>(n1,n2,n3)),
  fft_ssc(cat::tvector<int,3>(n1,n2,n3)),
  sfft_s(cat::tvector<int,3>(n1,n2,n3)),
  sfft_c(cat::tvector<int,3>(n1,n2,n3))
{

  //define aspect ratios
  Real ar1=2*M_PI/l1;
  Real ar2=2*M_PI/l2;
  Real ar3=M_PI/l3;

  //evaluate wave vectors
  for (int k=0;k<n3;++k)
    for (int j=0;j<n2/2+1;++j)
      {
        for (int i=0;i<n1/2+1;++i)
          wv(i,j,k)=cat::tvector<Real,3>(ar1*i,ar2*j,ar3*k);
        for (int i=n1/2+1;i<n1;++i)
          wv(i,j,k)=cat::tvector<Real,3>(ar1*(i-n1),ar2*j,ar3*k);
      }
  wv(0,0,0)=1e-30;

  //evaluate the square of norm of wavevectors
  wv2=norm_sq(wv);

  dealiasing_limit=4./9.*max(wv2);
  //dealiasing_limit=4./9.*pow(max(wv[2]),2.);
  //dealiasing_limit=4./9.*pow(max(wv[0]),2.);


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
      (*dealiasing_mask_iterator)=0;
    else
      (*dealiasing_mask_iterator)=1;

}


//Desctructor
spectral_fourier_layer::~spectral_fourier_layer()
{}


//Dealiasing
//performs dealiasing for second order non-linearities
//scalar fields
void spectral_fourier_layer::dealias(CSF& field) const
{
  //dealiasing
  field*=dealiasing_mask;
}
//vector fields
void spectral_fourier_layer::dealias(CVF& field) const
{
  //dealiasing
  field*=dealiasing_mask;
}




//Solve lap(f)=g in fourier space
CSF spectral_fourier_layer::poisson_hat(const CSF & field)
{
  CSF out(field.shape());
  out=field;
  wv2(0,0,0)=1;
  out/=(-wv2);
  wv2(0,0,0)=1e-30;
  return out;
}

//Solve lap(f)=g in fourier space
CVF spectral_fourier_layer::poisson_hat(const CVF & field)
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
CSF spectral_fourier_layer::d_dhorizontal_hat(const CSF & field,const int index)
{
  assert(index==0||index==1);
  CSF out(field.shape());
  out=field;
  out*=wv[index];
  out*=I;
  return out;
}
//deivatives in horizontal directions - vectors
CVF spectral_fourier_layer::d_dhorizontal_hat(const CVF & field, const int index)
{
  CVF out(field.shape());
  out[0]=d_dhorizontal_hat(field[0],index);
  out[1]=d_dhorizontal_hat(field[1],index);
  out[2]=d_dhorizontal_hat(field[2],index);
  return out;
}



//Gradient of scalar field
CVF spectral_fourier_layer::grad_hat(const CSF & field,const bool kind)
{
//   CVF out(field.shape());
//   out=cat::tvector<Complex,3>(I,I,(kind?-1:1));
//   out*=wv;
//   out*=field;
//   return out;

  return CVF(cat::tvector<Complex,3>(I,I,(kind?-1:1))*wv*field);

}
//Divergence of vector field
CSF spectral_fourier_layer::div_hat(const CVF & field,const bool kind)
{
//   CVF outc(field.shape());
//   outc=cat::tvector<Complex,3>(I,I,(kind?-1:1));
//   outc*=wv;
//   return CSF(dot_product(outc,field));

  return CSF(dot_product(cat::tvector<Complex,3>(I,I,(kind?-1:1))*wv,field));  

}

//Curl of vector field
CVF spectral_fourier_layer::curl_hat(const CVF & field,const bool kind)
{
//   CVF outc(field.shape());
//   //Note that the z derivative acts on the 1st 2 components
//   outc=cat::tvector<Complex,3>(I,I,(kind?1:-1));
//   outc*=wv;
//   return CVF(cross_product(outc,field));
  
  //Note that the z derivative acts on the 1st 2 components
  return CVF(cross_product(cat::tvector<Complex,3>(I,I,(kind?1:-1))*wv,field));
}


//Laplacian of scalar field
CSF spectral_fourier_layer::lap_hat(const CSF & field)
{
//   CSF out(field.shape());
//   cat::array_iterator<Complex,3> out_iterator(out);
//   cat::array_const_iterator<Complex,3> field_iterator(field);
//   cat::array_const_iterator<Real,3> wv2_iterator(wv2);
//   for(out_iterator=out.begin(),
// 	field_iterator=field.begin(),
// 	wv2_iterator=wv2.begin();
//       out_iterator!=out.end(),
// 	field_iterator!=field.end(),
// 	wv2_iterator!=wv2.end();
//       ++out_iterator,
// 	++field_iterator,
// 	++wv2_iterator)
//     {
//       (*out_iterator)=(*field_iterator);
//       (*out_iterator)*=(-(*wv2_iterator));
//     }


// #if 0
//   CSF out(field.shape());
//   for (int i=0;i<field.shape()[0];++i)
//     for (int j=0;j<field.shape()[1];++j)
//       for (int k=0;k<field.shape()[2];++k)
// 	{
// 	  out(i,j,k)=field(i,j,k);
// 	  out(i,j,k)*=(-wv2(i,j,k));
// 	}
// #endif

//   return out;

  return CSF(-wv2*field);

}
//Laplacian of vector field
CVF spectral_fourier_layer::lap_hat(const CVF & field)
{
//   CVF out(field.shape());
//   cat::array_iterator<cat::tvector<Complex,3>,3> out_iterator(out);
//   cat::array_const_iterator<cat::tvector<Complex,3>,3> field_iterator(field);
//   cat::array_const_iterator<Real,3> wv2_iterator(wv2);
//   for(out_iterator=out.begin(),
// 	field_iterator=field.begin(),
// 	wv2_iterator=wv2.begin();
//       out_iterator!=out.end(),
// 	field_iterator!=field.end(),
// 	wv2_iterator!=wv2.end();
//       ++out_iterator,
// 	++field_iterator,
// 	++wv2_iterator)
//     {
//       (*out_iterator)=(*field_iterator);
//       (*out_iterator)*=(-(*wv2_iterator));
//     }



// #if 0
//   CVF out(field.shape());
//   for (int i=0;i<field.shape()[0];++i)
//     for (int j=0;j<field.shape()[1];++j)
//       for (int k=0;k<field.shape()[2];++k)
//         for(int m=0;m<3;++m)
//           {
//             out(i,j,k)[m]=field(i,j,k)[m];
//             out(i,j,k)[m]*=(-wv2(i,j,k));
//           }
// #endif

//   return out;

return CVF(-wv2*field);

}


//Extract gradient part 
//3D vector field
CVF spectral_fourier_layer::remove_gradient(CVF & field,const bool kind)
{
//   CVF out(field.shape());
//   CSF sout(field.shape());
//   sout=div_hat(field,kind);
//   wv2(0,0,0)=1;
//   sout/=(-wv2);
//   wv2(0,0,0)=1e-30;
//   out=grad_hat(sout,!kind);
//   field-=out;
//   return out;

  CVF out(field.shape());
  out=grad_hat(-div_hat(field,kind)/wv2,!kind);
  field-=out;
  return out;

}



//L2 scalar product
//scalars
Real spectral_fourier_layer::scalar_prod(const CSF & x,
                                   const CSF & y,
				   const bool kind) const
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

  Real out=0;
  for(int i=0;i<n1;++i)
    for(int k=1;k<n3;++k)
      out+=(x(i,0,k)*conj(y(i,0,k))).real();
  out*=-.5;
  out+=sum(real(x*conj(y)));
  return out;
#endif

}
//scalars kind s
Real spectral_fourier_layer::scalar_prod(const CSF & x,
                                   const CSF & y) const
{
  return Real(scalar_prod(x,y,0));
}
//vectors kind ccs
Real spectral_fourier_layer::scalar_prod(const CVF & x,
                                   const CVF & y) const
{
  return Real(scalar_prod(x[0],y[0],1)+
		scalar_prod(x[1],y[1],1)+
		scalar_prod(x[2],y[2],0));
}



//Evaluate energy spectrum
//scalar fields
cat::array<Real,1> spectral_fourier_layer::eval_energ_spec(const CSF & field)
{
  int mx12=n1>n2?n1:n2;
  int mx=mx12>n3?mx12:n3;
  return eval_energ_spec(field,mx);
}
//vector fields
cat::array<Real,1> spectral_fourier_layer::eval_energ_spec(const CVF & field)
{
  int mx12=n1>n2?n1:n2;
  int mx=mx12>n3?mx12:n3;
  return eval_energ_spec(field,mx);
}
//Evaluate energy spectrum
//dividing the sphere in npoints shells 
//scalar fields
cat::array<Real,1> spectral_fourier_layer::eval_energ_spec(const CSF & field,const int & npoints)
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
cat::array<Real,1> spectral_fourier_layer::eval_energ_spec(const CVF & field,const int & npoints)
{
  cat::array<Real,1> out(npoints);
  out=0;
  Real wv2step=max(wv2)/(npoints-1);
  cat::array_const_iterator<cat::tvector<Complex,3>,3> 
    field_iterator(field);
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
void spectral_fourier_layer::pnvh(const CSF & field)
{
  for(int i=0;i<field.shape()[0];i++)
    for(int j=0;j<field.shape()[1];j++)
      for(int k=0;k<field.shape()[2];k++)
        if(abs(field(i,j,k))>1e-10)
          cout << i << " " << j << " " << k << " " << setprecision(20) << field(i,j,k) << endl;
}
//vector
void spectral_fourier_layer::pnvh(const CVF & field)
{
  for(int i=0;i<field.shape()[0];i++)
    for(int j=0;j<field.shape()[1];j++)
      for(int k=0;k<field.shape()[2];k++)
        if(norm(field(i,j,k))>1e-10)
          cout << i << " " << j << " " << k << " " << setprecision(20) << field(i,j,k) << endl;
}


void spectral_fourier_layer::pnvh(const RSF & field)
{
  CSF field_hat(field.shape());
  sfft_s.direct_transform(field_hat,field);
  pnvh(field_hat);
}

void spectral_fourier_layer::pnvh(const RVF & field)
{
  CVF field_hat(field.shape());
  fft_ccs.direct_transform(field_hat,field);
  pnvh(field_hat);  
}
