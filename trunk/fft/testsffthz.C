/*

Copyright 2004,2005,2006 Manuel Baptista

This file is part of GOOPS

GOOPS is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

GOOPS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/


#include "fft.h"

#include "../goops.h"

#include <cat.h>

#include <iostream>
#include <complex>

#include<cmath>

using namespace cat;
using namespace std;
//using namespace goops;


int main()
{
    
    
    int n1,n2,n3;
    n1=32;
    n2=32;
    n3=16;
    double l1,l2,l3;
    l1=2.*M_PI;
    l2=2.*M_PI;
    l3=M_PI;
    
    FFT<Array<double,1>,Array<double,1> > fft_s("sin");
    FFT<Array<double,1>,Array<double,1> > fft_c("cos");
    s_rfft_h_z sfft_s("sin");
    s_rfft_h_z sfft_c("cos");
    v_rfft_h_z vfft_ssc("sin","sin","cos");
    v_rfft_h_z vfft_ccs("cos","cos","sin");
    
    
#if 0
    Array<double,1> rf(n1);
    Array<double,1> trf(n1);
    Array<double,1> ff(n1);
    
    for(int i=0;i<n1;++i)
    {
	double x=1.*i*M_PI/(n1-1);
	trf(i)=cos(2*x);
    }
    
    ff=0;
    ff(2)=1;	
    
    fftc.inverse_transform(rf,ff);
    
    cout << rf << endl;
    cout << "\n\n\n" << endl;
    cout << trf << endl;
    
    ff=0;
    
    fftc.direct_transform(ff,trf);
    
    cout << "\n\n\n" << endl;
    cout << ff << endl;
    
    cout << sum(rf-trf) << endl;
    
#endif
    
    
#if 1 
    const double alpha=1;
    
    Array<double,3> sfield(n1,n2,n3);
    Array<Tvector<double,3>,3> gsfield(n1,n2,n3);
    Array<double,3> dgsfield(n1,n2,n3);
    
    for(int i=0;i<n1;++i)
    {
	double x=l1/n1*i;
	for(int j=0;j<n2;++j)
	{
	    double y=l2/n2*j;			
	    for(int k=0;k<n3;++k)
	    {
		double z=l3/(n3-1)*k;
		sfield(i,j,k)=sin(alpha*x+y)*cos(alpha*z);
		gsfield(i,j,k)=Tvector<double,3>(alpha*cos(alpha*x+y)*cos(alpha*z),cos(alpha*x+y)*cos(alpha*z),-alpha*sin(alpha*x+y)*sin(alpha*z));
	    }
	}
    }
    dgsfield=-alpha*(2*alpha+1)*sfield;
    
    
    SpectralFourierLayer so(n1,n2,n3,l1,l2,l3);
    
    Array<complex<double>,3> tsfield_hat(n1,n2/2+1,n3);
    Array<Tvector<complex<double>,3>,3> tgsfield_hat(n1,n2/2+1,n3);
    Array<Tvector<complex<double>,3>,3> tcgsfield_hat(n1,n2/2+1,n3);
    Array<complex<double>,3> tdgsfield_hat(n1,n2/2+1,n3);
    
    Array<double,3> tsfield(n1,n2,n3);
    Array<Tvector<double,3>,3> tgsfield(n1,n2,n3);
    Array<Tvector<double,3>,3> tcgsfield(n1,n2,n3);
    Array<double,3> tdgsfield(n1,n2,n3);
    
    sfft_c.direct_transform(tsfield_hat,sfield);
    tgsfield_hat=so.grad_hat(tsfield_hat,1);
    vfft_ccs.inverse_transform(tgsfield,tgsfield_hat);
    tcgsfield_hat=so.curl_hat(tgsfield_hat,0);
    vfft_ssc.inverse_transform(tcgsfield,tcgsfield_hat);
    tdgsfield_hat=so.div_hat(tgsfield_hat,0);
    sfft_c.inverse_transform(tdgsfield,tdgsfield_hat);
    
    cout << sum(norm_sq(gsfield)-norm_sq(tgsfield)) <<endl;
    cout << sum(norm_sq(tcgsfield)) <<endl;
    cout << sum(dgsfield*dgsfield-tdgsfield*tdgsfield) <<endl;
    
    cout << "\n" << endl;
    
    
    for(int i=0;i<n1;++i)
    {
	double x=l1/n1*i;
	for(int j=0;j<n2;++j)
	{
	    double y=l2/n2*j;			
	    for(int k=0;k<n3;++k)
	    {
		cout << gsfield(i,j,k) << " " << tgsfield(i,j,k) << endl;
	    }
	}
    }
    
    
    
    
    
#endif
    
    return 0;
    
    
}
