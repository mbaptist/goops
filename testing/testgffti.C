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


#include "../goops.h"

#include <cat.h>

#include <iostream>
#include <complex>

using namespace cat;
using namespace std;



int main()
{
  int n1,n2,n3;
  n1=8;
  n2=n1;
  n3=n2;

#if 0

  array<tvector<int,2>,2> fint(n1,n2);
  for(int i=0;i<n1;++i)
    for(int j=0;j<n2;++j)
      fint(i,j)=tvector<int,2>(i,j);

  cout << fint << endl;
  exit (0);
#endif

#if 0

  array<double,2> field2(n1,n2);
  array<complex<double>,2> field2_hat(n1,n2/2+1);

  for (int i=0;i<n1;++i)
    for (int j=0;j<n2;++j)
      field2(i,j)=sin(2*M_PI/n2*j);

  s_rfft<2> sfft2(tvector<int,2>(n1,n2));
  sfft2.direct_transform(field2_hat,field2);

#if 1

  for(int i=0;i<n1;++i)
    for(int j=0;j<n2/2+1;++j)
      {
        if (norm(field2_hat(i,j))>1e-10)
          cout << i << " " << j << " " << field2_hat(i,j) << endl;
      }
#endif
#endif

#if 0
  array<tvector<int,3>,3> fint3(n1,n2,n3);
  for(int i=0;i<n1;++i)
    for(int j=0;j<n2;++j)
      for(int k=0;k<n3;++k)
        fint3(i,j,k)=tvector<int,3>(i,j,k);

  cout << fint3 << endl;
  
  for(int i=0;i<n1;++i)
    for(int j=0;j<n2;++j)
      for(int k=0;k<n3;++k)
        cout << i << " " << j << " " << k << " " << fint3(i,j,k) << endl;
  
  exit (0);
#endif


#if 0

  array<double,3> field(n1,n2,n3);
  array<complex<double>,3> field_hat(n1,n2,n3/2+1);

  for (int i=0;i<n1;++i)
    for (int j=0;j<n2;++j)
      for (int k=0;k<n3;++k)
        field(i,j,k)=sin(2*M_PI/n1*i);

  s_rfft<3> sfft(tvector<int,3>(n1,n2,n3));
  sfft.direct_transform(field_hat,field);

#if 1

  for(int i=0;i<n1;++i)
    for(int j=0;j<n2;++j)
      for(int k=0;k<n3/2+1;++k)
        {
          if (norm(field_hat(i,j,k))>1e-10)
            cout << i << " " << j << " " << k << " " << field_hat(i,j,k) << endl;
        }
#endif

#endif


#if 1
int n=15;

array<double,1> field(n);
array<double,1> field_hat(n);

s_sinfft_1d sinfft(n);
s_cosfft_1d cosfft(n);

cout << "cos " << endl;

 field_hat=0;
// field_hat(5)=1;
// field_hat(3)=2;
 field_hat(0)=1;


cosfft.inverse_transform(field,field_hat);

cout << field << endl;

field_hat=0;

cosfft.direct_transform(field_hat,field);

cout << " " << endl;
cout << field_hat << endl;



cout << " " << endl;
cout << "sin " << endl;

field_hat=0;
field_hat(1)=1;

sinfft.inverse_transform(field,field_hat);

cout << field << endl;

field_hat=0;

sinfft.direct_transform(field_hat,field);

cout << " " << endl;
cout << field_hat << endl;
#endif



#if 0

 cout << n1 << n2 << n3 << endl; 

  array<double,3> field(n1,n2,n3);
  array<complex<double>,3> field_hat(n1,n2/2+1,n3);

  field_hat=0;
  field_hat(1,0,0)=1;
  field_hat(7,0,0)=1;

  s_rfft_h_z<s_cosfft_1d> sfft(tvector<int,3>(n1,n2,n3));
  sfft.inverse_transform(field,field_hat);

	field_hat=0;
	sfft.direct_transform(field_hat,field);

  for(int i=0;i<n1;++i)
    for(int j=0;j<n2/2+1;++j)
      for(int k=0;k<n3;++k)
        {
          if (norm(field_hat(i,j,k))>1e-10)
            cout << i << " " << j << " " << k << " " << field_hat(i,j,k) << endl;
        }

#endif





  return(0);
}
