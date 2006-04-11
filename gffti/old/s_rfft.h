// -*- C++ -*-

/////////////////////////////////////////////////
// Real TO COMPLEX TRANSFORM FOR SCALAR FIELDS //
/////////////////////////////////////////////////

#ifndef S_RFFT_H
#define S_RFFT_H

#include <complex>
using std::complex;

#include <fftw3.h>

#include <cat.h>
using namespace cat;

#include "fft_base.h"

//CLASS S_RFFT
//DECLARATION
template <int D>
class s_rfft : public real_fft_base<D>
  {
    //Types for scalar fields
    typedef cat::array<Real,D> RT;
    typedef cat::array<complex<Real>,D> CT;
  private:
    using real_fft_base<D>::size;//Size of the vector to be transformed, along each dimension
    using real_fft_base<D>::total_size;//Total size(with last dimension not padded)
    using real_fft_base<D>::sizep;//Size along the padded dimension
    using real_fft_base<D>::pad_size;//Size of padding
    using real_fft_base<D>::total_size_Dm1;//Total size without last dimension
    using real_fft_base<D>::total_sizep;//Total size padded (with last dimension padded)
    //plan variables
    using fft_base<D>::direct_plan;//direct transforms
    using fft_base<D>::inverse_plan;//inverse transforms

  public:
    //Constructor
    s_rfft(cat::tvector<int,D> size_);
    s_rfft(cat::tvector<int,D> size_,int num_threads);
    //Destructor
    ~s_rfft();

    //pointers to local data
    double * real_data;//Pointer to real data
    fftw_complex * complex_data;//Pointer to complex data
    //Transformation functions (with references)
    void direct_transform(CT & u_hat,const RT & u);
    void inverse_transform(RT & u,const CT & u_hat);
    //Transformation functions (with output temporary)
    CT direct_transform(const RT & u);
    RT inverse_transform(const CT & u_hat);

  private:
    
    void initialise();

  };

//IMPLEMENTATION
//constructor
template <int D>
s_rfft<D>::s_rfft(cat::tvector<int,D> size_):
  real_fft_base<D>(size_)
{
  initialise();
}
template <int D>
s_rfft<D>::s_rfft(cat::tvector<int,D> size_,int num_threads):
  real_fft_base<D>(size_,num_threads)
{
  initialise();
}
//Destructor
template <int D>
s_rfft<D>::~s_rfft()
{
  fftw_destroy_plan(direct_plan);
  fftw_destroy_plan(inverse_plan);
  complex_data=0;
  //fftw_free(real_data);
  delete[] real_data;
  delete[] size;
}

template <int D>
void s_rfft<D>::initialise()
{
  //allocate temporary space for inline transform (as a real variable)
  //  real_data=static_cast<Real *>(fftw_malloc(sizeof(Real) * total_sizep));
  real_data=new double[total_sizep];
  //initialize temporary space with zeros
  for(int i=0;i<total_sizep;i++)
    real_data[i]=0;
  //create a complex pointer to real data
  complex_data=reinterpret_cast<fftw_complex *>(real_data);
  //create plan for direct (real-to-complex) transforms
  direct_plan=fftw_plan_dft_r2c(D,size,real_data,complex_data,FFTW_ESTIMATE);
  //create plan for inverse (complex-to-real) transforms
  inverse_plan=fftw_plan_dft_c2r(D,size,complex_data,real_data,FFTW_ESTIMATE);
}

//direct transform (scalar)
template <int D>
void s_rfft<D>::direct_transform(CT& u_hat,const RT& u)
{

#if 0

  int indpos=0;
  for (int i=0;i<size[0];++i)
    for (int j=0;j<size[1];++j)
      for (int k=0;k<size[2];++k)
        {
          cout << u(i,j,k)-u.data()[indpos] << endl;
          ++indpos;
        }
  exit(0);

#endif





  //read real data from input field via pointer
#if 1

  int depad=0;
  int jb=0,je=0;
  for(int i=0;i<total_size_Dm1;++i)
    {
      je=(i+1)*sizep;
      for(int j=jb;j<je-pad_size;++j)
        real_data[j]=(u.data()[j-depad]);
      for(int j=je-pad_size;j<je;++j)
        real_data[j]=0;
      jb=je;
      depad+=pad_size;
    }
#endif

#if 0
  int index_data=0;
  int index_temp=0;
  while(index_temp<total_sizep)
    {

      //cout << index_data << " " << index_temp << endl;

      real_data[index_temp]=(u.data()[index_data]);

      ++index_data;
      ++index_temp;

      //if ((index_data%total_size_Dm1)==0)

      if ((index_data%size[D-1])==size[D-1]-1)
        {
          real_data[index_temp]=(u.data()[index_data]);
          ++index_data;
          ++index_temp;
          for (int i=0;i<(sizep-size[D-1]);++i)
            {
              real_data[index_temp]=0;
              ++index_temp;
            }
        }

    }
#endif





#if 0

  cat::array<Real,3> uu(size[0],size[1],sizep);
  cat::array_iterator<Real,3> uuit(uu);
  for(uuit=uu.begin();uuit!=uu.end();++uuit)
    *uuit=real_data[uuit.pos()];
  cout << u << "\n\n" << uu << endl;
  exit(0);

#endif




  //     cout << total_size << " " << total_sizep << endl;
  //
  //     cout << endl << u << endl;
  //
  //     for (int i=0;i<total_sizep;++i)
  //       cout << real_data[i] << endl;


  //perform the transform
  fftw_execute(direct_plan);

  //Normalize output
  for (int i=0;i<total_sizep;++i)
    real_data[i] /= total_size;


  //write complex data in output field via pointer
  for(int i=0;i<total_sizep/2;++i)
    (u_hat.data()[i])=complex<Real>
                      (complex_data[i][0],complex_data[i][1]);
                      
                    //  cout << total_sizep << total_size << endl;
                      
                      
}

//inverse transform (scalar)
template <int D>
void s_rfft<D>::inverse_transform(RT& u,const CT& u_hat)
{

  //read complex data from input field via pointer
  for(int i=0;i<total_sizep/2;i++)
    {
      complex_data[i][0]=(u_hat.data()[i]).real();
      complex_data[i][1]=(u_hat.data()[i]).imag();
    }

  //perform the transform
  fftw_execute(inverse_plan);


  //write real data to output field via pointer
  
  #if 1
  int depad=0;
  int jb=0,je=0;
  for(int i=0;i<total_size_Dm1;i++)
    {
      je=(i+1)*sizep;
      for(int j=jb;j<je-pad_size;j++)
        (u.data()[j-depad])=real_data[j];
      jb=je;
      depad+=pad_size;
    }
#endif
   

//to adapt 
    #if 0
  int index_data=0;
  int index_temp=0;
  while(index_temp<total_sizep)
    {

      //cout << index_data << " " << index_temp << endl;

      real_data[index_temp]=(u.data()[index_data]);

      ++index_data;
      ++index_temp;

      //if ((index_data%total_size_Dm1)==0)

      if ((index_data%size[D-1])==size[D-1]-1)
        {
          real_data[index_temp]=(u.data()[index_data]);
          ++index_data;
          ++index_temp;
          for (int i=0;i<pad_size;++i)
            {
              real_data[index_temp]=0;
              ++index_temp;
            }
        }

    }
#endif

    
    
}





// template <int D>
// CT s_rfft<D>::direct_transform(const RT & u)
// {
//   //  cat::tvector<int,D> size_hat(size
//   //CT u_hat(cat::tvector<int,D>(
// }




// template <int D>
// RT s_rfft<D>::inverse_transform(const CT & u_hat)
// {

// }




#endif
