// -*- C++ -*-


////////////////////////////
// BASE CLASS TO FFT   
///////////////////////////


#ifndef FFT_BASE_H
#define FFT_BASE_H



template <int D>
class fft_base
{
protected:
  //sizes
  int * size;//Size of the vector to be transformed, along each dimension
  int total_size;//Total size (with last dimension not padded)    
  //plan variables
  fftw_plan direct_plan;//direct transforms
  fftw_plan inverse_plan;//inverse transforms
protected:
  //Constructors
  fft_base(cat::tvector<int,D> size_);
  fft_base(cat::tvector<int,D> size_,int num_threads);
  //Destructor
  ~fft_base();
public:
  void resize(cat::tvector<int,D> size_);
private:
  void initialise(cat::tvector<int,D> size_);
};

template <int D>
class real_fft_base : public fft_base<D>
{
protected:
  //sizes
  using fft_base<D>::size;//Size of the vector to be transformed, along each dimension
  using fft_base<D>::total_size;//Total size (with last dimension not padded) 
  int sizep;//Size along the padded dimension
  int pad_size;//Size of padding
  int total_size_Dm1;//Total size without last dimension
  int total_sizep;//Total size padded (with last dimension padded)
  //plan variables
  using fft_base<D>::direct_plan;//direct transforms
  using fft_base<D>::inverse_plan;//inverse transforms
protected:
  //Constructors
  real_fft_base(cat::tvector<int,D> size_);
  real_fft_base(cat::tvector<int,D> size_,int num_threads);
  //Destructor
  ~real_fft_base();
private:
  void initialise(cat::tvector<int,D> size_);
};

//#endif



//IMPLEMENTATION

//class fft_base

//contructor from size
template <int D>
fft_base<D>::fft_base(cat::tvector<int,D> size_)
{
//	fftw_init_threads();
  //fftw_plan_with_nthreads(2);


  //initialise sizes
  initialise(size_);
}
//constructor from size and number of threads
template <int D>
fft_base<D>::fft_base(cat::tvector<int,D> size_,int num_threads)
{
  //initialise fftw for threads
  fftw_init_threads();
  fftw_plan_with_nthreads(num_threads);
  //initialise sizes
  initialise(size_);
}

//destructor
template <int D>
fft_base<D>::~fft_base()
{
}

template <int D>
void fft_base<D>::initialise(cat::tvector<int,D> size_)
{
  //initialise size
  size=new int[D];
  for (int i=0;i<D;i++)
    size[i]=size_[i];
  //evaluate total size
  total_size=1;
  for (int i=0;i<D;i++)
    total_size *= size[i];
}



//class real_fft_base

//contructor from size
template <int D>
real_fft_base<D>::real_fft_base(cat::tvector<int,D> size_):
  fft_base<D>(size_)
{
  initialise(size_);
}
//constructor from size and number of threads
template <int D>
real_fft_base<D>::real_fft_base(cat::tvector<int,D> size_,int num_threads):
  fft_base<D>(size_,num_threads)
{
  initialise(size_);
}


//destructor
template <int D>
real_fft_base<D>::~real_fft_base()
{
}

template <int D>
void real_fft_base<D>::initialise(cat::tvector<int,D> size_)
{
    //initialize size
  size=new int[D];
  for (int i=0;i<D;i++)
    size[i]=size_[i];
  //evaluate pad_size
  pad_size=2-size[D-1]%2;
  //evaluate sizep
  sizep=2*(size[D-1]/2+1);
  //evaluate total size without last dimension
  total_size_Dm1=total_size/size[D-1];
  //evaluate total size padded
  total_sizep=total_size_Dm1*sizep;
}

#endif
