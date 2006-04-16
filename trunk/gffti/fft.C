


// namespace goops
// {

// #include <fftw3.h>
// 
// #include <cmath>

using namespace std;

//Generic Ctors
template <class RealType,class FourierType>
FFT_BASE<RealType,FourierType>::FFT():
direct_plan(),
inverse_plan()
{
}
template <class RealType,class FourierType>
FFT_BASE<RealType,FourierType>::FFT(const string & direction):
direct_plan(direction),
inverse_plan(direction)
{
}
template <class RealType,class FourierType>
FFT_BASE<RealType,FourierType>::FFT(const string & subtype,const string & direction):
direct_plan(subtype,direction),
inverse_plan(subtype,direction)
{
}

//Generic Dtor
template <class RealType,class FourierType>
FFT<RealType,FourierType>::~FFT()
{
}

//Generic direct transform
template <class RealType,class FourierType>
void FFT<RealType,FourierType>::direct_transform(FourierType & fourierfield, const RealType & realfield)
{
	direct_plan.switch_data(realfield,fourierfield);
	direct_plan.execute();
	fourierfield/=realfield.size();
}

//Generic inverse transform
template <class RealType,class FourierType>
void FFT<RealType,FourierType>::inverse_transform(RealType & realfield,const FourierType & fourierfield)
{
	inverse_plan.switch_data(realfield,fourierfield);
	inverse_plan.execute();
}


///////////////////////////////

//Complex to Complex Ctor
template <>
template<int D>
FFT<cat::array<CS,D>,cat::array<CS,D> >::FFT(const string & direction):
FFT_BASE(direction)
{
}

//1D Real to Real Ctor
template <>
FFT<cat::array<RS,1>,cat::array<RS,1> >::FFT(const string & subtype,const string & direction):
FFT_BASE(const string & subtype,const string & direction)
{
}

//1D Real to Real direct transform
template <>
void FFT<cat::array<RS,1>,cat::array<RS,1> >::direct_transform(cat::array<RS,1> & fourierfield, const cat::array<RS,1> & realfield)
{
	direct_plan.switch_data(realfield,fourierfield);
	direct_plan.execute();
	direct_plan.normalise();
}

//1D Real to Real inverse transform
template <>
void FFT<cat::array<RS,1>,cat::array<RS,1> >::inverse_transform(cat::array<RS,1> & realfield,const cat::array<RS,1> & fourierfield)
{
	inverse_plan.switch_data(realfield,fourierfield);
	inverse_plan.execute();
	inverse_plan.normalise();
}

