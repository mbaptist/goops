


// namespace goops
// {

// #include <fftw3.h>
// 
// #include <cmath>

using namespace std;

//Generic Ctor
template <class RealType,class FourierType>
FFT<RealType,FourierType>::FFT()
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
	direct_plan.normalise();
}

//Generic inverse transform
template <class RealType,class FourierType>
void FFT<RealType,FourierType>::inverse_transform(RealType & realfield,const FourierType & fourierfield)
{
	inverse_plan.switch_data(realfield,fourierfield);
	inverse_plan.execute();
	inverse_plan.normalise();
}
