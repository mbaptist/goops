
//Generic Ctor
template <class RealType,class FourierType>
FFT<RealType,FourierType>::FFT():
realdata(0),
fourierdata(0)
{
}

//Generic Dtor
template <class RealType,class FourierType>
FFT<RealType,FourierType>::~FFT()
{
	realdata=0;
	fourierdata=0;
}

//Generic direct transform
template <class RealType,class FourierType>
void FFT<RealType,FourierType>::direct_transform(FourierType & fourierfield,RealType & realfield)
{
	if(realdata!=realfield.data()||fourierdata!=fourierfield.data())
	{
		destroy_plans();
		realdata=realfield.data();
		fourierdata=fourierfield.data();
		size=fourierfield.size();
		create_plans();
	}
	fftw_execute(direct_plan);
	fourierfield/=size;
}

//Generic inverse transform
template <class RealType,class FourierType>
void FFT<RealType,FourierType>::(RealType & realfield,FourierType & fourierfield)
{
	if(realdata!=realfield.data()||fourierdata!=fourierfield.data())
	{
		destroy_plans();
		realdata=realfield.data();
		fourierdata=fourierfield.data();
		size=fourierfield.size();
		create_plans();
	}
	fftw_execute(inverse_plan);
	realfield/=size;
}

