//
// C++ Interface: fft
//
// Description: 
//
//
// Author: Manuel Baptista <mbaptist@fc.up.pt>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "fft_types.h"

template <class TR,calss TF>
class FFT
{
private:
	//Types	
	typedef typename FFT_TYPES<TR,TF>::Real RT;
	typedef typename FFT_TYPES<TR,TF>::Complex FT;
	//Members
	RT * realdata;
	FT * fourierdata;
	typename FFT_TYPES<TR,TF>::Direct_plan direct_plan;
	typename FFT_TYPES<TR,TF>::Inverse_plan inverse_plan;

public:
	//Ctors
	FFT(typename FFT_TYPES<TR,TF>::ST size);//Ctor from size
	FFT(TR & realfield,TF	& fourierfield);//Ctor from fields
	//Dtor
	~FFT();

private:
	//Forbidden Ctors
	FFT();//Default ctor
	FFT(const & FFT);//Copy ctor
	
private:
	//Private methods 	
	void set_data(TR & realfield,TF	& fourierfield);
	void create_plans();
	void destroy_plans();

public:
	//Public methods
	//void switchmode();
	
	
};
