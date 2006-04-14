/////////////////////////////////////////////////
// Real TO COMPLEX TRANSFORM FOR VECTOR FIELDS //
/////////////////////////////////////////////////

// namespace goops
// {

#include "v_rfft_h_z_ccs.h"

#include <cat.h>
using namespace cat;

//IMPLEMENTATION
//direct transform (vector)
void v_rfft_h_z_ccs::direct_transform(VCT& u_hat,const VRT& u)
{
  //First two components use cos along z
	for(int comp=0;comp<2;++comp)
	{				
		
		RT s_u(u.shape());
		CT s_u_hat(u_hat.shape());
		
		for(int i=0;i<s_u.size();++i)
			s_u.data()[i]=(u.data()[i])[comp];
		
		s_rfft_cos_obj.direct_transform(s_u_hat,s_u);
		
		for(int i=0;i<s_u_hat.size();++i)
			(u_hat.data()[i])[comp]=s_u_hat.data()[i];
	}
  //Third component uses sin along z
	RT s_u(u.shape());
	CT s_u_hat(u_hat.shape());
	for(int i=0;i<s_u.size();++i)
		s_u.data()[i]=(u.data()[i])[2];
	s_rfft_sin_obj.direct_transform(s_u_hat,s_u);
	for(int i=0;i<s_u_hat.size();++i)
		(u_hat.data()[i])[2]=s_u_hat.data()[i];
}

//inverse transform (vector)
void v_rfft_h_z_ccs::inverse_transform(VRT& u,const VCT& u_hat)
{
  //First two components use cos along z
	for(int comp=0;comp<2;++comp)
	{
		RT s_u(u.shape());
		CT s_u_hat(u_hat.shape());
		for(int i=0;i<s_u_hat.size();++i)
			s_u_hat.data()[i]=(u_hat.data()[i])[comp];		
		s_rfft_cos_obj.inverse_transform(s_u,s_u_hat);
		for(int i=0;i<s_u.size();++i)
			(u.data()[i])[comp]=s_u.data()[i];
	}
  //Third component uses sin along z
	RT s_u(u.shape());
	CT s_u_hat(u_hat.shape());
	for(int i=0;i<s_u_hat.size();++i)
		s_u_hat.data()[i]=(u_hat.data()[i])[2];
	s_rfft_sin_obj.inverse_transform(s_u,s_u_hat);
	for(int i=0;i<s_u.size();++i)
		(u.data()[i])[2]=s_u.data()[i];
}

//}
