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
	RT u_x(u[0]);
	CT u_hat_x(u_hat[0]);
	s_rfft_x_obj.direct_transform(u_hat_x,u_x);
	RT u_y(u[1]);
	CT u_hat_y(u_hat[1]);
	s_rfft_y_obj.direct_transform(u_hat_y,u_y);
	RT u_z(u[2]);
	CT u_hat_z(u_hat[2]);
	s_rfft_z_obj.direct_transform(u_hat_z,u_z);
	
//   //First two components use cos along z
// 	for(int comp=0;comp<2;++comp)
// 	{				
// 		
// 		RT s_u(u.shape());
// 		CT s_u_hat(u_hat.shape());
// 		
// 		for(int i=0;i<s_u.size();++i)
// 			s_u.data()[i]=(u.data()[i])[comp];
// 		
// 		s_rfft_cos_obj.direct_transform(s_u_hat,s_u);
// 		
// 		for(int i=0;i<s_u_hat.size();++i)
// 			(u_hat.data()[i])[comp]=s_u_hat.data()[i];
// 	}
//   //Third component uses sin along z
// 	RT s_u(u.shape());
// 	CT s_u_hat(u_hat.shape());
// 	for(int i=0;i<s_u.size();++i)
// 		s_u.data()[i]=(u.data()[i])[2];
// 	s_rfft_sin_obj.direct_transform(s_u_hat,s_u);
// 	for(int i=0;i<s_u_hat.size();++i)
// 		(u_hat.data()[i])[2]=s_u_hat.data()[i];
}

//inverse transform (vector)
void v_rfft_h_z_ccs::inverse_transform(VRT& u,const VCT& u_hat)
{
	RT u_x(u[0]);
	CT u_hat_x(u_hat[0]);
	s_rfft_x_obj.inverse_transform(u_x,u_hat_x);
	RT u_y(u[1]);
	CT u_hat_y(u_hat[1]);
	s_rfft_y_obj.inverse_transform(u_y,u_hat_y);
	RT u_z(u[2]);
	CT u_hat_z(u_hat[2]);
	s_rfft_z_obj.inverse_transform(u_z,u_hat_z);
	
//   //First two components use cos along z
// 	for(int comp=0;comp<2;++comp)
// 	{
// 		RT s_u(u.shape());
// 		CT s_u_hat(u_hat.shape());
// 		for(int i=0;i<s_u_hat.size();++i)
// 			s_u_hat.data()[i]=(u_hat.data()[i])[comp];		
// 		s_rfft_cos_obj.inverse_transform(s_u,s_u_hat);
// 		for(int i=0;i<s_u.size();++i)
// 			(u.data()[i])[comp]=s_u.data()[i];
// 	}
//   //Third component uses sin along z
// 	RT s_u(u.shape());
// 	CT s_u_hat(u_hat.shape());
// 	for(int i=0;i<s_u_hat.size();++i)
// 		s_u_hat.data()[i]=(u_hat.data()[i])[2];
// 	s_rfft_sin_obj.inverse_transform(s_u,s_u_hat);
// 	for(int i=0;i<s_u.size();++i)
// 		(u.data()[i])[2]=s_u.data()[i];
}

//}
