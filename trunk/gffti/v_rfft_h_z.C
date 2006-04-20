/////////////////////////////////////////////////
// Real TO COMPLEX TRANSFORM FOR VECTOR FIELDS //
/////////////////////////////////////////////////

// namespace goops
// {

#include "v_rfft_h_z.h"

#include <cat.h>
using namespace cat;

//IMPLEMENTATION
//direct transform (vector)
void v_rfft_h_z::direct_transform(VCT& u_hat,const VRT& u)
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
	}

//inverse transform (vector)
void v_rfft_h_z::inverse_transform(VRT& u,const VCT& u_hat)
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
}

//}
