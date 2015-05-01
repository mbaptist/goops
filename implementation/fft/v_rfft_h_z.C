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
    /*	u_hat[0]=u_hat_x;
	u_hat[1]=u_hat_y;
	u_hat[2]=u_hat_z;*/
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
    // 	u[0]=u_x;
    // 	u[1]=u_y;
    // 	u[2]=u_z;
}

//}
