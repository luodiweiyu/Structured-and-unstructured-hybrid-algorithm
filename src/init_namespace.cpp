#include"include/const.h"
#include"include/shockwave.h"
#include"include/Prandtl-Meyer.h"
#include"include/patition.h"
using namespace MeshPara;
using ConstPara::t_end;


namespace Init//≥ı ºÕ¯∏Ò
{
	double rho0 = 1;
	double v0 = 0;
	double p0 = 1;
	double u0 = 0 * sqrt(ConstPara::gama * p0 / rho0);
}
namespace Normal
{
	//double rho1 = 1;
	//double v1 = 0;
	//double p1 = 1;
	//double u1 = 1.1 * sqrt(ConstPara::gama * p1 / rho1);
	//double rho2;
	//double v2;
	//double p2;
	//double u2;
	double rho1 = 5;
	double v1 = 0;
	double p1 = 4;
	double u1 = 3 * sqrt(ConstPara::gama * p1 / rho1);
	double rho2 = 3.0919952786742946;
	double v2 = 0;
	double p2 = 5.5743027276648576;
	double u2 = 2.2088635344975316;

}
namespace Couette
{
	double rho1 = 5;
	double v1 = 0;
	double p1 = 4;
	double u1 = 0.08 * sqrt(ConstPara::gama * p1 / rho1);


}

namespace Oblique//–±º§≤®
{
	double beta = 50* ConstPara::pi / 180;
	double rho1 = 2;
	double v1 = 0;
	double p1 = 3;
	double u1 = 1.01* sqrt(ConstPara::gama * p1 / rho1);
	int startpoint = MeshPara::Xnum * 3 / 10;
	double Ma1 = get_Ma(u1, v1, rho1, p1);
	double Ma2 = get_Ma2(Ma1, beta);
	double p2 = get_p2(Ma1, beta, p1);
	double rho2 = get_rho2(Ma1, beta, rho1);
	double delta = get_delta(Ma1, beta);
	double u2 = get_ufromMa2(Ma2, rho2, p2, delta);
	double v2 = get_vfromMa2(Ma2, rho2, p2, delta);
}
namespace Prandtl_Meyer//∆’¿ Ãÿ¬Û“Æ∂˚¡˜∂Ø
{
	double rho1 = 2;
	double v1 = 0;
	double p1 = 3;
	double u1 = 1.3 * sqrt(ConstPara::gama * p1 / rho1);
	double Ma1 = get_Ma(u1, v1, rho1, p1);
	double lambda1 = getlambdafromMa(Ma1);
	double delta1 = getdeltafromlambda(lambda1);
	double mu1 = getmufromMa(Ma1);
	double theta1 = getthetafromMa(Ma1);
	double delta2 = 10.0 / 180 * ConstPara::pi;
	double lambda2 = getlambdafromdelta(delta2);
	double Ma2 = getMafromlambda(lambda2);
	double mu2 = getmufromMa(Ma2);
	double theta2 = getthetafromMa(Ma2);
	double p0 = getp0fromlambdaandp(lambda1, p1);
	double rho0 = getrho0fromlambdaandrho(lambda1, rho1);
	double p2 = getpfromlambda(lambda2, p0);
	double rho2 = getrhofromlambda(lambda2, rho0);
	double u2 = Ma2 * sqrt(p2 * ConstPara::gama / rho2) * cos(delta2);
	double v2 = -Ma2 * sqrt(p2 * ConstPara::gama / rho2) * sin(delta2);

}
namespace ShockwaveCross//º§≤®œ‡Ωª
{
	double rho1 = 8;
	double v1 = 0;
	double p1 = 5;
	double u1 = 10* sqrt(ConstPara::gama * p1 / rho1);
	double beta2 = -30/ 180.0 * ConstPara::pi;
	double beta3 = 30/ 180.0 * ConstPara::pi;
	double Ma1 = get_Ma(u1, v1, rho1, p1);

	double Ma2 = get_Ma2(Ma1, beta2);
	double p2 = get_p2(Ma1, beta2, p1);
	double rho2 = get_rho2(Ma1, beta2, rho1);
	double delta2 = get_delta(Ma1, beta2);
	double u2 = get_ufromMa2(Ma2, rho2, p2, delta2);
	double v2 = get_vfromMa2(Ma2, rho2, p2, delta2);

	double Ma3 = get_Ma2(Ma1, beta3);
	double p3 = get_p2(Ma1, beta3, p1);
	double rho3 = get_rho2(Ma1, beta3, rho1);
	double delta3 = get_delta(Ma1, beta3);
	double u3 = get_ufromMa2(Ma3, rho3, p3, delta3);
	double v3 = get_vfromMa2(Ma3, rho3, p3, delta3);

	double p4 = 1;
	double rho4 = 1;
	double delta4 = 1;
	double u4 = 1;
	double v4 = 1;
	double beta4 = 1;

	double p5 = 1;
	double rho5 = 1;
	double delta5 = 1;
	double u5 = 1;
	double v5 = 1;
	double beta5 = 1;
}
