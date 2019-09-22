#include"/Structured-and-unstructured-hybrid-algorithm/include/const.h"
#include"/Structured-and-unstructured-hybrid-algorithm/include/shockwave.h"
#include"/Structured-and-unstructured-hybrid-algorithm/include/Prandtl-Meyer.h"
#include"/Structured-and-unstructured-hybrid-algorithm/include/patition.h"
using namespace MeshPara;
using ConstPara::t_end;


namespace Init//��ʼ����
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

namespace Oblique//б����
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
	double �� = get_��(Ma1, beta);
	double u2 = get_ufromMa2(Ma2, rho2, p2, ��);
	double v2 = get_vfromMa2(Ma2, rho2, p2, ��);
}
namespace Prandtl_Meyer//��������Ү������
{
	double rho1 = 2;
	double v1 = 0;
	double p1 = 3;
	double u1 = 1.3 * sqrt(ConstPara::gama * p1 / rho1);
	double Ma1 = get_Ma(u1, v1, rho1, p1);
	double ��1 = get��fromMa(Ma1);
	double ��1 = get��from��(��1);
	double ��1 = get��fromMa(Ma1);
	double theta1 = getthetafromMa(Ma1);
	double ��2 = 10.0 / 180 * ConstPara::pi;
	double ��2 = get��from��(��2);
	double Ma2 = getMafrom��(��2);
	double ��2 = get��fromMa(Ma2);
	double theta2 = getthetafromMa(Ma2);
	double p0 = getp0from��andp(��1, p1);
	double rho0 = getrho0from��andrho(��1, rho1);
	double p2 = getpfrom��(��2, p0);
	double rho2 = getrhofrom��(��2, rho0);
	double u2 = Ma2 * sqrt(p2 * ConstPara::gama / rho2) * cos(��2);
	double v2 = -Ma2 * sqrt(p2 * ConstPara::gama / rho2) * sin(��2);

}
namespace ShockwaveCross//�����ཻ
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
	double ��2 = get_��(Ma1, beta2);
	double u2 = get_ufromMa2(Ma2, rho2, p2, ��2);
	double v2 = get_vfromMa2(Ma2, rho2, p2, ��2);

	double Ma3 = get_Ma2(Ma1, beta3);
	double p3 = get_p2(Ma1, beta3, p1);
	double rho3 = get_rho2(Ma1, beta3, rho1);
	double ��3 = get_��(Ma1, beta3);
	double u3 = get_ufromMa2(Ma3, rho3, p3, ��3);
	double v3 = get_vfromMa2(Ma3, rho3, p3, ��3);

	double p4 = 1;
	double rho4 = 1;
	double ��4 = 1;
	double u4 = 1;
	double v4 = 1;
	double beta4 = 1;

	double p5 = 1;
	double rho5 = 1;
	double ��5 = 1;
	double u5 = 1;
	double v5 = 1;
	double beta5 = 1;
}
