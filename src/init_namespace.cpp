#include"const.h"
#include"shockwave.h"
#include"Prandtl-Meyer.h"
#include"patition.h"
using namespace MeshPara;
using ConstPara::t_end;


namespace Init//場宎厙跡
{
	double 老0 = 1;
	double v0 = 0;
	double p0 = 1;
	double u0 = 0 * sqrt(ConstPara::gama * p0 / 老0);
}
namespace Normal
{
	//double 老1 = 1;
	//double v1 = 0;
	//double p1 = 1;
	//double u1 = 1.1 * sqrt(ConstPara::gama * p1 / 老1);
	//double 老2;
	//double v2;
	//double p2;
	//double u2;
	double 老1 = 5;
	double v1 = 0;
	double p1 = 4;
	double u1 = 3 * sqrt(ConstPara::gama * p1 / 老1);
	double 老2 = 3.0919952786742946;
	double v2 = 0;
	double p2 = 5.5743027276648576;
	double u2 = 2.2088635344975316;

}
namespace Couette
{
	double 老1 = 5;
	double v1 = 0;
	double p1 = 4;
	double u1 = 0.08 * sqrt(ConstPara::gama * p1 / 老1);


}

namespace Oblique//訇慾疏
{
	double beta = 50* ConstPara::pi / 180;
	double 老1 = 2;
	double v1 = 0;
	double p1 = 3;
	double u1 = 1.01* sqrt(ConstPara::gama * p1 / 老1);
	int startpoint = MeshPara::Xnum * 3 / 10;
	double Ma1 = get_Ma(u1, v1, 老1, p1);
	double Ma2 = get_Ma2(Ma1, beta);
	double p2 = get_p2(Ma1, beta, p1);
	double 老2 = get_老2(Ma1, beta, 老1);
	double 汛 = get_汛(Ma1, beta);
	double u2 = get_ufromMa2(Ma2, 老2, p2, 汛);
	double v2 = get_vfromMa2(Ma2, 老2, p2, 汛);
}
namespace Prandtl_Meyer//ぱ檄杻闔珖嫌霜雄
{
	double 老1 = 2;
	double v1 = 0;
	double p1 = 3;
	double u1 = 1.3 * sqrt(ConstPara::gama * p1 / 老1);
	double Ma1 = get_Ma(u1, v1, 老1, p1);
	double 竹1 = get竹fromMa(Ma1);
	double 汛1 = get汛from竹(竹1);
	double 米1 = get米fromMa(Ma1);
	double 牟1 = get牟fromMa(Ma1);
	double 汛2 = 10.0 / 180 * ConstPara::pi;
	double 竹2 = get竹from汛(汛2);
	double Ma2 = getMafrom竹(竹2);
	double 米2 = get米fromMa(Ma2);
	double 牟2 = get牟fromMa(Ma2);
	double p0 = getp0from竹andp(竹1, p1);
	double 老0 = get老0from竹and老(竹1, 老1);
	double p2 = getpfrom竹(竹2, p0);
	double 老2 = get老from竹(竹2, 老0);
	double u2 = Ma2 * sqrt(p2 * ConstPara::gama / 老2) * cos(汛2);
	double v2 = -Ma2 * sqrt(p2 * ConstPara::gama / 老2) * sin(汛2);

}
namespace ShockwaveCross//慾疏眈蝠
{
	double 老1 = 8;
	double v1 = 0;
	double p1 = 5;
	double u1 = 10* sqrt(ConstPara::gama * p1 / 老1);
	double beta2 = -30/ 180.0 * ConstPara::pi;
	double beta3 = 30/ 180.0 * ConstPara::pi;
	double Ma1 = get_Ma(u1, v1, 老1, p1);

	double Ma2 = get_Ma2(Ma1, beta2);
	double p2 = get_p2(Ma1, beta2, p1);
	double 老2 = get_老2(Ma1, beta2, 老1);
	double 汛2 = get_汛(Ma1, beta2);
	double u2 = get_ufromMa2(Ma2, 老2, p2, 汛2);
	double v2 = get_vfromMa2(Ma2, 老2, p2, 汛2);

	double Ma3 = get_Ma2(Ma1, beta3);
	double p3 = get_p2(Ma1, beta3, p1);
	double 老3 = get_老2(Ma1, beta3, 老1);
	double 汛3 = get_汛(Ma1, beta3);
	double u3 = get_ufromMa2(Ma3, 老3, p3, 汛3);
	double v3 = get_vfromMa2(Ma3, 老3, p3, 汛3);

	double p4 = 1;
	double 老4 = 1;
	double 汛4 = 1;
	double u4 = 1;
	double v4 = 1;
	double beta4 = 1;

	double p5 = 1;
	double 老5 = 1;
	double 汛5 = 1;
	double u5 = 1;
	double v5 = 1;
	double beta5 = 1;
}
