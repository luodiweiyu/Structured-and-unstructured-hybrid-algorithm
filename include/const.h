#pragma once
#include<cmath>
#include<string>
#include<vector>
using std::vector;
using std::string;
namespace ConstPara//常数
{
	const double CFL = 0.3;
	const double t_end = 300;
	const double pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628;
	const double gama = 1.4;
	const double random = 0.495;
}
namespace Init//初始网格
{
	extern double rho0;
	extern double v0;
	extern double p0;
	extern double u0;
}
namespace Normal//初始网格
{
	extern double rho1;
	extern double v1;
	extern double p1;
	extern double u1;
	extern double rho2;
	extern double v2;
	extern double p2;
	extern double u2;
}
namespace Couette//初始网格
{
	extern double rho1;
	extern double v1;
	extern double p1;
	extern double u1;
}

namespace MeshPara
{
	const string FlowType = "cylinder";
	//目前有，膨胀波Prandtl-Meyer;正激波normal，斜激波oblique ，激波相交intersection,库塔流couette，圆柱绕流cylinder， 其他other
	const int meshType = 20;
	//10正六边形网格，11扰动正六边形网格，12砖块网格，
	//20长方形结构网格
	const string methodType = "C";
	//计算方法，激波捕捉C-Capturing，激波装配F-Fitting
	const int method[12][4] = {
{1,2,0,0},//方法3
{0,1,2,2},
{2,0,1,1},

{0,1,1,2},//方法2
{1,2,2,0},
{2,0,0,1},


{0,0,1,2},//方法1
{1,1,2,0},
{2,2,0,1},

{0,1,2,0},//方法4
{2,0,1,2},
{1,2,0,1},
	};

	const int Xnum = 101;//x方向点的个数20,40,150
	const int Ynum = 101;//y方向点的个数43,87,215
	const int Pnum = Xnum * Ynum;//一共点的个数
	const double dx = 0.0002;
	const double dy = dx;
}
struct Flux
{
	double f1;
	double f2;
	double f3;
	double f4;
};

struct mesh
{
	double x;
	double y;
	int id;
	double rho;//密度
	double u;//水平速度
	double v;//竖直速度
	double p;//压强
	double x0;//动网格的初始网格
	double y0;
	double um;//格点运动速度
	double vm;
	int section = 1;//分区
	int sec_num;//分区后计算的参数1或0
	int step = 0;//表明该点在何步骤更新的
	int neiborsec = -1;//相邻分区，如果该值为负，则表示该点不是激波点，反之为相邻分区共用的激波点
	int neiborsec_ad = -1;//相邻分区的地址
	vector <double> xix;
	vector <double> xiy;
	vector <double> xit;
	vector <double> etax;
	vector <double> etay;
	vector <double> etat;
	vector <double> xxi;
	vector <double> xeta;
	vector <double> xtau;
	vector <double> yxi;
	vector <double> yeta;
	vector <double> ytau;
	vector <double> J;//雅可比行列式
	vector <mesh*>neibor;//记录该格点的相邻格点位置信息
	vector <int>moveConnct;//运动关联点，即本点运动与该moveConnect点有关
	string type;//格点类型，分为上下左右边界以及内部，U，D,L,R,IN,激波点SHOCK，接触间断点DISCON（contact discontinuity）,激波相交点CENTER
	string change = "Y";//周围网格物理量是否改变，没有“N”，改变“Y”；用于计算加速
	Flux alpha = { 0 };
	Flux beta = { 0 };


};
struct Coor//表示坐标
{
	double x;
	double y;
};
struct Line//表示直线Ax+Bx+C=0
{
	double A;
	double B;
	double C;
};
struct polygon_mesh
{
	//都只记录点的序号，用于输出多边形网格
	vector <int> node;
	vector <int>face_start;
	vector <int>face_end;
	int section;
};
namespace Oblique//斜激波
{
	extern double beta;
	extern double delta;
	extern int startpoint;
	extern double  rho1, p1, Ma1, u1, v1;
	extern double  rho2, p2, Ma2, u2, v2;
}
namespace Prandtl_Meyer//普朗特麦耶尔流动
{
	extern double  rho0, p0, theta0, lambda0, delta0, mu0, Ma0, u0, v0;
	extern double  rho1, p1, theta1, lambda1, delta1, mu1, Ma1, u1, v1;
	extern double  rho2, p2, theta2, lambda2, delta2, mu2, Ma2, u2, v2;
}
namespace ShockwaveCross//激波相交
{

	extern double  rho1, p1, beta1, delta1, Ma1, u1, v1;
	extern double  rho2, p2, beta2, delta2, Ma2, u2, v2;
	extern double  rho3, p3, beta3, delta3, Ma3, u3, v3;
	extern double  rho4, p4, beta4, delta4, Ma4, u4, v4;
	extern double  rho5, p5, beta5, delta5, Ma5, u5, v5;
}

