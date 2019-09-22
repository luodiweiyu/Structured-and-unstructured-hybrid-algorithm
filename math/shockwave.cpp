#include"const.h"
#include<cmath>
#include"shockwave.h"
//计算激波相关参数钱翼稷《空气动力学》p228-241
using namespace ConstPara;
double get_c(double ρ, double p)//声速公式
{
	return sqrt(γ * p / ρ);
}
double get_Ma(double u, double v, double ρ, double p)//求马赫数
{
	double a = get_c(ρ, p);
	double velocity = sqrt(u * u + v * v);
	return velocity / a;
}
//正激波关系式
//7-99
double get_Mas(double p1, double p2)
{
	double Mas = 1 + (γ + 1) * (p2 / p1 - 1) / (2 * γ);
	return sqrt(Mas);
}
//7-105
double get_p2p1(double Ma1)
{
	return  2 * γ* Ma1* Ma1 / (γ + 1) - (γ - 1) / (γ + 1);
}
//7-106
double get_ρ2ρ1(double Ma1)
{
	return  (γ + 1) * Ma1 * Ma1 / ((γ - 1) * Ma1 * Ma1 + 2);
}
double get_Ma2(double Ma1)
{
	return sqrt((Ma1 * Ma1 + 2 / (γ - 1)) / (2 * γ * Ma1 * Ma1 / (γ - 1) - 1));
}
//求上游马赫数  邹东阳,博士论文p57 3.5
double get_Mu(mesh & U, mesh & D, double θ)
{
	using namespace ConstPara;
	double δ = 1e-10;
	double Mu = 2, Mu1 = 20;
	double udn, uun, udt, uut;
	uun = get_un(U, θ);
	udn = get_un(D, θ);
	uut = get_ut(U, θ);
	udt = get_ut(D, θ);
	double cu = sqrt(γ * U.p / U.ρ);
	double J, Fx, fx;
	double a, b;
	while (abs(Mu - Mu1) >= δ)
	{
		Mu = Mu1;
		double J = 2 * sqrt(γ * D.p / D.ρ) / (γ - 1) - udn;
		a = 2 * γ * Mu1 * Mu1 - (γ - 1);
		b = (γ - 1) * Mu1 * Mu1 + 2;
		Fx = sqrt(a * b) / (Mu1 * (γ - 1)) + (Mu1 * Mu1 - 1) / Mu1 - (γ + 1) * (J + uun) / (2 * cu);
		fx = (4 * γ * Mu1 * b + 2 * (γ - 1) * Mu1 * a) / (2 * Mu1 * (γ - 1) * sqrt(a * b)) - sqrt(a * b) / (Mu1 * Mu1 * (γ - 1)) + 1 + 1 / (Mu1 * Mu1);
		Mu1 = Mu1 - Fx / fx;
	}
	return Mu;
}
//3.4
double get_udn(mesh & U, mesh & D, double Ma1, double Vs, double θ)
{
	double λ1 = getλfromMa(Ma1);
	double λ2 = 1 / λ1;
	double Ma2 = getMafromλ(λ2);
	double c2 = sqrt(γ * D.p / D.ρ);
	double V2 = Ma2 * c2;
	return V2 + Vs;
	//double a = ((γ - 1) * Ma1 + 2) / (2 * γ * Ma1 - (γ - 1));
	//double b = (a * γ * D.p) / D.ρ;
	//double cu = sqrt(γ * U.p / U.ρ);
	//double uun = U.u * sin(θ) - U.v * cos(θ);
	//double udn = sqrt(b) - uun +cu * Ma1 ;
	//return udn;
}

double getλfromMa(double Ma)
{
	double a = 2 + (γ - 1) * Ma * Ma;
	double b = 0;
	double c = -(γ + 1) * Ma * Ma;
	return (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
}
double getMafromλ(double λ)
{
	double λ2 = λ * λ;
	return sqrt((2 * λ2 / (γ + 1)) / (1 - (γ - 1) / (γ + 1) * λ2));
}

//以下为斜激波,经验证只适合波前速度水平，v=0,不具有通用性
double get_Ma2(double Ma1, double β)//斜激波波后马赫数，钱翼稷《空气动力学》p241,7-131
{
	double Ma2 = (Ma1 * Ma1 + 2 / (γ - 1)) / (2 * γ / (γ - 1) * Ma1 * Ma1 * sin(β) * sin(β) - 1) + 2 / (γ - 1) * Ma1 * Ma1 * cos(β) * cos(β) / (Ma1 * Ma1 * sin(β) * sin(β) + 2 / (γ - 1));
	Ma2 = sqrt(Ma2);
	return Ma2;
}
//double get_δ(double Ma1, double β)//气流折射角
//{
//
//	double tanδ = (Ma1 * Ma1 * sin(β) * sin(β) - 1) / ((Ma1 * Ma1 * ((γ + 1) / 2 - sin(β) * sin(β)) + 1) * tan(β));
//	double δ = atan(tanδ);
//	//if (δ < 0)
//	//	δ = δ + pi;
//	return δ;
//}
double get_βfromδ(double Ma1, double δ)
{
	//弦截法 计算机科学与技术p132
	double β = 40 * pi / 180, β1 = 50 * pi / 180, β0 = 60 * pi / 180;
	double fβ1, fβ0;
	if (δ == 0.92846675631531839)
		δ = 0.92846675631531839;
	while (abs(β - β1) > 1e-10)
	{
		β0 = β1;
		β1 = β;
		while (β1 > pi / 2 || β1 < -pi / 2)
		{
			if (β1 > pi / 2)
				β1 = β1 - pi / 2;
			else if (β1 < -pi / 2)
				β1 = β1 + pi / 2;
		}
		while (β0 > pi / 2 || β0 < -pi / 2)
		{
			if (β0 > pi / 2)
				β0 = β0 - pi / 2;
			else if (β0 < -pi / 2)
				β0 = β0 + pi / 2;
		}

		fβ1 = tan(δ) - (Ma1 * Ma1 * sin(β1) * sin(β1) - 1) / ((Ma1 * Ma1 * ((γ + 1) / 2 - sin(β1) * sin(β1)) + 1) * tan(β1));
		fβ0 = tan(δ) - (Ma1 * Ma1 * sin(β0) * sin(β0) - 1) / ((Ma1 * Ma1 * ((γ + 1) / 2 - sin(β0) * sin(β0)) + 1) * tan(β0));
		β = β1 - fβ1 / (fβ1 - fβ0) * (β1 - β0);
	}
	return β;
}
double get_ufromMa2(double Ma2, double ρ2, double p2, double δ)
{
	double a2 = sqrt(γ * p2 / ρ2);
	double velocity = Ma2 * a2;
	double u = velocity * cos(δ);
	return u;
}
double get_vfromMa2(double Ma2, double ρ2, double p2, double δ)
{
	double a2 = sqrt(γ * p2 / ρ2);
	double velocity = Ma2 * a2;
	double v = velocity * sin(δ);
	return v;
}
double get_p2(double Ma1, double β, double p1)
{
	double p2 = (2 * γ / (γ + 1) * Ma1 * Ma1 * sin(β) * sin(β) - (γ - 1) / (γ + 1)) * p1;
	return p2;
}
double get_ρ2(double Ma1, double β, double ρ1)
{
	double ρ2 = ((γ + 1) * Ma1 * Ma1 * sin(β) * sin(β) / ((γ - 1) * Ma1 * Ma1 * sin(β) * sin(β) + 2)) * ρ1;
	return ρ2;
}


//以下为通用型斜激波算法
double get_un(mesh & A, double β)
{
	if (β >= 0)
		return A.u * sin(β) - A.v * cos(β);
	else
		return A.u * sin(-β) + A.v * cos(-β);
}
double get_ut(mesh & A, double β)
{
	if (β >= 0)
		return A.u * cos(β) + A.v * sin(β);
	else
		return A.u * cos(-β) - A.v * sin(-β);
}

void get_down(mesh & U, mesh & D, double β)//下游参数
{
	double udn, uun, udt, uut;
	uun = get_un(U, β);
	uut = get_ut(U, β);
	double au = get_c(U.ρ, U.p);
	double Mu = uun / au;
	D.p = U.p * get_p2p1(Mu);
	D.ρ = U.ρ * get_ρ2ρ1(Mu);
	double Md = get_Ma2(Mu);
	double ad = get_c(D.ρ, D.p);
	udn = Md * ad;
	udt = uut;
	if (β >= 0)
	{
		D.u = udt * cos(β) + udn * sin(β);
		D.v = udt * sin(β) - udn * cos(β);
	}
	else
	{
		D.u = udt * cos(-β) + udn * sin(-β);
		D.v = -udt * sin(-β) + udn * cos(-β);
	}
}
double get_δ(double u, double v)
{
	double δ = atan(abs(v) / abs(u));
	if (u > 0)
	{
		if (v > 0)
			return δ;
		else if (v < 0)
			return -δ;
		else
			return 0;
	}
	else if (u < 0)
	{
		if (v > 0)
			return -δ;
		else if (v < 0)
			return δ;
		else
			return 0;
	}
	else
	{
		if (v > 0)
			return pi / 2;
		else if (v < 0)
			return -pi / 2;
		else
			return 0;
	}
}
double get_β(mesh & U, double p2, int type)//type=-1表示求出的β角为负，否则为正
{
	double β = 40 * pi / 180, β1 = 50 * pi / 180, β0 = 60 * pi / 180;
	double uun0, uun1;
	double Mu0, Mu1;
	double p20, p21;
	double fβ1, fβ0;
	while (abs(β - β1) > 1e-15)
	{
		while (β > pi / 2 || β < -pi / 2)
		{
			if (β > pi / 2)
				β = β - pi / 2;
			else if (β < -pi / 2)
				β = β + pi / 2;
		}
		β0 = β1;
		β1 = β;
		if (type == -1)
		{
			if (β0 > 0)
				β0 = -β0;
			if (β1 > 0)
				β1 = -β1;
		}
		else if (type == 1)
		{
			if (β0 < 0)
				β0 = -β0;
			if (β1 < 0)
				β1 = -β1;
		}
		uun0 = get_un(U, β0);
		uun1 = get_un(U, β1);
		Mu0 = get_Ma(uun0, 0, U.ρ, U.p);
		Mu1 = get_Ma(uun1, 0, U.ρ, U.p);
		p20 = U.p * get_p2p1(Mu0);
		p21 = U.p * get_p2p1(Mu1);
		fβ0 = p20 - p2;
 		fβ1 = p21 - p2;
		if (fβ0 == fβ1)
			break;
		else
			β = β1 - fβ1 / (fβ1 - fβ0) * (β1 - β0);
	}
	return β;

}