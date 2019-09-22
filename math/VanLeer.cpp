#include<cmath>
#include"const.h"
#include<iostream>
#include<omp.h>
#include"functions.h"
using namespace std;
using namespace ConstPara;
Flux VanLeerA(mesh& C, double ¦Îx, double ¦Îy, double ¦Ît,double J)
{
	double u = C.u;
	double v = C.v;
	double p = C.p;
	double ¦Ñ = C.¦Ñ;
	double c = sqrt(¦Ã * p / ¦Ñ);
	double D¦Î = sqrt(¦Îx * ¦Îx + ¦Îy * ¦Îy);
	double ub = u * ¦Îx + v * ¦Îy + ¦Ît;
	double uc = ub / D¦Î;
	double ¦Î1 = ¦Îx / D¦Î;
	double ¦Î2 = ¦Îy / D¦Î;
	double ¦Î3 = ¦Ît / D¦Î;
	Flux FL;
	if (uc + c <= 0)
	{
		FL = { 0 };
	}
	else if (uc - c >= 0)
	{
		double E = 0.5 * ¦Ñ * (u * u + v * v) + p / (¦Ã - 1);
		FL.f1 = ¦Ñ * uc;
		FL.f2 = ¦Ñ * uc * u + p * ¦Î1;
		FL.f3 = ¦Ñ * uc * v + p * ¦Î2;
		FL.f4 = uc * (E + p) - p * ¦Î3;
	}
	double M = uc / c;
	if (abs(M) <= 1)
	{
		double FM = 0.25 * ¦Ñ * c * (M + 1) * (M + 1);
		FL.f1 = FM;
		FL.f2 = FM * (¦Î1 * (-uc + 2 * c) / ¦Ã + u);
		FL.f3 = FM * (¦Î2 * (-uc + 2 * c) / ¦Ã + v);
		FL.f4 = FM * (uc * (-uc + 2 * c) / (¦Ã + 1) + 2 * c * c / (¦Ã * ¦Ã - 1) + 0.5 * (u * u + v * v) - ¦Î3 * (-uc + 2 * c) / ¦Ã);
	}
	FL.f1 = FL.f1 * D¦Î / J;
	FL.f2 = FL.f2 * D¦Î / J;
	FL.f3 = FL.f3 * D¦Î / J;
	FL.f4 = FL.f4 * D¦Î / J;
	return FL;
}
Flux VanLeerB(mesh & C, double ¦Îx, double ¦Îy, double ¦Ît, double J)
{
	double u = C.u;
	double v = C.v;
	double p = C.p;
	double ¦Ñ = C.¦Ñ;
	double c = sqrt(¦Ã * p / ¦Ñ);
	double D¦Î = sqrt(¦Îx * ¦Îx + ¦Îy * ¦Îy);
	double ub = u * ¦Îx + v * ¦Îy + ¦Ît;
	double uc = ub / D¦Î;
	double ¦Î1 = ¦Îx / D¦Î;
	double ¦Î2 = ¦Îy / D¦Î;
	double ¦Î3 = ¦Ît / D¦Î;
	Flux FR;
	if (uc - c >= 0)
	{
		FR = { 0 };
	}
	else if (uc + c <=0)
	{
		double E = 0.5 * ¦Ñ * (u * u + v * v) + p / (¦Ã - 1);
		FR.f1 = ¦Ñ * uc;
		FR.f2 = ¦Ñ * uc * u + p * ¦Î1;
		FR.f3 = ¦Ñ * uc * v + p * ¦Î2;
		FR.f4 = uc * (E + p) - p * ¦Î3;
	}
	double M = uc / c;
	if (abs(M) <= 1)
	{
		double FM = -0.25 * ¦Ñ * c * (M - 1) * (M - 1);
		FR.f1 = FM;
		FR.f2 = FM * (¦Î1 * (-uc - 2 * c) / ¦Ã + u);
		FR.f3 = FM * (¦Î2 * (-uc - 2 * c) / ¦Ã + v);
		FR.f4 = FM * (uc * (-uc - 2 * c) / (¦Ã + 1) + 2 * c * c / (¦Ã * ¦Ã - 1) + 0.5 * (u * u + v * v) - ¦Î3 * (-uc - 2 * c) / ¦Ã);
	}
	FR.f1 = FR.f1 * D¦Î / J;
	FR.f2 = FR.f2 * D¦Î / J;
	FR.f3 = FR.f3 * D¦Î / J;
	FR.f4 = FR.f4 * D¦Î / J;
	return FR;
}
