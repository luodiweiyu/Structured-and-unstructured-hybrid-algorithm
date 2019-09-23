#include<cmath>
#include"include/const.h"
#include<iostream>
#include<omp.h>
#include"include/functions.h"
using namespace std;
using namespace ConstPara;
Flux VanLeerA(mesh& C, double xix, double xiy, double xit,double J)
{
	double u = C.u;
	double v = C.v;
	double p = C.p;
	double rho = C.rho;
	double c = sqrt(gama * p / rho);
	double Dxi = sqrt(xix * xix + xiy * xiy);
	double ub = u * xix + v * xiy + xit;
	double uc = ub / Dxi;
	double xi1 = xix / Dxi;
	double xi2 = xiy / Dxi;
	double xi3 = xit / Dxi;
	Flux FL;
	if (uc + c <= 0)
	{
		FL = { 0 };
	}
	else if (uc - c >= 0)
	{
		double E = 0.5 * rho * (u * u + v * v) + p / (gama - 1);
		FL.f1 = rho * uc;
		FL.f2 = rho * uc * u + p * xi1;
		FL.f3 = rho * uc * v + p * xi2;
		FL.f4 = uc * (E + p) - p * xi3;
	}
	double M = uc / c;
	if (abs(M) <= 1)
	{
		double FM = 0.25 * rho * c * (M + 1) * (M + 1);
		FL.f1 = FM;
		FL.f2 = FM * (xi1 * (-uc + 2 * c) / gama + u);
		FL.f3 = FM * (xi2 * (-uc + 2 * c) / gama + v);
		FL.f4 = FM * (uc * (-uc + 2 * c) / (gama + 1) + 2 * c * c / (gama * gama - 1) + 0.5 * (u * u + v * v) - xi3 * (-uc + 2 * c) / gama);
	}
	FL.f1 = FL.f1 * Dxi / J;
	FL.f2 = FL.f2 * Dxi / J;
	FL.f3 = FL.f3 * Dxi / J;
	FL.f4 = FL.f4 * Dxi / J;
	return FL;
}
Flux VanLeerB(mesh & C, double xix, double xiy, double xit, double J)
{
	double u = C.u;
	double v = C.v;
	double p = C.p;
	double rho = C.rho;
	double c = sqrt(gama * p / rho);
	double Dxi = sqrt(xix * xix + xiy * xiy);
	double ub = u * xix + v * xiy + xit;
	double uc = ub / Dxi;
	double xi1 = xix / Dxi;
	double xi2 = xiy / Dxi;
	double xi3 = xit / Dxi;
	Flux FR;
	if (uc - c >= 0)
	{
		FR = { 0 };
	}
	else if (uc + c <=0)
	{
		double E = 0.5 * rho * (u * u + v * v) + p / (gama - 1);
		FR.f1 = rho * uc;
		FR.f2 = rho * uc * u + p * xi1;
		FR.f3 = rho * uc * v + p * xi2;
		FR.f4 = uc * (E + p) - p * xi3;
	}
	double M = uc / c;
	if (abs(M) <= 1)
	{
		double FM = -0.25 * rho * c * (M - 1) * (M - 1);
		FR.f1 = FM;
		FR.f2 = FM * (xi1 * (-uc - 2 * c) / gama + u);
		FR.f3 = FM * (xi2 * (-uc - 2 * c) / gama + v);
		FR.f4 = FM * (uc * (-uc - 2 * c) / (gama + 1) + 2 * c * c / (gama * gama - 1) + 0.5 * (u * u + v * v) - xi3 * (-uc - 2 * c) / gama);
	}
	FR.f1 = FR.f1 * Dxi / J;
	FR.f2 = FR.f2 * Dxi / J;
	FR.f3 = FR.f3 * Dxi / J;
	FR.f4 = FR.f4 * Dxi / J;
	return FR;
}
