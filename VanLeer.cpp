#include<cmath>
#include"const.h"
#include<iostream>
#include<omp.h>
#include"functions.h"
using namespace std;
using namespace ConstPara;
Flux VanLeerA(mesh& C, double ��x, double ��y, double ��t,double J)
{
	double u = C.u;
	double v = C.v;
	double p = C.p;
	double �� = C.��;
	double c = sqrt(�� * p / ��);
	double D�� = sqrt(��x * ��x + ��y * ��y);
	double ub = u * ��x + v * ��y + ��t;
	double uc = ub / D��;
	double ��1 = ��x / D��;
	double ��2 = ��y / D��;
	double ��3 = ��t / D��;
	Flux FL;
	if (uc + c <= 0)
	{
		FL = { 0 };
	}
	else if (uc - c >= 0)
	{
		double E = 0.5 * �� * (u * u + v * v) + p / (�� - 1);
		FL.f1 = �� * uc;
		FL.f2 = �� * uc * u + p * ��1;
		FL.f3 = �� * uc * v + p * ��2;
		FL.f4 = uc * (E + p) - p * ��3;
	}
	double M = uc / c;
	if (abs(M) <= 1)
	{
		double FM = 0.25 * �� * c * (M + 1) * (M + 1);
		FL.f1 = FM;
		FL.f2 = FM * (��1 * (-uc + 2 * c) / �� + u);
		FL.f3 = FM * (��2 * (-uc + 2 * c) / �� + v);
		FL.f4 = FM * (uc * (-uc + 2 * c) / (�� + 1) + 2 * c * c / (�� * �� - 1) + 0.5 * (u * u + v * v) - ��3 * (-uc + 2 * c) / ��);
	}
	FL.f1 = FL.f1 * D�� / J;
	FL.f2 = FL.f2 * D�� / J;
	FL.f3 = FL.f3 * D�� / J;
	FL.f4 = FL.f4 * D�� / J;
	return FL;
}
Flux VanLeerB(mesh & C, double ��x, double ��y, double ��t, double J)
{
	double u = C.u;
	double v = C.v;
	double p = C.p;
	double �� = C.��;
	double c = sqrt(�� * p / ��);
	double D�� = sqrt(��x * ��x + ��y * ��y);
	double ub = u * ��x + v * ��y + ��t;
	double uc = ub / D��;
	double ��1 = ��x / D��;
	double ��2 = ��y / D��;
	double ��3 = ��t / D��;
	Flux FR;
	if (uc - c >= 0)
	{
		FR = { 0 };
	}
	else if (uc + c <=0)
	{
		double E = 0.5 * �� * (u * u + v * v) + p / (�� - 1);
		FR.f1 = �� * uc;
		FR.f2 = �� * uc * u + p * ��1;
		FR.f3 = �� * uc * v + p * ��2;
		FR.f4 = uc * (E + p) - p * ��3;
	}
	double M = uc / c;
	if (abs(M) <= 1)
	{
		double FM = -0.25 * �� * c * (M - 1) * (M - 1);
		FR.f1 = FM;
		FR.f2 = FM * (��1 * (-uc - 2 * c) / �� + u);
		FR.f3 = FM * (��2 * (-uc - 2 * c) / �� + v);
		FR.f4 = FM * (uc * (-uc - 2 * c) / (�� + 1) + 2 * c * c / (�� * �� - 1) + 0.5 * (u * u + v * v) - ��3 * (-uc - 2 * c) / ��);
	}
	FR.f1 = FR.f1 * D�� / J;
	FR.f2 = FR.f2 * D�� / J;
	FR.f3 = FR.f3 * D�� / J;
	FR.f4 = FR.f4 * D�� / J;
	return FR;
}
