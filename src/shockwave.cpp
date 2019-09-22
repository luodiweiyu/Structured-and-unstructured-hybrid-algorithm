#include"/Structured-and-unstructured-hybrid-algorithm/include/const.h"
#include<cmath>
#include"/Structured-and-unstructured-hybrid-algorithm/include/shockwave.h"
//���㼤����ز���Ǯ��𢡶��������ѧ��p228-241
using namespace ConstPara;
double get_c(double ��, double p)//���ٹ�ʽ
{
	return sqrt(gama * p / ��);
}
double get_Ma(double u, double v, double ��, double p)//��������
{
	double a = get_c(��, p);
	double velocity = sqrt(u * u + v * v);
	return velocity / a;
}
//��������ϵʽ
//7-99
double get_Mas(double p1, double p2)
{
	double Mas = 1 + (gama + 1) * (p2 / p1 - 1) / (2 * gama);
	return sqrt(Mas);
}
//7-105
double get_p2p1(double Ma1)
{
	return  2 * gama* Ma1* Ma1 / (gama + 1) - (gama - 1) / (gama + 1);
}
//7-106
double get_��2��1(double Ma1)
{
	return  (gama + 1) * Ma1 * Ma1 / ((gama - 1) * Ma1 * Ma1 + 2);
}
double get_Ma2(double Ma1)
{
	return sqrt((Ma1 * Ma1 + 2 / (gama - 1)) / (2 * gama * Ma1 * Ma1 / (gama - 1) - 1));
}
//������������  �޶���,��ʿ����p57 3.5
double get_Mu(mesh & U, mesh & D, double ��)
{
	using namespace ConstPara;
	double �� = 1e-10;
	double Mu = 2, Mu1 = 20;
	double udn, uun, udt, uut;
	uun = get_un(U, ��);
	udn = get_un(D, ��);
	uut = get_ut(U, ��);
	udt = get_ut(D, ��);
	double cu = sqrt(gama * U.p / U.��);
	double J, Fx, fx;
	double a, b;
	while (abs(Mu - Mu1) >= ��)
	{
		Mu = Mu1;
		double J = 2 * sqrt(gama * D.p / D.��) / (gama - 1) - udn;
		a = 2 * gama * Mu1 * Mu1 - (gama - 1);
		b = (gama - 1) * Mu1 * Mu1 + 2;
		Fx = sqrt(a * b) / (Mu1 * (gama - 1)) + (Mu1 * Mu1 - 1) / Mu1 - (gama + 1) * (J + uun) / (2 * cu);
		fx = (4 * gama * Mu1 * b + 2 * (gama - 1) * Mu1 * a) / (2 * Mu1 * (gama - 1) * sqrt(a * b)) - sqrt(a * b) / (Mu1 * Mu1 * (gama - 1)) + 1 + 1 / (Mu1 * Mu1);
		Mu1 = Mu1 - Fx / fx;
	}
	return Mu;
}
//3.4
double get_udn(mesh & U, mesh & D, double Ma1, double Vs, double ��)
{
	double ��1 = get��fromMa(Ma1);
	double ��2 = 1 / ��1;
	double Ma2 = getMafrom��(��2);
	double c2 = sqrt(gama * D.p / D.��);
	double V2 = Ma2 * c2;
	return V2 + Vs;
	//double a = ((gama - 1) * Ma1 + 2) / (2 * gama * Ma1 - (gama - 1));
	//double b = (a * gama * D.p) / D.��;
	//double cu = sqrt(gama * U.p / U.��);
	//double uun = U.u * sin(��) - U.v * cos(��);
	//double udn = sqrt(b) - uun +cu * Ma1 ;
	//return udn;
}

double get��fromMa(double Ma)
{
	double a = 2 + (gama - 1) * Ma * Ma;
	double b = 0;
	double c = -(gama + 1) * Ma * Ma;
	return (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
}
double getMafrom��(double ��)
{
	double ��2 = �� * ��;
	return sqrt((2 * ��2 / (gama + 1)) / (1 - (gama - 1) / (gama + 1) * ��2));
}

//����Ϊб����,����ֻ֤�ʺϲ�ǰ�ٶ�ˮƽ��v=0,������ͨ����
double get_Ma2(double Ma1, double beta)//б����������������Ǯ��𢡶��������ѧ��p241,7-131
{
	double Ma2 = (Ma1 * Ma1 + 2 / (gama - 1)) / (2 * gama / (gama - 1) * Ma1 * Ma1 * sin(beta) * sin(beta) - 1) + 2 / (gama - 1) * Ma1 * Ma1 * cos(beta) * cos(beta) / (Ma1 * Ma1 * sin(beta) * sin(beta) + 2 / (gama - 1));
	Ma2 = sqrt(Ma2);
	return Ma2;
}
//double get_��(double Ma1, double beta)//���������
//{
//
//	double tan�� = (Ma1 * Ma1 * sin(beta) * sin(beta) - 1) / ((Ma1 * Ma1 * ((gama + 1) / 2 - sin(beta) * sin(beta)) + 1) * tan(beta));
//	double �� = atan(tan��);
//	//if (�� < 0)
//	//	�� = �� + pi;
//	return ��;
//}
double get_betafrom��(double Ma1, double ��)
{
	//�ҽط� �������ѧ�뼼��p132
	double beta = 40 * pi / 180, beta1 = 50 * pi / 180, beta0 = 60 * pi / 180;
	double fbeta1, fbeta0;
	if (�� == 0.92846675631531839)
		�� = 0.92846675631531839;
	while (abs(beta - beta1) > 1e-10)
	{
		beta0 = beta1;
		beta1 = beta;
		while (beta1 > pi / 2 || beta1 < -pi / 2)
		{
			if (beta1 > pi / 2)
				beta1 = beta1 - pi / 2;
			else if (beta1 < -pi / 2)
				beta1 = beta1 + pi / 2;
		}
		while (beta0 > pi / 2 || beta0 < -pi / 2)
		{
			if (beta0 > pi / 2)
				beta0 = beta0 - pi / 2;
			else if (beta0 < -pi / 2)
				beta0 = beta0 + pi / 2;
		}

		fbeta1 = tan(��) - (Ma1 * Ma1 * sin(beta1) * sin(beta1) - 1) / ((Ma1 * Ma1 * ((gama + 1) / 2 - sin(beta1) * sin(beta1)) + 1) * tan(beta1));
		fbeta0 = tan(��) - (Ma1 * Ma1 * sin(beta0) * sin(beta0) - 1) / ((Ma1 * Ma1 * ((gama + 1) / 2 - sin(beta0) * sin(beta0)) + 1) * tan(beta0));
		beta = beta1 - fbeta1 / (fbeta1 - fbeta0) * (beta1 - beta0);
	}
	return beta;
}
double get_ufromMa2(double Ma2, double ��2, double p2, double ��)
{
	double a2 = sqrt(gama * p2 / ��2);
	double velocity = Ma2 * a2;
	double u = velocity * cos(��);
	return u;
}
double get_vfromMa2(double Ma2, double ��2, double p2, double ��)
{
	double a2 = sqrt(gama * p2 / ��2);
	double velocity = Ma2 * a2;
	double v = velocity * sin(��);
	return v;
}
double get_p2(double Ma1, double beta, double p1)
{
	double p2 = (2 * gama / (gama + 1) * Ma1 * Ma1 * sin(beta) * sin(beta) - (gama - 1) / (gama + 1)) * p1;
	return p2;
}
double get_��2(double Ma1, double beta, double ��1)
{
	double ��2 = ((gama + 1) * Ma1 * Ma1 * sin(beta) * sin(beta) / ((gama - 1) * Ma1 * Ma1 * sin(beta) * sin(beta) + 2)) * ��1;
	return ��2;
}


//����Ϊͨ����б�����㷨
double get_un(mesh & A, double beta)
{
	if (beta >= 0)
		return A.u * sin(beta) - A.v * cos(beta);
	else
		return A.u * sin(-beta) + A.v * cos(-beta);
}
double get_ut(mesh & A, double beta)
{
	if (beta >= 0)
		return A.u * cos(beta) + A.v * sin(beta);
	else
		return A.u * cos(-beta) - A.v * sin(-beta);
}

void get_down(mesh & U, mesh & D, double beta)//���β���
{
	double udn, uun, udt, uut;
	uun = get_un(U, beta);
	uut = get_ut(U, beta);
	double au = get_c(U.��, U.p);
	double Mu = uun / au;
	D.p = U.p * get_p2p1(Mu);
	D.�� = U.�� * get_��2��1(Mu);
	double Md = get_Ma2(Mu);
	double ad = get_c(D.��, D.p);
	udn = Md * ad;
	udt = uut;
	if (beta >= 0)
	{
		D.u = udt * cos(beta) + udn * sin(beta);
		D.v = udt * sin(beta) - udn * cos(beta);
	}
	else
	{
		D.u = udt * cos(-beta) + udn * sin(-beta);
		D.v = -udt * sin(-beta) + udn * cos(-beta);
	}
}
double get_��(double u, double v)
{
	double �� = atan(abs(v) / abs(u));
	if (u > 0)
	{
		if (v > 0)
			return ��;
		else if (v < 0)
			return -��;
		else
			return 0;
	}
	else if (u < 0)
	{
		if (v > 0)
			return -��;
		else if (v < 0)
			return ��;
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
double get_beta(mesh & U, double p2, int type)//type=-1��ʾ�����beta��Ϊ��������Ϊ��
{
	double beta = 40 * pi / 180, beta1 = 50 * pi / 180, beta0 = 60 * pi / 180;
	double uun0, uun1;
	double Mu0, Mu1;
	double p20, p21;
	double fbeta1, fbeta0;
	while (abs(beta - beta1) > 1e-15)
	{
		while (beta > pi / 2 || beta < -pi / 2)
		{
			if (beta > pi / 2)
				beta = beta - pi / 2;
			else if (beta < -pi / 2)
				beta = beta + pi / 2;
		}
		beta0 = beta1;
		beta1 = beta;
		if (type == -1)
		{
			if (beta0 > 0)
				beta0 = -beta0;
			if (beta1 > 0)
				beta1 = -beta1;
		}
		else if (type == 1)
		{
			if (beta0 < 0)
				beta0 = -beta0;
			if (beta1 < 0)
				beta1 = -beta1;
		}
		uun0 = get_un(U, beta0);
		uun1 = get_un(U, beta1);
		Mu0 = get_Ma(uun0, 0, U.��, U.p);
		Mu1 = get_Ma(uun1, 0, U.��, U.p);
		p20 = U.p * get_p2p1(Mu0);
		p21 = U.p * get_p2p1(Mu1);
		fbeta0 = p20 - p2;
 		fbeta1 = p21 - p2;
		if (fbeta0 == fbeta1)
			break;
		else
			beta = beta1 - fbeta1 / (fbeta1 - fbeta0) * (beta1 - beta0);
	}
	return beta;

}