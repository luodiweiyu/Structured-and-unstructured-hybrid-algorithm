#include<cmath>
#include"const.h"
#include"shockwave.h"
using ConstPara::gama;
double getδfromλ(double λ)
{
	double k = (gama - 1) / (gama + 1);
	double λ2 = λ * λ;
	return sqrt(1 / k)*atan(sqrt((k*(λ2 - 1)) / (1 - k * λ2))) - atan(sqrt((λ2 - 1) / (1 - k * λ2)));
}
double getμfromMa(double Ma)
{
	return asin(1 / Ma);
}
double getθfromMa(double Ma)
{
	double λ = getλfromMa(Ma);
	return sqrt((gama + 1) / (gama - 1))*asin(sqrt((gama - 1)*(λ*λ - 1) / 2));
}
double getλfromδ(double δ)
{
	double Δ = 1e-6;
	double λ1 = 1, λ2 = sqrt((gama + 1) / (gama - 1));
	double λM = (λ1 + λ2) / 2;
	double δ1, δM=0, δ2;
	while (abs(δM - δ) > 1e-6)
	{
		λM = (λ1 + λ2) / 2;
		δ1 = getδfromλ(λ1);
		δM = getδfromλ(λM);
		δ2 = getδfromλ(λ2);
		if ((δ1 - δ)*(δM - δ) < 0)
			λ1 = λ1, λ2 = λM;
		else
			λ1 = λM, λ2 = λ2;
	}
	return (λ1 + λ2) / 2;
}
//以下为等熵流公式，可以用于膨胀波计算
double getTfromλ(double λ,double T0)
{
	double temp = 1 - (gama - 1)*λ*λ / (gama + 1);
	return temp*T0;
}
double getT0fromλandT(double λ, double T)
{
	double temp = 1 - (gama - 1)*λ*λ / (gama + 1);
	return T / temp;
}
double getpfromλ(double λ, double p0)
{
	double temp = 1 - (gama - 1)*λ*λ / (gama + 1);
	temp = pow(temp, gama / (gama - 1));
	return temp * p0;
}
double getp0fromλandp(double λ, double p)
{
	double temp = 1 - (gama - 1)*λ*λ / (gama + 1);
	temp = pow(temp, gama / (gama - 1));
	return p / temp;
}
double getρfromλ(double λ, double ρ0)
{
	double temp = 1 - (gama - 1)*λ*λ / (gama + 1);
	temp = pow(temp, 1 / (gama - 1)); 
	return temp * ρ0;
}
double getρ0fromλandρ(double λ, double ρ)
{
	double temp = 1 - (gama - 1)*λ*λ / (gama + 1);
	temp = pow(temp, 1 / (gama - 1));
	return ρ / temp;
}
