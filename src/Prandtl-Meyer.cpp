#include<cmath>
#include"/Structured-and-unstructured-hybrid-algorithm/include/const.h"
#include"/Structured-and-unstructured-hybrid-algorithm/include/shockwave.h"
using ConstPara::gama;
double getdeltafromlambda(double lambda)
{
	double k = (gama - 1) / (gama + 1);
	double lambda2 = lambda * lambda;
	return sqrt(1 / k)*atan(sqrt((k*(lambda2 - 1)) / (1 - k * lambda2))) - atan(sqrt((lambda2 - 1) / (1 - k * lambda2)));
}
double getmufromMa(double Ma)
{
	return asin(1 / Ma);
}
double getthetafromMa(double Ma)
{
	double lambda = getlambdafromMa(Ma);
	return sqrt((gama + 1) / (gama - 1))*asin(sqrt((gama - 1)*(lambda*lambda - 1) / 2));
}
double getlambdafromdelta(double delta)
{
	double m = 1e-6;
	double lambda1 = 1, lambda2 = sqrt((gama + 1) / (gama - 1));
	double lambdaM = (lambda1 + lambda2) / 2;
	double delta1, deltaM=0, delta2;
	while (abs(deltaM - delta) > 1e-6)
	{
		lambdaM = (lambda1 + lambda2) / 2;
		delta1 = getdeltafromlambda(lambda1);
		deltaM = getdeltafromlambda(lambdaM);
		delta2 = getdeltafromlambda(lambda2);
		if ((delta1 - delta)*(deltaM - delta) < 0)
			lambda1 = lambda1, lambda2 = lambdaM;
		else
			lambda1 = lambdaM, lambda2 = lambda2;
	}
	return (lambda1 + lambda2) / 2;
}
//以下为等熵流公式，可以用于膨胀波计算
double getTfromlambda(double lambda,double T0)
{
	double temp = 1 - (gama - 1)*lambda*lambda / (gama + 1);
	return temp*T0;
}
double getT0fromlambdaandT(double lambda, double T)
{
	double temp = 1 - (gama - 1)*lambda*lambda / (gama + 1);
	return T / temp;
}
double getpfromlambda(double lambda, double p0)
{
	double temp = 1 - (gama - 1)*lambda*lambda / (gama + 1);
	temp = pow(temp, gama / (gama - 1));
	return temp * p0;
}
double getp0fromlambdaandp(double lambda, double p)
{
	double temp = 1 - (gama - 1)*lambda*lambda / (gama + 1);
	temp = pow(temp, gama / (gama - 1));
	return p / temp;
}
double getrhofromlambda(double lambda, double rho0)
{
	double temp = 1 - (gama - 1)*lambda*lambda / (gama + 1);
	temp = pow(temp, 1 / (gama - 1)); 
	return temp * rho0;
}
double getrho0fromlambdaandrho(double lambda, double rho)
{
	double temp = 1 - (gama - 1)*lambda*lambda / (gama + 1);
	temp = pow(temp, 1 / (gama - 1));
	return rho / temp;
}
