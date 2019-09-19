#include<cmath>
#include"const.h"
#include<iostream>
#include<omp.h>
#include"functions.h"
using namespace std;
using namespace ConstPara;
//2维欧拉方程组离散：张德良《计算流体力学教程》p433 13.3.7
Flux HLLC_Χ(mesh& CL, mesh& CR, mesh& C, int method)//半点左侧右侧格点，坐标变换参考点
{
	double ξx = C.ξx[method];
	double ξy = C.ξy[method];
	double ξt = C.ξt[method];
	double J = C.J[method];
	double Dξ = sqrt(ξx * ξx + ξy * ξy);
	double ξ1 = ξx / Dξ;
	double ξ2 = ξy / Dξ;
	double ξ3 = ξt / Dξ;
	double ξcL = CL.u * ξ1 + CL.v * ξ2 + ξ3;
	double ξcR = CR.u * ξ1 + CR.v * ξ2 + ξ3;
	double aL = sqrt(γ * CL.p / CL.ρ);
	double aR = sqrt(γ * CR.p / CR.ρ);
	double a = sqrt(CL.ρ);
	double b = sqrt(CR.ρ);

	double SL = ξcL - aL;
	double SR = ξcR + aR;

	double EL = CL.p / (γ - 1) + 0.5 * CL.ρ * CL.u * CL.u + 0.5 * CL.ρ * CL.v * CL.v;
	double ER = CR.p / (γ - 1) + 0.5 * CR.ρ * CR.u * CR.u + 0.5 * CR.ρ * CR.v * CR.v;
	Flux FL, FR;
	FL.f1 = (CL.ρ * ξcL);
	FL.f2 = (CL.ρ * CL.u * ξcL + ξ1 * CL.p);
	FL.f3 = (CL.ρ * CL.v * ξcL + ξ2 * CL.p);
	FL.f4 = (ξcL * (EL + CL.p) - ξ3 * CL.p);
	FR.f1 = (CR.ρ * ξcR);
	FR.f2 = (CR.ρ * CR.u * ξcR + ξ1 * CR.p);
	FR.f3 = (CR.ρ * CR.v * ξcR + ξ2 * CR.p);
	FR.f4 = (ξcR * (ER + CR.p) - ξ3 * CR.p);
	double SM = (CL.p - CR.p + CL.ρ * ξcL * (ξcL - SL) + CR.ρ * ξcR * (SR - ξcR)) / (CR.ρ * (SR - ξcR) + CL.ρ * (ξcL - SL));
	double pM = 0.5 * (CL.ρ * (ξcL - SL) * (ξcL - SM) + CR.ρ * (ξcR - SR) * (ξcR - SM) + CR.p + CL.p);
	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
	double par1 = (ξcL - SL) / (SM - SL);
	double par2 = (pM - CL.p) / (SM - SL);
	Flux  FML;
	FML.f1 = (SM * par1 * CL.ρ);
	FML.f2 = (SM * par1 * CL.ρ * CL.u - par2 * ξ1 + ξ1 * pM);
	FML.f3 = (SM * par1 * CL.ρ * CL.v - par2 * ξ2 + ξ2 * pM);
	FML.f4 = (SM * par1 * (EL + CL.p) + par2 * ξ3 - ξ3 * pM);
	double par3 = (ξcR - SR) / (SM - SR);
	double par4 = (pM - CR.p) / (SM - SR);
	Flux FMR;
	FMR.f1 = (SM * par3 * CR.ρ);
	FMR.f2 = (SM * par3 * CR.ρ * CR.u - par4 * ξ1 + ξ1 * pM);
	FMR.f3 = (SM * par3 * CR.ρ * CR.v - par4 * ξ2 + ξ2 * pM);
	FMR.f4 = (SM * par3 * (ER + CR.p) + par4 * ξ3 - ξ3 * pM);

	//double ΩL = 1 / (SL - SM);
	//double ΩR = 1 / (SR - SM);

	//double ρLM = ΩL * CL.ρ * (SL - ξcL);
	//double ρuLS = ΩL * ((SL - ξcL) * (CL.ρ * CL.u) + (pM - CL.p) * ξ1);
	//double ρvLS = ΩL * ((SL - ξcL) * (CL.ρ * CL.v) + (pM - CL.p) * ξ2);
	//double eLS = ΩL * ((SL - ξcL) * EL - CL.p * ξcL + pM * SM - (pM - CL.p) * ξ3);

	//double ρRM = ΩR * CR.ρ * (SR - ξcR);
	//double ρuRS = ΩR * ((SR - ξcR) * (CR.ρ * CR.u) + (pM - CR.p) * ξ1);
	//double ρvRS = ΩR * ((SR - ξcR) * (CR.ρ * CR.v) + (pM - CR.p) * ξ2);
	//double eRS = ΩR * ((SR - ξcR) * ER - CR.p * ξcR + pM * SM - (pM - CR.p) * ξ3);

	//Flux FML, FMR;
	//FML.f1 = SM * ρLM;
	//FML.f2 = ρuLS * SM + pM * ξ1;
	//FML.f3 = ρvLS * SM + pM * ξ2;
	//FML.f4 = (eLS + pM) * SM - pM * ξ3;
	//FMR.f1 = SM * ρRM;
	//FMR.f2 = ρuRS * SM + pM * ξ1;
	//FMR.f3 = ρvRS * SM + pM * ξ2;
	//FMR.f4 = (eRS + pM) * SM - pM * ξ3;

	Flux F_HLLC;
	if (SL > 0)
		F_HLLC = FL;
	else if (SL <= 0 && SM > 0)
		F_HLLC = FML;
	else if (SM <= 0 && SR >= 0)
		F_HLLC = FMR;
	else
		F_HLLC = FR;
	F_HLLC.f1 = F_HLLC.f1 * Dξ;
	F_HLLC.f2 = F_HLLC.f2 * Dξ;
	F_HLLC.f3 = F_HLLC.f3 * Dξ;
	F_HLLC.f4 = F_HLLC.f4 * Dξ;

	return F_HLLC;
}

Flux HLLC_Υ(mesh & CD, mesh & CU, mesh & C, int method)
{
	double ηx = C.ηx[method];
	double ηy = C.ηy[method];
	double ηt = C.ηt[method];
	double J = C.J[method];

	double Dη = sqrt(ηx * ηx + ηy * ηy);
	double η1 = ηx / Dη;
	double η2 = ηy / Dη;
	double η3 = ηt / Dη;
	double ηcU = CU.u * η1 + CU.v * η2;
	double ηcD = CD.u * η1 + CD.v * η2;
	double aU = sqrt(γ * CU.p / CU.ρ);
	double aD = sqrt(γ * CD.p / CD.ρ);
	double EU = CU.p / (γ - 1) + 0.5 * CU.ρ * CU.u * CU.u + 0.5 * CU.ρ * CU.v * CU.v;
	double ED = CD.p / (γ - 1) + 0.5 * CD.ρ * CD.u * CD.u + 0.5 * CD.ρ * CD.v * CD.v;
	Flux GD, GU;
	GD.f1 = (CD.ρ * ηcD);
	GD.f2 = (CD.ρ * CD.u * ηcD + η1 * CD.p);
	GD.f3 = (CD.ρ * CD.v * ηcD + η2 * CD.p);
	GD.f4 = (ηcD * (ED + CD.p) - η3 * CD.p);

	GU.f1 = (CU.ρ * ηcU);
	GU.f2 = (CU.ρ * CU.u * ηcU + η1 * CD.p);
	GU.f3 = (CU.ρ * CU.v * ηcU + η2 * CD.p);
	GU.f4 = (ηcU * (EU + CU.p) - η3 * CD.p);

	double SD = ηcD - aD;
	double SU = ηcU + aU;
	double SM = (CD.p - CU.p - CD.ρ * ηcD * (SD - ηcD) + CU.ρ * ηcU * (SU - ηcU)) / (CU.ρ * (SU - ηcU) - CD.ρ * (SD - ηcD));
	double pM = 0.5 * (CD.ρ * (ηcD - SD) * (ηcD - SM) + CU.ρ * (ηcU - SU) * (ηcU - SM) + CU.p + CD.p);
	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
	double par1 = (ηcD - SD) / (SM - SD);
	double par2 = (pM - CD.p) / (SM - SD);
	Flux  FMD;
	FMD.f1 = (SM * par1 * CD.ρ);
	FMD.f2 = (SM * par1 * CD.ρ * CD.u - par2 * η1 + η1 * pM);
	FMD.f3 = (SM * par1 * CD.ρ * CD.v - par2 * η2 + η2 * pM);
	FMD.f4 = (SM * par1 * (ED + CD.p) + par2 * η3 - η3 * pM);
	double par3 = (ηcU - SU) / (SM - SU);
	double par4 = (pM - CU.p) / (SM - SU);
	Flux FMU;
	FMU.f1 = (SM * par3 * CU.ρ);
	FMU.f2 = (SM * par3 * CU.ρ * CU.u - par4 * η1 + η1 * pM);
	FMU.f3 = (SM * par3 * CU.ρ * CU.v - par4 * η2 + η2 * pM);
	FMU.f4 = (SM * par3 * (EU + CU.p) + par4 * η3 - η3 * pM);

	//double ΩD = 1 / (SD - SM);
	//double ΩU = 1 / (SU - SM);

	//double ρDM = ΩD * CD.ρ * (SD - ηcD);
	//double ρuDS = ΩD * ((SD - ηcD) * (CD.ρ * CD.u) + (pM - CD.p) * η1);
	//double ρvDS = ΩD * ((SD - ηcD) * (CD.ρ * CD.v) + (pM - CD.p) * η2);
	//double eDS = ΩD * ((SD - ηcD) * ED - CD.p * ηcD + pM * SM - (pM - CD.p) * η3);

	//double ρUM = ΩU * CU.ρ * (SU - ηcU);
	//double ρuUS = ΩU * ((SU - ηcU) * (CU.ρ * CU.u) + (pM - CU.p) * η1);
	//double ρvUS = ΩU * ((SU - ηcU) * (CU.ρ * CU.v) + (pM - CU.p) * η2);
	//double eUS = ΩU * ((SU - ηcU) * EU - CU.p * ηcU + pM * SM - (pM - CU.p) * η3);

	//Flux  FMD, FMU;

	//FMD.f1 = SM * ρDM;
	//FMD.f2 = ρuDS * SM + pM * η1;
	//FMD.f3 = ρvDS * SM + pM * η2;
	//FMD.f4 = (eDS + pM) * SM - pM * η3;
	//FMU.f1 = SM * ρUM;
	//FMU.f2 = ρuUS * SM + pM * η1;
	//FMU.f3 = ρvUS * SM + pM * η2;
	//FMU.f4 = (eUS + pM) * SM - pM * η3;

	Flux G_HLLC;
	if (SD >= 0)
		G_HLLC = GD;
	else if (SD <= 0 && SM >= 0)
		G_HLLC = FMD;
	else if (SM <= 0 && SU >= 0)
		G_HLLC = FMU;
	else
		G_HLLC = GU;
	G_HLLC.f1 = G_HLLC.f1 * Dη;
	G_HLLC.f2 = G_HLLC.f2 * Dη;
	G_HLLC.f3 = G_HLLC.f3 * Dη;
	G_HLLC.f4 = G_HLLC.f4 * Dη;

	return G_HLLC;
}

Flux HLLC_Χ2(mesh& CL, mesh& CR, mesh& C, int method)//半点左侧右侧格点，坐标变换参考点
{
	double ξx = C.ξx[method];
	double ξy = C.ξy[method];
	double ξt = C.ξt[method];
	double J = C.J[method];
	double Dξ = sqrt(ξx * ξx + ξy * ξy);
	double ξ1 = ξx / Dξ;
	double ξ2 = ξy / Dξ;
	double ξ3 = ξt / Dξ;
	double aL = sqrt(γ * CL.p / CL.ρ);
	double aR = sqrt(γ * CR.p / CR.ρ);
	double a = sqrt(CL.ρ);
	double b = sqrt(CR.ρ);
	double ubl = CL.u * ξ1 + CL.v * ξ2 + ξ3;
	double ubr = CR.u * ξ1 + CR.v * ξ2 + ξ3;
	double ubp = (a * ubl + b * ubr) / (a + b);
	double aM = (a * aL * aL + b * aR * aR) / (a + b) + 0.5 * (a * b / ((a + b) * (a + b))) * ((CL.u - CR.u) * (CL.u - CR.u) + (CL.v - CR.v) * (CL.v - CR.v));


	double conl = min(ubl - aL, ubp - aM);
	double conr = max(ubr + aR, ubp + aM);
	double s1 = CR.ρ * (conr - ubr);
	double s2 = CL.ρ * (ubl - conl);
	double ss = (s1 * ubr + s2 * ubl + CL.p - CR.p) / (s1 + s2);
	double ps = s2 * (ubl - ss) + CL.p;
	Flux FL, FR;
	FL.f1 = CL.ρ;
	FL.f2 = CL.ρ * CL.u;
	FL.f3 = CL.ρ * CL.v;
	FL.f4 = CL.p * γ / (γ - 1) + 0.5 * CL.ρ * CL.u * CL.u + 0.5 * CL.ρ * CL.v * CL.v; 
	FR.f1 = CR.ρ;
	FR.f2 = CR.ρ * CR.u;
	FR.f3 = CR.ρ * CR.v;
	FR.f4 = CR.p * γ / (γ - 1) + 0.5 * CR.ρ * CR.u * CR.u + 0.5 * CR.ρ * CR.v * CR.v;
	double par1 = (ubl - conl) / (ss - conl);
	double par2 = (ps - CL.p) / (ss - conl);

	Flux FML, FMR;
	FML.f1 = par1 * FL.f1;
	FML.f2 = par1 * FL.f2 - par2 * ξ1;
	FML.f3 = par1 * FL.f3 - par2 * ξ2;
	FML.f4 = par1 * FL.f4 + par2 * ξ3;
	double par3 = (ubr - conr) / (ss - conr);
	double par4 = (ps - CR.p) / (ss - conr);

	FMR.f1 = par3 * FR.f1;
	FMR.f2 = par3 * FR.f2 - par4 * ξ1;
	FMR.f3 = par3 * FR.f3 - par4 * ξ2;
	FMR.f4 = par3 * FR.f4 + par4 * ξ3;

	Flux F_HLLC;
	if (conl > 0)
	{
		F_HLLC.f1 = FL.f1 * ubl;
		F_HLLC.f2 = FL.f2 * ubl + ξ1 * CL.p;
		F_HLLC.f3 = FL.f3 * ubl + ξ2 * CL.p;
		F_HLLC.f4 = FL.f4 * ubl - ξ3 * CL.p;
	}
	else if (conl <= 0 && ss > 0)
	{
		F_HLLC.f1 = FML.f1 * ss;
		F_HLLC.f2 = FML.f2 * ss + ξ1 * ps;
		F_HLLC.f3 = FML.f3 * ss + ξ2 * ps;
		F_HLLC.f4 = FML.f4 * ss - ξ3 * ps;
	}
	else if (ss <= 0 && conr >= 0)
	{
		F_HLLC.f1 = FMR.f1 * ss;
		F_HLLC.f2 = FMR.f2 * ss + ξ1 * ps;
		F_HLLC.f3 = FMR.f3 * ss + ξ2 * ps;
		F_HLLC.f4 = FMR.f4 * ss - ξ3 * ps;
	}
	else
	{
		F_HLLC.f1 = FR.f1 * ubr;
		F_HLLC.f2 = FR.f2 * ubr + ξ1 * CR.p;
		F_HLLC.f3 = FR.f3 * ubr + ξ2 * CR.p;
		F_HLLC.f4 = FR.f4 * ubr - ξ3 * CR.p;
	}
	F_HLLC.f1 = F_HLLC.f1 * Dξ;
	F_HLLC.f2 = F_HLLC.f2 * Dξ;
	F_HLLC.f3 = F_HLLC.f3 * Dξ;
	F_HLLC.f4 = F_HLLC.f4 * Dξ;

	return F_HLLC;
}
Flux HLLC_Υ2(mesh& CL, mesh& CR, mesh& C, int method)
{
	double ηx = C.ηx[method];
	double ηy = C.ηy[method];
	double ηt = C.ηt[method];
	double J = C.J[method];
	double Dη = sqrt(ηx * ηx + ηy * ηy);
	double η1 = ηx / Dη;
	double η2 = ηy / Dη;
	double η3 = ηt / Dη;
	double aL = sqrt(γ * CL.p / CL.ρ);
	double aR = sqrt(γ * CR.p / CR.ρ);
	double a = sqrt(CL.ρ);
	double b = sqrt(CR.ρ);
	double ubl = CL.u * η1 + CL.v * η2 + η3;
	double ubr = CR.u * η1 + CR.v * η2 + η3;
	double ubp = (a * ubl + b * ubr) / (a + b);
	double aM = (a * aL * aL + b * aR * aR) / (a + b) + 0.5 * (a * b / ((a + b) * (a + b))) * ((CL.u - CR.u) * (CL.u - CR.u) + (CL.v - CR.v) * (CL.v - CR.v));


	double conl = min(ubl - aL, ubp - aM);
	double conr = max(ubr + aR, ubp + aM);
	double s1 = CR.ρ * (conr - ubr);
	double s2 = CL.ρ * (ubl - conl);
	double ss = (s1 * ubr + s2 * ubl + CL.p - CR.p) / (s1 + s2);
	double ps = s2 * (ubl - ss) + CL.p;
	Flux FL, FR;
	FL.f1 = CL.ρ;
	FL.f2 = CL.ρ * CL.u;
	FL.f3 = CL.ρ * CL.v;
	FL.f4 = CL.p * γ / (γ - 1) + 0.5 * CL.ρ * CL.u * CL.u + 0.5 * CL.ρ * CL.v * CL.v; 
	FR.f1 = CR.ρ;
	FR.f2 = CR.ρ * CR.u;
	FR.f3 = CR.ρ * CR.v;
	FR.f4 = CR.p * γ / (γ - 1) + 0.5 * CR.ρ * CR.u * CR.u + 0.5 * CR.ρ * CR.v * CR.v;
	double par1 = (ubl - conl) / (ss - conl);
	double par2 = (ps - CL.p) / (ss - conl);

	Flux FML, FMR;
	FML.f1 = par1 * FL.f1;
	FML.f2 = par1 * FL.f2 - par2 * η1;
	FML.f3 = par1 * FL.f3 - par2 * η2;
	FML.f4 = par1 * FL.f4 + par2 * η3;
	double par3 = (ubr - conr) / (ss - conr);
	double par4 = (ps - CR.p) / (ss - conr);

	FMR.f1 = par3 * FR.f1;
	FMR.f2 = par3 * FR.f2 - par4 * η1;
	FMR.f3 = par3 * FR.f3 - par4 * η2;
	FMR.f4 = par3 * FR.f4 + par4 * η3;

	Flux F_HLLC;
	if (conl > 0)
	{
		F_HLLC.f1 = FL.f1 * ubl;
		F_HLLC.f2 = FL.f2 * ubl + η1 * CL.p;
		F_HLLC.f3 = FL.f3 * ubl + η2 * CL.p;
		F_HLLC.f4 = FL.f4 * ubl - η3 * CL.p;
	}
	else if (conl <= 0 && ss > 0)
	{
		F_HLLC.f1 = FML.f1 * ss;
		F_HLLC.f2 = FML.f2 * ss + η1 * ps;
		F_HLLC.f3 = FML.f3 * ss + η2 * ps;
		F_HLLC.f4 = FML.f4 * ss - η3 * ps;
	}
	else if (ss <= 0 && conr >= 0)
	{
		F_HLLC.f1 = FMR.f1 * ss;
		F_HLLC.f2 = FMR.f2 * ss + η1 * ps;
		F_HLLC.f3 = FMR.f3 * ss + η2 * ps;
		F_HLLC.f4 = FMR.f4 * ss - η3 * ps;
	}
	else
	{
		F_HLLC.f1 = FR.f1 * ubr;
		F_HLLC.f2 = FR.f2 * ubr + η1 * CR.p;
		F_HLLC.f3 = FR.f3 * ubr + η2 * CR.p;
		F_HLLC.f4 = FR.f4 * ubr - η3 * CR.p;
	}
	F_HLLC.f1 = F_HLLC.f1 * Dη;
	F_HLLC.f2 = F_HLLC.f2 * Dη;
	F_HLLC.f3 = F_HLLC.f3 * Dη;
	F_HLLC.f4 = F_HLLC.f4 * Dη;

	return F_HLLC;
}
