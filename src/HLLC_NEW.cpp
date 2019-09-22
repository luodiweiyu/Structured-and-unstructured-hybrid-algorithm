#include<cmath>
#include"const.h"
#include<iostream>
#include<omp.h>
#include"functions.h"
using namespace std;
using namespace ConstPara;
//2άŷ����������ɢ���ŵ���������������ѧ�̡̳�p433 13.3.7
Flux HLLC_��(mesh& CL, mesh& CR, mesh& C, int method)//�������Ҳ��㣬����任�ο���
{
	double xix = C.xix[method];
	double xiy = C.xiy[method];
	double xit = C.xit[method];
	double J = C.J[method];
	double Dxi = sqrt(xix * xix + xiy * xiy);
	double xi1 = xix / Dxi;
	double xi2 = xiy / Dxi;
	double xi3 = xit / Dxi;
	double xicL = CL.u * xi1 + CL.v * xi2 + xi3;
	double xicR = CR.u * xi1 + CR.v * xi2 + xi3;
	double aL = sqrt(gama * CL.p / CL.��);
	double aR = sqrt(gama * CR.p / CR.��);
	double a = sqrt(CL.��);
	double b = sqrt(CR.��);

	double SL = xicL - aL;
	double SR = xicR + aR;

	double EL = CL.p / (gama - 1) + 0.5 * CL.�� * CL.u * CL.u + 0.5 * CL.�� * CL.v * CL.v;
	double ER = CR.p / (gama - 1) + 0.5 * CR.�� * CR.u * CR.u + 0.5 * CR.�� * CR.v * CR.v;
	Flux FL, FR;
	FL.f1 = (CL.�� * xicL);
	FL.f2 = (CL.�� * CL.u * xicL + xi1 * CL.p);
	FL.f3 = (CL.�� * CL.v * xicL + xi2 * CL.p);
	FL.f4 = (xicL * (EL + CL.p) - xi3 * CL.p);
	FR.f1 = (CR.�� * xicR);
	FR.f2 = (CR.�� * CR.u * xicR + xi1 * CR.p);
	FR.f3 = (CR.�� * CR.v * xicR + xi2 * CR.p);
	FR.f4 = (xicR * (ER + CR.p) - xi3 * CR.p);
	double SM = (CL.p - CR.p + CL.�� * xicL * (xicL - SL) + CR.�� * xicR * (SR - xicR)) / (CR.�� * (SR - xicR) + CL.�� * (xicL - SL));
	double pM = 0.5 * (CL.�� * (xicL - SL) * (xicL - SM) + CR.�� * (xicR - SR) * (xicR - SM) + CR.p + CL.p);
	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
	double par1 = (xicL - SL) / (SM - SL);
	double par2 = (pM - CL.p) / (SM - SL);
	Flux  FML;
	FML.f1 = (SM * par1 * CL.��);
	FML.f2 = (SM * par1 * CL.�� * CL.u - par2 * xi1 + xi1 * pM);
	FML.f3 = (SM * par1 * CL.�� * CL.v - par2 * xi2 + xi2 * pM);
	FML.f4 = (SM * par1 * (EL + CL.p) + par2 * xi3 - xi3 * pM);
	double par3 = (xicR - SR) / (SM - SR);
	double par4 = (pM - CR.p) / (SM - SR);
	Flux FMR;
	FMR.f1 = (SM * par3 * CR.��);
	FMR.f2 = (SM * par3 * CR.�� * CR.u - par4 * xi1 + xi1 * pM);
	FMR.f3 = (SM * par3 * CR.�� * CR.v - par4 * xi2 + xi2 * pM);
	FMR.f4 = (SM * par3 * (ER + CR.p) + par4 * xi3 - xi3 * pM);

	//double ��L = 1 / (SL - SM);
	//double ��R = 1 / (SR - SM);

	//double ��LM = ��L * CL.�� * (SL - xicL);
	//double ��uLS = ��L * ((SL - xicL) * (CL.�� * CL.u) + (pM - CL.p) * xi1);
	//double ��vLS = ��L * ((SL - xicL) * (CL.�� * CL.v) + (pM - CL.p) * xi2);
	//double eLS = ��L * ((SL - xicL) * EL - CL.p * xicL + pM * SM - (pM - CL.p) * xi3);

	//double ��RM = ��R * CR.�� * (SR - xicR);
	//double ��uRS = ��R * ((SR - xicR) * (CR.�� * CR.u) + (pM - CR.p) * xi1);
	//double ��vRS = ��R * ((SR - xicR) * (CR.�� * CR.v) + (pM - CR.p) * xi2);
	//double eRS = ��R * ((SR - xicR) * ER - CR.p * xicR + pM * SM - (pM - CR.p) * xi3);

	//Flux FML, FMR;
	//FML.f1 = SM * ��LM;
	//FML.f2 = ��uLS * SM + pM * xi1;
	//FML.f3 = ��vLS * SM + pM * xi2;
	//FML.f4 = (eLS + pM) * SM - pM * xi3;
	//FMR.f1 = SM * ��RM;
	//FMR.f2 = ��uRS * SM + pM * xi1;
	//FMR.f3 = ��vRS * SM + pM * xi2;
	//FMR.f4 = (eRS + pM) * SM - pM * xi3;

	Flux F_HLLC;
	if (SL > 0)
		F_HLLC = FL;
	else if (SL <= 0 && SM > 0)
		F_HLLC = FML;
	else if (SM <= 0 && SR >= 0)
		F_HLLC = FMR;
	else
		F_HLLC = FR;
	F_HLLC.f1 = F_HLLC.f1 * Dxi;
	F_HLLC.f2 = F_HLLC.f2 * Dxi;
	F_HLLC.f3 = F_HLLC.f3 * Dxi;
	F_HLLC.f4 = F_HLLC.f4 * Dxi;

	return F_HLLC;
}

Flux HLLC_��(mesh & CD, mesh & CU, mesh & C, int method)
{
	double etax = C.etax[method];
	double etay = C.etay[method];
	double etat = C.etat[method];
	double J = C.J[method];

	double Deta = sqrt(etax * etax + etay * etay);
	double eta1 = etax / Deta;
	double eta2 = etay / Deta;
	double eta3 = etat / Deta;
	double etacU = CU.u * eta1 + CU.v * eta2;
	double etacD = CD.u * eta1 + CD.v * eta2;
	double aU = sqrt(gama * CU.p / CU.��);
	double aD = sqrt(gama * CD.p / CD.��);
	double EU = CU.p / (gama - 1) + 0.5 * CU.�� * CU.u * CU.u + 0.5 * CU.�� * CU.v * CU.v;
	double ED = CD.p / (gama - 1) + 0.5 * CD.�� * CD.u * CD.u + 0.5 * CD.�� * CD.v * CD.v;
	Flux GD, GU;
	GD.f1 = (CD.�� * etacD);
	GD.f2 = (CD.�� * CD.u * etacD + eta1 * CD.p);
	GD.f3 = (CD.�� * CD.v * etacD + eta2 * CD.p);
	GD.f4 = (etacD * (ED + CD.p) - eta3 * CD.p);

	GU.f1 = (CU.�� * etacU);
	GU.f2 = (CU.�� * CU.u * etacU + eta1 * CD.p);
	GU.f3 = (CU.�� * CU.v * etacU + eta2 * CD.p);
	GU.f4 = (etacU * (EU + CU.p) - eta3 * CD.p);

	double SD = etacD - aD;
	double SU = etacU + aU;
	double SM = (CD.p - CU.p - CD.�� * etacD * (SD - etacD) + CU.�� * etacU * (SU - etacU)) / (CU.�� * (SU - etacU) - CD.�� * (SD - etacD));
	double pM = 0.5 * (CD.�� * (etacD - SD) * (etacD - SM) + CU.�� * (etacU - SU) * (etacU - SM) + CU.p + CD.p);
	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
	double par1 = (etacD - SD) / (SM - SD);
	double par2 = (pM - CD.p) / (SM - SD);
	Flux  FMD;
	FMD.f1 = (SM * par1 * CD.��);
	FMD.f2 = (SM * par1 * CD.�� * CD.u - par2 * eta1 + eta1 * pM);
	FMD.f3 = (SM * par1 * CD.�� * CD.v - par2 * eta2 + eta2 * pM);
	FMD.f4 = (SM * par1 * (ED + CD.p) + par2 * eta3 - eta3 * pM);
	double par3 = (etacU - SU) / (SM - SU);
	double par4 = (pM - CU.p) / (SM - SU);
	Flux FMU;
	FMU.f1 = (SM * par3 * CU.��);
	FMU.f2 = (SM * par3 * CU.�� * CU.u - par4 * eta1 + eta1 * pM);
	FMU.f3 = (SM * par3 * CU.�� * CU.v - par4 * eta2 + eta2 * pM);
	FMU.f4 = (SM * par3 * (EU + CU.p) + par4 * eta3 - eta3 * pM);

	//double ��D = 1 / (SD - SM);
	//double ��U = 1 / (SU - SM);

	//double ��DM = ��D * CD.�� * (SD - etacD);
	//double ��uDS = ��D * ((SD - etacD) * (CD.�� * CD.u) + (pM - CD.p) * eta1);
	//double ��vDS = ��D * ((SD - etacD) * (CD.�� * CD.v) + (pM - CD.p) * eta2);
	//double eDS = ��D * ((SD - etacD) * ED - CD.p * etacD + pM * SM - (pM - CD.p) * eta3);

	//double ��UM = ��U * CU.�� * (SU - etacU);
	//double ��uUS = ��U * ((SU - etacU) * (CU.�� * CU.u) + (pM - CU.p) * eta1);
	//double ��vUS = ��U * ((SU - etacU) * (CU.�� * CU.v) + (pM - CU.p) * eta2);
	//double eUS = ��U * ((SU - etacU) * EU - CU.p * etacU + pM * SM - (pM - CU.p) * eta3);

	//Flux  FMD, FMU;

	//FMD.f1 = SM * ��DM;
	//FMD.f2 = ��uDS * SM + pM * eta1;
	//FMD.f3 = ��vDS * SM + pM * eta2;
	//FMD.f4 = (eDS + pM) * SM - pM * eta3;
	//FMU.f1 = SM * ��UM;
	//FMU.f2 = ��uUS * SM + pM * eta1;
	//FMU.f3 = ��vUS * SM + pM * eta2;
	//FMU.f4 = (eUS + pM) * SM - pM * eta3;

	Flux G_HLLC;
	if (SD >= 0)
		G_HLLC = GD;
	else if (SD <= 0 && SM >= 0)
		G_HLLC = FMD;
	else if (SM <= 0 && SU >= 0)
		G_HLLC = FMU;
	else
		G_HLLC = GU;
	G_HLLC.f1 = G_HLLC.f1 * Deta;
	G_HLLC.f2 = G_HLLC.f2 * Deta;
	G_HLLC.f3 = G_HLLC.f3 * Deta;
	G_HLLC.f4 = G_HLLC.f4 * Deta;

	return G_HLLC;
}

Flux HLLC_��2(mesh& CL, mesh& CR, mesh& C, int method)//�������Ҳ��㣬����任�ο���
{
	double xix = C.xix[method];
	double xiy = C.xiy[method];
	double xit = C.xit[method];
	double J = C.J[method];
	double Dxi = sqrt(xix * xix + xiy * xiy);
	double xi1 = xix / Dxi;
	double xi2 = xiy / Dxi;
	double xi3 = xit / Dxi;
	double aL = sqrt(gama * CL.p / CL.��);
	double aR = sqrt(gama * CR.p / CR.��);
	double a = sqrt(CL.��);
	double b = sqrt(CR.��);
	double ubl = CL.u * xi1 + CL.v * xi2 + xi3;
	double ubr = CR.u * xi1 + CR.v * xi2 + xi3;
	double ubp = (a * ubl + b * ubr) / (a + b);
	double aM = (a * aL * aL + b * aR * aR) / (a + b) + 0.5 * (a * b / ((a + b) * (a + b))) * ((CL.u - CR.u) * (CL.u - CR.u) + (CL.v - CR.v) * (CL.v - CR.v));


	double conl = min(ubl - aL, ubp - aM);
	double conr = max(ubr + aR, ubp + aM);
	double s1 = CR.�� * (conr - ubr);
	double s2 = CL.�� * (ubl - conl);
	double ss = (s1 * ubr + s2 * ubl + CL.p - CR.p) / (s1 + s2);
	double ps = s2 * (ubl - ss) + CL.p;
	Flux FL, FR;
	FL.f1 = CL.��;
	FL.f2 = CL.�� * CL.u;
	FL.f3 = CL.�� * CL.v;
	FL.f4 = CL.p * gama / (gama - 1) + 0.5 * CL.�� * CL.u * CL.u + 0.5 * CL.�� * CL.v * CL.v; 
	FR.f1 = CR.��;
	FR.f2 = CR.�� * CR.u;
	FR.f3 = CR.�� * CR.v;
	FR.f4 = CR.p * gama / (gama - 1) + 0.5 * CR.�� * CR.u * CR.u + 0.5 * CR.�� * CR.v * CR.v;
	double par1 = (ubl - conl) / (ss - conl);
	double par2 = (ps - CL.p) / (ss - conl);

	Flux FML, FMR;
	FML.f1 = par1 * FL.f1;
	FML.f2 = par1 * FL.f2 - par2 * xi1;
	FML.f3 = par1 * FL.f3 - par2 * xi2;
	FML.f4 = par1 * FL.f4 + par2 * xi3;
	double par3 = (ubr - conr) / (ss - conr);
	double par4 = (ps - CR.p) / (ss - conr);

	FMR.f1 = par3 * FR.f1;
	FMR.f2 = par3 * FR.f2 - par4 * xi1;
	FMR.f3 = par3 * FR.f3 - par4 * xi2;
	FMR.f4 = par3 * FR.f4 + par4 * xi3;

	Flux F_HLLC;
	if (conl > 0)
	{
		F_HLLC.f1 = FL.f1 * ubl;
		F_HLLC.f2 = FL.f2 * ubl + xi1 * CL.p;
		F_HLLC.f3 = FL.f3 * ubl + xi2 * CL.p;
		F_HLLC.f4 = FL.f4 * ubl - xi3 * CL.p;
	}
	else if (conl <= 0 && ss > 0)
	{
		F_HLLC.f1 = FML.f1 * ss;
		F_HLLC.f2 = FML.f2 * ss + xi1 * ps;
		F_HLLC.f3 = FML.f3 * ss + xi2 * ps;
		F_HLLC.f4 = FML.f4 * ss - xi3 * ps;
	}
	else if (ss <= 0 && conr >= 0)
	{
		F_HLLC.f1 = FMR.f1 * ss;
		F_HLLC.f2 = FMR.f2 * ss + xi1 * ps;
		F_HLLC.f3 = FMR.f3 * ss + xi2 * ps;
		F_HLLC.f4 = FMR.f4 * ss - xi3 * ps;
	}
	else
	{
		F_HLLC.f1 = FR.f1 * ubr;
		F_HLLC.f2 = FR.f2 * ubr + xi1 * CR.p;
		F_HLLC.f3 = FR.f3 * ubr + xi2 * CR.p;
		F_HLLC.f4 = FR.f4 * ubr - xi3 * CR.p;
	}
	F_HLLC.f1 = F_HLLC.f1 * Dxi;
	F_HLLC.f2 = F_HLLC.f2 * Dxi;
	F_HLLC.f3 = F_HLLC.f3 * Dxi;
	F_HLLC.f4 = F_HLLC.f4 * Dxi;

	return F_HLLC;
}
Flux HLLC_��2(mesh& CL, mesh& CR, mesh& C, int method)
{
	double etax = C.etax[method];
	double etay = C.etay[method];
	double etat = C.etat[method];
	double J = C.J[method];
	double Deta = sqrt(etax * etax + etay * etay);
	double eta1 = etax / Deta;
	double eta2 = etay / Deta;
	double eta3 = etat / Deta;
	double aL = sqrt(gama * CL.p / CL.��);
	double aR = sqrt(gama * CR.p / CR.��);
	double a = sqrt(CL.��);
	double b = sqrt(CR.��);
	double ubl = CL.u * eta1 + CL.v * eta2 + eta3;
	double ubr = CR.u * eta1 + CR.v * eta2 + eta3;
	double ubp = (a * ubl + b * ubr) / (a + b);
	double aM = (a * aL * aL + b * aR * aR) / (a + b) + 0.5 * (a * b / ((a + b) * (a + b))) * ((CL.u - CR.u) * (CL.u - CR.u) + (CL.v - CR.v) * (CL.v - CR.v));


	double conl = min(ubl - aL, ubp - aM);
	double conr = max(ubr + aR, ubp + aM);
	double s1 = CR.�� * (conr - ubr);
	double s2 = CL.�� * (ubl - conl);
	double ss = (s1 * ubr + s2 * ubl + CL.p - CR.p) / (s1 + s2);
	double ps = s2 * (ubl - ss) + CL.p;
	Flux FL, FR;
	FL.f1 = CL.��;
	FL.f2 = CL.�� * CL.u;
	FL.f3 = CL.�� * CL.v;
	FL.f4 = CL.p * gama / (gama - 1) + 0.5 * CL.�� * CL.u * CL.u + 0.5 * CL.�� * CL.v * CL.v; 
	FR.f1 = CR.��;
	FR.f2 = CR.�� * CR.u;
	FR.f3 = CR.�� * CR.v;
	FR.f4 = CR.p * gama / (gama - 1) + 0.5 * CR.�� * CR.u * CR.u + 0.5 * CR.�� * CR.v * CR.v;
	double par1 = (ubl - conl) / (ss - conl);
	double par2 = (ps - CL.p) / (ss - conl);

	Flux FML, FMR;
	FML.f1 = par1 * FL.f1;
	FML.f2 = par1 * FL.f2 - par2 * eta1;
	FML.f3 = par1 * FL.f3 - par2 * eta2;
	FML.f4 = par1 * FL.f4 + par2 * eta3;
	double par3 = (ubr - conr) / (ss - conr);
	double par4 = (ps - CR.p) / (ss - conr);

	FMR.f1 = par3 * FR.f1;
	FMR.f2 = par3 * FR.f2 - par4 * eta1;
	FMR.f3 = par3 * FR.f3 - par4 * eta2;
	FMR.f4 = par3 * FR.f4 + par4 * eta3;

	Flux F_HLLC;
	if (conl > 0)
	{
		F_HLLC.f1 = FL.f1 * ubl;
		F_HLLC.f2 = FL.f2 * ubl + eta1 * CL.p;
		F_HLLC.f3 = FL.f3 * ubl + eta2 * CL.p;
		F_HLLC.f4 = FL.f4 * ubl - eta3 * CL.p;
	}
	else if (conl <= 0 && ss > 0)
	{
		F_HLLC.f1 = FML.f1 * ss;
		F_HLLC.f2 = FML.f2 * ss + eta1 * ps;
		F_HLLC.f3 = FML.f3 * ss + eta2 * ps;
		F_HLLC.f4 = FML.f4 * ss - eta3 * ps;
	}
	else if (ss <= 0 && conr >= 0)
	{
		F_HLLC.f1 = FMR.f1 * ss;
		F_HLLC.f2 = FMR.f2 * ss + eta1 * ps;
		F_HLLC.f3 = FMR.f3 * ss + eta2 * ps;
		F_HLLC.f4 = FMR.f4 * ss - eta3 * ps;
	}
	else
	{
		F_HLLC.f1 = FR.f1 * ubr;
		F_HLLC.f2 = FR.f2 * ubr + eta1 * CR.p;
		F_HLLC.f3 = FR.f3 * ubr + eta2 * CR.p;
		F_HLLC.f4 = FR.f4 * ubr - eta3 * CR.p;
	}
	F_HLLC.f1 = F_HLLC.f1 * Deta;
	F_HLLC.f2 = F_HLLC.f2 * Deta;
	F_HLLC.f3 = F_HLLC.f3 * Deta;
	F_HLLC.f4 = F_HLLC.f4 * Deta;

	return F_HLLC;
}
