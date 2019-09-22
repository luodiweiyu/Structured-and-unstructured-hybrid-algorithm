//#include<cmath>
//#include"const.h"
//#include<iostream>
//#include<omp.h>
//
//using namespace std;
////2维欧拉方程组离散：张德良《计算流体力学教程》p433 13.3.7
//double * HLLC_Χ(mesh CL, mesh CR, mesh C, int method)//半点左侧右侧格点，坐标变换参考点
//{
//	double xix = C.xix[method];//此处不是xix,是xix*J-1=yeta
//	double xiy = C.xiy[method];
//	double J = C.J[method];
//	double Dxi = sqrt(xix*xix + xiy * xiy);
//	double xicL = (xix * CL.u + xiy * CL.v) / Dxi;
//	double xicR = (xix * CR.u + xiy * CR.v) / Dxi;
//	double aL = sqrt(gama*CL.p / CL.ρ);
//	double aR = sqrt(gama*CR.p / CR.ρ);
//	double EL = CL.p / (gama - 1) + 0.5*CL.ρ * CL.u * CL.u + 0.5*CL.ρ * CL.v * CL.v;
//	double ER = CR.p / (gama - 1) + 0.5*CR.ρ * CR.u * CR.u + 0.5*CR.ρ * CR.v * CR.v;
//	double FL[4], FR[4];
//
//	FL[0] = (CL.ρ*xicL);
//	FL[1] = (CL.ρ*CL.u*xicL + xix * CL.p);
//	FL[2] = (CL.ρ*CL.v*xicL + xiy * CL.p);
//	FL[3] = (xicL * (EL + CL.p));
//	FR[0] = (CR.ρ*xicR);
//	FR[1] = (CR.ρ*CR.u*xicR + xix * CR.p);
//	FR[2] = (CR.ρ*CR.v*xicR + xiy * CR.p);
//	FR[3] = (xicR * (ER + CR.p));
//	double SL = xicL - aL;
//	double SR = xicR + aR;
//	double SM = (CL.p - CR.p + CL.ρ * xicL * (xicL - SL) + CR.ρ * xicR* (SR - xicR)) / (CR.ρ * (SR - xicR) - CL.ρ * (SL - xicL));
//	double pM = 0.5*(CL.ρ * (xicL - SL)*(xicL - SM) + CR.ρ * (xicR - SR)*(xicR - SM) + CR.p + CL.p);
//	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
//	double ΩL = 1 / (SL - SM);
//	double ΩR = 1 / (SR - SM);
//
//	double ρLM = ΩL * CL.ρ*(SL - xicL);
//	double ρuLS = ΩL * ((SL - xicL)*(CL.ρ*CL.u) + (pM - CL.p)*xix);
//	double ρvLS = ΩL * ((SL - xicL)*(CL.ρ*CL.v) + (pM - CL.p)*xiy);
//	double eLS = ΩL * ((SL - xicL)*EL - CL.p*xicL + pM * SM);
//
//	double ρRM = ΩR * CR.ρ*(SR - xicR);
//	double ρuRS = ΩR * ((SR - xicR)*(CR.ρ*CR.u) + (pM - CR.p)*xix);
//	double ρvRS = ΩR * ((SR - xicR)*(CR.ρ*CR.v) + (pM - CR.p)*xiy);
//	double eRS = ΩR * ((SR - xicR)*ER - CR.p*xicR + pM * SM);
//
//	double  FML[4], FMR[4];
//	FML[0] = SM * ρLM;
//	FML[1] = ρuLS * SM + pM * xix;
//	FML[2] = ρvLS * SM + pM * xiy;
//	FML[3] = (eLS + pM)*SM;
//	FMR[0] = SM * ρRM;
//	FMR[1] = ρuRS * SM + pM * xix;
//	FMR[2] = ρvRS * SM + pM * xiy;
//	FMR[3] = (eRS + pM)*SM;
//
//	double *F_HLLC = new double[4];
//	if (SL > 0)
//		for (int i = 0; i < 4; i++)
//			F_HLLC[i] = FL[i];
//	else if (SL <= 0 && SM > 0)
//		for (int i = 0; i < 4; i++)
//			F_HLLC[i] = FML[i];
//	else if (SM <= 0 && SR >= 0)
//		for (int i = 0; i < 4; i++)
//			F_HLLC[i] = FMR[i];
//	else
//		for (int i = 0; i < 4; i++)
//			F_HLLC[i] = FR[i];
//	for (int i = 0; i < 4; i++)
//		F_HLLC[i] = F_HLLC[i] *Dxi;
//	return F_HLLC;
//}
//
//double * HLLC_Υ(mesh CD, mesh CU, mesh C, int method)
//{
//	double etax = C.etax[method];
//	double etay = C.etay[method];
//	double J = C.J[method];
//	double Deta = sqrt(etax*etax + etay * etay);
//	double etacU = (C.etax[method] * CU.u + C.etay[method] * CD.v) / Deta;
//	double etacD = (C.etax[method] * CD.u + C.etay[method] * CU.v) / Deta;
//	double aU = sqrt(gama*CU.p / CU.ρ);
//	double aD = sqrt(gama*CD.p / CD.ρ);
//	double EU = CU.p / (gama - 1) + 0.5*CU.ρ * CU.u * CU.u + 0.5*CU.ρ * CU.v * CU.v;
//	double ED = CD.p / (gama - 1) + 0.5*CD.ρ * CD.u * CD.u + 0.5*CD.ρ * CD.v * CD.v;
//	double GD[4], GU[4];
//
//	GD[0] = (CD.ρ*etacD);
//	GD[1] = (CD.ρ*CD.u*etacD + etax * CD.p);
//	GD[2] = (CD.ρ*CD.v*etacD + etay * CD.p);
//	GD[3] = (etacD * (ED + CD.p));
//
//	GU[0] = (CU.ρ*etacU);
//	GU[1] = (CU.ρ*CU.u*etacU + etax * CD.p);
//	GU[2] = (CU.ρ*CU.v*etacU + etay * CD.p);
//	GU[3] = (etacU * (EU + CU.p));
//
//	double SD = etacD - aD;
//	double SU = etacU + aU;
//	double SM = (CD.p - CU.p - CD.ρ * etax * (SD - etacD) + CU.ρ * etax* (SU - etacU)) / (CU.ρ * (SU - etacU) - CD.ρ * (SD - etacD));
//	double pM = 0.5*(CD.ρ * (etacD - SD)*(etacD - SM) + CU.ρ * (etacU - SU)*(etacU - SM) + CU.p + CD.p);
//	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
//	double ΩD = 1 / (SD - SM);
//	double ΩU = 1 / (SU - SM);
//
//	double ρDM = ΩD * CD.ρ*(SD - etacD);
//	double ρuDS = ΩD * ((SD - etacD)*(CD.ρ*CD.u) + (pM - CD.p)*etax);
//	double ρvDS = ΩD * ((SD - etacD)*(CD.ρ*CD.v) + (pM - CD.p)*etay);
//	double eDS = ΩD * ((SD - etacD)*ED - CD.p*etacD + pM * SM);
//
//	double ρUM = ΩU * CU.ρ*(SU - etacU);
//	double ρuUS = ΩU * ((SU - etacU)*(CU.ρ*CU.u) + (pM - CU.p)*etax);
//	double ρvUS = ΩU * ((SU - etacU)*(CU.ρ*CU.v) + (pM - CU.p)*etay);
//	double eUS = ΩU * ((SU - etacU)*EU - CU.p*etacU + pM * SM);
//
//	double  FMD[4], FMU[4];
//
//	FMD[0] = SM * ρDM;
//	FMD[1] = ρuDS * SM + pM * etax;
//	FMD[2] = ρvDS * SM + pM * etay;
//	FMD[3] = (eDS + pM)*SM;
//	FMU[0] = SM * ρUM;
//	FMU[1] = ρuUS * SM + pM * etax;
//	FMU[2] = ρvUS * SM + pM * etay;
//	FMU[3] = (eUS + pM)*SM;
//
//	double *G_HLLC = new double[4];
//	if (SD >= 0)
//		for (int i = 0; i < 4; i++)
//			G_HLLC[i] = GD[i];
//	else if (SD <= 0 && SM >= 0)
//		for (int i = 0; i < 4; i++)
//			G_HLLC[i] = FMD[i];
//	else if (SM <= 0 && SU >= 0)
//		for (int i = 0; i < 4; i++)
//			G_HLLC[i] = FMU[i];
//	else
//		for (int i = 0; i < 4; i++)
//			G_HLLC[i] = GU[i];
//	for (int i = 0; i < 4; i++)
//		G_HLLC[i] = G_HLLC[i] *Deta;
//
//	return G_HLLC;
//}
