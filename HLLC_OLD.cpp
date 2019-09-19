//#include<cmath>
//#include"const.h"
//#include<iostream>
//#include<omp.h>
//
//using namespace std;
////2维欧拉方程组离散：张德良《计算流体力学教程》p433 13.3.7
//double * HLLC_Χ(mesh CL, mesh CR, mesh C, int method)//半点左侧右侧格点，坐标变换参考点
//{
//	double ξx = C.ξx[method];//此处不是ξx,是ξx*J-1=yη
//	double ξy = C.ξy[method];
//	double J = C.J[method];
//	double Dξ = sqrt(ξx*ξx + ξy * ξy);
//	double ξcL = (ξx * CL.u + ξy * CL.v) / Dξ;
//	double ξcR = (ξx * CR.u + ξy * CR.v) / Dξ;
//	double aL = sqrt(γ*CL.p / CL.ρ);
//	double aR = sqrt(γ*CR.p / CR.ρ);
//	double EL = CL.p / (γ - 1) + 0.5*CL.ρ * CL.u * CL.u + 0.5*CL.ρ * CL.v * CL.v;
//	double ER = CR.p / (γ - 1) + 0.5*CR.ρ * CR.u * CR.u + 0.5*CR.ρ * CR.v * CR.v;
//	double FL[4], FR[4];
//
//	FL[0] = (CL.ρ*ξcL);
//	FL[1] = (CL.ρ*CL.u*ξcL + ξx * CL.p);
//	FL[2] = (CL.ρ*CL.v*ξcL + ξy * CL.p);
//	FL[3] = (ξcL * (EL + CL.p));
//	FR[0] = (CR.ρ*ξcR);
//	FR[1] = (CR.ρ*CR.u*ξcR + ξx * CR.p);
//	FR[2] = (CR.ρ*CR.v*ξcR + ξy * CR.p);
//	FR[3] = (ξcR * (ER + CR.p));
//	double SL = ξcL - aL;
//	double SR = ξcR + aR;
//	double SM = (CL.p - CR.p + CL.ρ * ξcL * (ξcL - SL) + CR.ρ * ξcR* (SR - ξcR)) / (CR.ρ * (SR - ξcR) - CL.ρ * (SL - ξcL));
//	double pM = 0.5*(CL.ρ * (ξcL - SL)*(ξcL - SM) + CR.ρ * (ξcR - SR)*(ξcR - SM) + CR.p + CL.p);
//	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
//	double ΩL = 1 / (SL - SM);
//	double ΩR = 1 / (SR - SM);
//
//	double ρLM = ΩL * CL.ρ*(SL - ξcL);
//	double ρuLS = ΩL * ((SL - ξcL)*(CL.ρ*CL.u) + (pM - CL.p)*ξx);
//	double ρvLS = ΩL * ((SL - ξcL)*(CL.ρ*CL.v) + (pM - CL.p)*ξy);
//	double eLS = ΩL * ((SL - ξcL)*EL - CL.p*ξcL + pM * SM);
//
//	double ρRM = ΩR * CR.ρ*(SR - ξcR);
//	double ρuRS = ΩR * ((SR - ξcR)*(CR.ρ*CR.u) + (pM - CR.p)*ξx);
//	double ρvRS = ΩR * ((SR - ξcR)*(CR.ρ*CR.v) + (pM - CR.p)*ξy);
//	double eRS = ΩR * ((SR - ξcR)*ER - CR.p*ξcR + pM * SM);
//
//	double  FML[4], FMR[4];
//	FML[0] = SM * ρLM;
//	FML[1] = ρuLS * SM + pM * ξx;
//	FML[2] = ρvLS * SM + pM * ξy;
//	FML[3] = (eLS + pM)*SM;
//	FMR[0] = SM * ρRM;
//	FMR[1] = ρuRS * SM + pM * ξx;
//	FMR[2] = ρvRS * SM + pM * ξy;
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
//		F_HLLC[i] = F_HLLC[i] *Dξ;
//	return F_HLLC;
//}
//
//double * HLLC_Υ(mesh CD, mesh CU, mesh C, int method)
//{
//	double ηx = C.ηx[method];
//	double ηy = C.ηy[method];
//	double J = C.J[method];
//	double Dη = sqrt(ηx*ηx + ηy * ηy);
//	double ηcU = (C.ηx[method] * CU.u + C.ηy[method] * CD.v) / Dη;
//	double ηcD = (C.ηx[method] * CD.u + C.ηy[method] * CU.v) / Dη;
//	double aU = sqrt(γ*CU.p / CU.ρ);
//	double aD = sqrt(γ*CD.p / CD.ρ);
//	double EU = CU.p / (γ - 1) + 0.5*CU.ρ * CU.u * CU.u + 0.5*CU.ρ * CU.v * CU.v;
//	double ED = CD.p / (γ - 1) + 0.5*CD.ρ * CD.u * CD.u + 0.5*CD.ρ * CD.v * CD.v;
//	double GD[4], GU[4];
//
//	GD[0] = (CD.ρ*ηcD);
//	GD[1] = (CD.ρ*CD.u*ηcD + ηx * CD.p);
//	GD[2] = (CD.ρ*CD.v*ηcD + ηy * CD.p);
//	GD[3] = (ηcD * (ED + CD.p));
//
//	GU[0] = (CU.ρ*ηcU);
//	GU[1] = (CU.ρ*CU.u*ηcU + ηx * CD.p);
//	GU[2] = (CU.ρ*CU.v*ηcU + ηy * CD.p);
//	GU[3] = (ηcU * (EU + CU.p));
//
//	double SD = ηcD - aD;
//	double SU = ηcU + aU;
//	double SM = (CD.p - CU.p - CD.ρ * ηx * (SD - ηcD) + CU.ρ * ηx* (SU - ηcU)) / (CU.ρ * (SU - ηcU) - CD.ρ * (SD - ηcD));
//	double pM = 0.5*(CD.ρ * (ηcD - SD)*(ηcD - SM) + CU.ρ * (ηcU - SU)*(ηcU - SM) + CU.p + CD.p);
//	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
//	double ΩD = 1 / (SD - SM);
//	double ΩU = 1 / (SU - SM);
//
//	double ρDM = ΩD * CD.ρ*(SD - ηcD);
//	double ρuDS = ΩD * ((SD - ηcD)*(CD.ρ*CD.u) + (pM - CD.p)*ηx);
//	double ρvDS = ΩD * ((SD - ηcD)*(CD.ρ*CD.v) + (pM - CD.p)*ηy);
//	double eDS = ΩD * ((SD - ηcD)*ED - CD.p*ηcD + pM * SM);
//
//	double ρUM = ΩU * CU.ρ*(SU - ηcU);
//	double ρuUS = ΩU * ((SU - ηcU)*(CU.ρ*CU.u) + (pM - CU.p)*ηx);
//	double ρvUS = ΩU * ((SU - ηcU)*(CU.ρ*CU.v) + (pM - CU.p)*ηy);
//	double eUS = ΩU * ((SU - ηcU)*EU - CU.p*ηcU + pM * SM);
//
//	double  FMD[4], FMU[4];
//
//	FMD[0] = SM * ρDM;
//	FMD[1] = ρuDS * SM + pM * ηx;
//	FMD[2] = ρvDS * SM + pM * ηy;
//	FMD[3] = (eDS + pM)*SM;
//	FMU[0] = SM * ρUM;
//	FMU[1] = ρuUS * SM + pM * ηx;
//	FMU[2] = ρvUS * SM + pM * ηy;
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
//		G_HLLC[i] = G_HLLC[i] *Dη;
//
//	return G_HLLC;
//}
