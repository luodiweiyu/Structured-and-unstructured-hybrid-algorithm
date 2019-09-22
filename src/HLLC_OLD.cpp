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
//	double aL = sqrt(gama*CL.p / CL.rho);
//	double aR = sqrt(gama*CR.p / CR.rho);
//	double EL = CL.p / (gama - 1) + 0.5*CL.rho * CL.u * CL.u + 0.5*CL.rho * CL.v * CL.v;
//	double ER = CR.p / (gama - 1) + 0.5*CR.rho * CR.u * CR.u + 0.5*CR.rho * CR.v * CR.v;
//	double FL[4], FR[4];
//
//	FL[0] = (CL.rho*xicL);
//	FL[1] = (CL.rho*CL.u*xicL + xix * CL.p);
//	FL[2] = (CL.rho*CL.v*xicL + xiy * CL.p);
//	FL[3] = (xicL * (EL + CL.p));
//	FR[0] = (CR.rho*xicR);
//	FR[1] = (CR.rho*CR.u*xicR + xix * CR.p);
//	FR[2] = (CR.rho*CR.v*xicR + xiy * CR.p);
//	FR[3] = (xicR * (ER + CR.p));
//	double SL = xicL - aL;
//	double SR = xicR + aR;
//	double SM = (CL.p - CR.p + CL.rho * xicL * (xicL - SL) + CR.rho * xicR* (SR - xicR)) / (CR.rho * (SR - xicR) - CL.rho * (SL - xicL));
//	double pM = 0.5*(CL.rho * (xicL - SL)*(xicL - SM) + CR.rho * (xicR - SR)*(xicR - SM) + CR.p + CL.p);
//	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
//	double ΩL = 1 / (SL - SM);
//	double ΩR = 1 / (SR - SM);
//
//	double rhoLM = ΩL * CL.rho*(SL - xicL);
//	double rhouLS = ΩL * ((SL - xicL)*(CL.rho*CL.u) + (pM - CL.p)*xix);
//	double rhovLS = ΩL * ((SL - xicL)*(CL.rho*CL.v) + (pM - CL.p)*xiy);
//	double eLS = ΩL * ((SL - xicL)*EL - CL.p*xicL + pM * SM);
//
//	double rhoRM = ΩR * CR.rho*(SR - xicR);
//	double rhouRS = ΩR * ((SR - xicR)*(CR.rho*CR.u) + (pM - CR.p)*xix);
//	double rhovRS = ΩR * ((SR - xicR)*(CR.rho*CR.v) + (pM - CR.p)*xiy);
//	double eRS = ΩR * ((SR - xicR)*ER - CR.p*xicR + pM * SM);
//
//	double  FML[4], FMR[4];
//	FML[0] = SM * rhoLM;
//	FML[1] = rhouLS * SM + pM * xix;
//	FML[2] = rhovLS * SM + pM * xiy;
//	FML[3] = (eLS + pM)*SM;
//	FMR[0] = SM * rhoRM;
//	FMR[1] = rhouRS * SM + pM * xix;
//	FMR[2] = rhovRS * SM + pM * xiy;
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
//	double aU = sqrt(gama*CU.p / CU.rho);
//	double aD = sqrt(gama*CD.p / CD.rho);
//	double EU = CU.p / (gama - 1) + 0.5*CU.rho * CU.u * CU.u + 0.5*CU.rho * CU.v * CU.v;
//	double ED = CD.p / (gama - 1) + 0.5*CD.rho * CD.u * CD.u + 0.5*CD.rho * CD.v * CD.v;
//	double GD[4], GU[4];
//
//	GD[0] = (CD.rho*etacD);
//	GD[1] = (CD.rho*CD.u*etacD + etax * CD.p);
//	GD[2] = (CD.rho*CD.v*etacD + etay * CD.p);
//	GD[3] = (etacD * (ED + CD.p));
//
//	GU[0] = (CU.rho*etacU);
//	GU[1] = (CU.rho*CU.u*etacU + etax * CD.p);
//	GU[2] = (CU.rho*CU.v*etacU + etay * CD.p);
//	GU[3] = (etacU * (EU + CU.p));
//
//	double SD = etacD - aD;
//	double SU = etacU + aU;
//	double SM = (CD.p - CU.p - CD.rho * etax * (SD - etacD) + CU.rho * etax* (SU - etacU)) / (CU.rho * (SU - etacU) - CD.rho * (SD - etacD));
//	double pM = 0.5*(CD.rho * (etacD - SD)*(etacD - SM) + CU.rho * (etacU - SU)*(etacU - SM) + CU.p + CD.p);
//	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
//	double ΩD = 1 / (SD - SM);
//	double ΩU = 1 / (SU - SM);
//
//	double rhoDM = ΩD * CD.rho*(SD - etacD);
//	double rhouDS = ΩD * ((SD - etacD)*(CD.rho*CD.u) + (pM - CD.p)*etax);
//	double rhovDS = ΩD * ((SD - etacD)*(CD.rho*CD.v) + (pM - CD.p)*etay);
//	double eDS = ΩD * ((SD - etacD)*ED - CD.p*etacD + pM * SM);
//
//	double rhoUM = ΩU * CU.rho*(SU - etacU);
//	double rhouUS = ΩU * ((SU - etacU)*(CU.rho*CU.u) + (pM - CU.p)*etax);
//	double rhovUS = ΩU * ((SU - etacU)*(CU.rho*CU.v) + (pM - CU.p)*etay);
//	double eUS = ΩU * ((SU - etacU)*EU - CU.p*etacU + pM * SM);
//
//	double  FMD[4], FMU[4];
//
//	FMD[0] = SM * rhoDM;
//	FMD[1] = rhouDS * SM + pM * etax;
//	FMD[2] = rhovDS * SM + pM * etay;
//	FMD[3] = (eDS + pM)*SM;
//	FMU[0] = SM * rhoUM;
//	FMU[1] = rhouUS * SM + pM * etax;
//	FMU[2] = rhovUS * SM + pM * etay;
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
