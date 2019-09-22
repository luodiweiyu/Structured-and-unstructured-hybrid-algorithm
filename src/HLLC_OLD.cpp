//#include<cmath>
//#include"const.h"
//#include<iostream>
//#include<omp.h>
//
//using namespace std;
////2άŷ����������ɢ���ŵ���������������ѧ�̡̳�p433 13.3.7
//double * HLLC_��(mesh CL, mesh CR, mesh C, int method)//�������Ҳ��㣬����任�ο���
//{
//	double xix = C.xix[method];//�˴�����xix,��xix*J-1=yeta
//	double xiy = C.xiy[method];
//	double J = C.J[method];
//	double Dxi = sqrt(xix*xix + xiy * xiy);
//	double xicL = (xix * CL.u + xiy * CL.v) / Dxi;
//	double xicR = (xix * CR.u + xiy * CR.v) / Dxi;
//	double aL = sqrt(gama*CL.p / CL.��);
//	double aR = sqrt(gama*CR.p / CR.��);
//	double EL = CL.p / (gama - 1) + 0.5*CL.�� * CL.u * CL.u + 0.5*CL.�� * CL.v * CL.v;
//	double ER = CR.p / (gama - 1) + 0.5*CR.�� * CR.u * CR.u + 0.5*CR.�� * CR.v * CR.v;
//	double FL[4], FR[4];
//
//	FL[0] = (CL.��*xicL);
//	FL[1] = (CL.��*CL.u*xicL + xix * CL.p);
//	FL[2] = (CL.��*CL.v*xicL + xiy * CL.p);
//	FL[3] = (xicL * (EL + CL.p));
//	FR[0] = (CR.��*xicR);
//	FR[1] = (CR.��*CR.u*xicR + xix * CR.p);
//	FR[2] = (CR.��*CR.v*xicR + xiy * CR.p);
//	FR[3] = (xicR * (ER + CR.p));
//	double SL = xicL - aL;
//	double SR = xicR + aR;
//	double SM = (CL.p - CR.p + CL.�� * xicL * (xicL - SL) + CR.�� * xicR* (SR - xicR)) / (CR.�� * (SR - xicR) - CL.�� * (SL - xicL));
//	double pM = 0.5*(CL.�� * (xicL - SL)*(xicL - SM) + CR.�� * (xicR - SR)*(xicR - SM) + CR.p + CL.p);
//	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
//	double ��L = 1 / (SL - SM);
//	double ��R = 1 / (SR - SM);
//
//	double ��LM = ��L * CL.��*(SL - xicL);
//	double ��uLS = ��L * ((SL - xicL)*(CL.��*CL.u) + (pM - CL.p)*xix);
//	double ��vLS = ��L * ((SL - xicL)*(CL.��*CL.v) + (pM - CL.p)*xiy);
//	double eLS = ��L * ((SL - xicL)*EL - CL.p*xicL + pM * SM);
//
//	double ��RM = ��R * CR.��*(SR - xicR);
//	double ��uRS = ��R * ((SR - xicR)*(CR.��*CR.u) + (pM - CR.p)*xix);
//	double ��vRS = ��R * ((SR - xicR)*(CR.��*CR.v) + (pM - CR.p)*xiy);
//	double eRS = ��R * ((SR - xicR)*ER - CR.p*xicR + pM * SM);
//
//	double  FML[4], FMR[4];
//	FML[0] = SM * ��LM;
//	FML[1] = ��uLS * SM + pM * xix;
//	FML[2] = ��vLS * SM + pM * xiy;
//	FML[3] = (eLS + pM)*SM;
//	FMR[0] = SM * ��RM;
//	FMR[1] = ��uRS * SM + pM * xix;
//	FMR[2] = ��vRS * SM + pM * xiy;
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
//double * HLLC_��(mesh CD, mesh CU, mesh C, int method)
//{
//	double etax = C.etax[method];
//	double etay = C.etay[method];
//	double J = C.J[method];
//	double Deta = sqrt(etax*etax + etay * etay);
//	double etacU = (C.etax[method] * CU.u + C.etay[method] * CD.v) / Deta;
//	double etacD = (C.etax[method] * CD.u + C.etay[method] * CU.v) / Deta;
//	double aU = sqrt(gama*CU.p / CU.��);
//	double aD = sqrt(gama*CD.p / CD.��);
//	double EU = CU.p / (gama - 1) + 0.5*CU.�� * CU.u * CU.u + 0.5*CU.�� * CU.v * CU.v;
//	double ED = CD.p / (gama - 1) + 0.5*CD.�� * CD.u * CD.u + 0.5*CD.�� * CD.v * CD.v;
//	double GD[4], GU[4];
//
//	GD[0] = (CD.��*etacD);
//	GD[1] = (CD.��*CD.u*etacD + etax * CD.p);
//	GD[2] = (CD.��*CD.v*etacD + etay * CD.p);
//	GD[3] = (etacD * (ED + CD.p));
//
//	GU[0] = (CU.��*etacU);
//	GU[1] = (CU.��*CU.u*etacU + etax * CD.p);
//	GU[2] = (CU.��*CU.v*etacU + etay * CD.p);
//	GU[3] = (etacU * (EU + CU.p));
//
//	double SD = etacD - aD;
//	double SU = etacU + aU;
//	double SM = (CD.p - CU.p - CD.�� * etax * (SD - etacD) + CU.�� * etax* (SU - etacU)) / (CU.�� * (SU - etacU) - CD.�� * (SD - etacD));
//	double pM = 0.5*(CD.�� * (etacD - SD)*(etacD - SM) + CU.�� * (etacU - SU)*(etacU - SM) + CU.p + CD.p);
//	//Average-State Jacobians and Implicit Methods for Compressible Viscous and Turbulent Flows,(7)(8)
//	double ��D = 1 / (SD - SM);
//	double ��U = 1 / (SU - SM);
//
//	double ��DM = ��D * CD.��*(SD - etacD);
//	double ��uDS = ��D * ((SD - etacD)*(CD.��*CD.u) + (pM - CD.p)*etax);
//	double ��vDS = ��D * ((SD - etacD)*(CD.��*CD.v) + (pM - CD.p)*etay);
//	double eDS = ��D * ((SD - etacD)*ED - CD.p*etacD + pM * SM);
//
//	double ��UM = ��U * CU.��*(SU - etacU);
//	double ��uUS = ��U * ((SU - etacU)*(CU.��*CU.u) + (pM - CU.p)*etax);
//	double ��vUS = ��U * ((SU - etacU)*(CU.��*CU.v) + (pM - CU.p)*etay);
//	double eUS = ��U * ((SU - etacU)*EU - CU.p*etacU + pM * SM);
//
//	double  FMD[4], FMU[4];
//
//	FMD[0] = SM * ��DM;
//	FMD[1] = ��uDS * SM + pM * etax;
//	FMD[2] = ��vDS * SM + pM * etay;
//	FMD[3] = (eDS + pM)*SM;
//	FMU[0] = SM * ��UM;
//	FMU[1] = ��uUS * SM + pM * etax;
//	FMU[2] = ��vUS * SM + pM * etay;
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
