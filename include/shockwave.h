#pragma once
double get_c(double rho, double p);//声速公式

double get_Mas(double p1, double p2);
double get_p2p1(double Ma1);
double get_rho2rho1(double Ma1);
double get_Ma2(double Ma1);
double get_un(mesh& A, double beta);
double get_ut(mesh& A, double beta);
void get_down(mesh& U, mesh& D, double beta);//下游参数

double get_Ma(double u, double v, double rho, double p);//求马赫数
double get_Ma2(double Ma1, double beta);//斜激波波后马赫数，钱翼稷《空气动力学》p241,7-131
//double get_delta(double Ma1, double beta);//气流折射角
double get_delta(double u, double v);
double get_beta(mesh& U, double p2, int type);

double get_betafromdelta(double Ma1, double delta);
double get_ufromMa2(double Ma2, double rho2, double p2, double delta);
double get_vfromMa2(double Ma2, double rho2, double p2, double delta);
double get_p2(double Ma1, double beta,double p1);
double get_rho2(double Ma1, double beta,double rho1);
double get_Mu(mesh& U, mesh& D, double theta);
double get_udn(mesh& U, mesh& D, double Ma1, double Vs, double theta);
double getlambdafromMa(double Ma);
double getMafromlambda(double lambda);
