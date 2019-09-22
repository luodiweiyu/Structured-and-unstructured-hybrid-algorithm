#pragma once
double get_c(double ρ, double p);//声速公式

double get_Mas(double p1, double p2);
double get_p2p1(double Ma1);
double get_ρ2ρ1(double Ma1);
double get_Ma2(double Ma1);
double get_un(mesh& A, double β);
double get_ut(mesh& A, double β);
void get_down(mesh& U, mesh& D, double β);//下游参数

double get_Ma(double u, double v, double ρ, double p);//求马赫数
double get_Ma2(double Ma1, double β);//斜激波波后马赫数，钱翼稷《空气动力学》p241,7-131
//double get_δ(double Ma1, double β);//气流折射角
double get_δ(double u, double v);
double get_β(mesh& U, double p2, int type);

double get_βfromδ(double Ma1, double δ);
double get_ufromMa2(double Ma2, double ρ2, double p2, double δ);
double get_vfromMa2(double Ma2, double ρ2, double p2, double δ);
double get_p2(double Ma1, double β,double p1);
double get_ρ2(double Ma1, double β,double ρ1);
double get_Mu(mesh& U, mesh& D, double θ);
double get_udn(mesh& U, mesh& D, double Ma1, double Vs, double θ);
double getλfromMa(double Ma);
double getMafromλ(double λ);
