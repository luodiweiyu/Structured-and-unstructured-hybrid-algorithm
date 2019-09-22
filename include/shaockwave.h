#pragma once
double get_Ma(double u, double v, double rho, double p);//求马赫数
double get_Ma2(double Ma1, double beta);//斜激波波后马赫数，钱翼稷《空气动力学》p241,7-131
double get_delta(double Ma1, double beta);//气流折射角
double get_ufromMa2(double Ma2, double rho2, double p2, double delta);
double get_vfromMa2(double Ma2, double rho2, double p2, double delta);
double get_p2(double Ma1, double beta, double p1);
double get_rho2(double Ma1, double beta, double rho1);
