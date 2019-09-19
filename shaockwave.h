#pragma once
double get_Ma(double u, double v, double ρ, double p);//求马赫数
double get_Ma2(double Ma1, double β);//斜激波波后马赫数，钱翼稷《空气动力学》p241,7-131
double get_δ(double Ma1, double β);//气流折射角
double get_ufromMa2(double Ma2, double ρ2, double p2, double δ);
double get_vfromMa2(double Ma2, double ρ2, double p2, double δ);
double get_p2(double Ma1, double β, double p1);
double get_ρ2(double Ma1, double β, double ρ1);
