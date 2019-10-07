#pragma once
#include"const.h"
double distance(mesh& a, mesh& b);
double distance(double x1, double y1, double x2, double y2);
void get_dt();
//void update_AfromU();
//void choose_U(int i);
//void get_F();
//void get_G();
void coordinate_trans();
Flux HLLC_Χ(mesh& CL, mesh& CR, mesh& C, int method);
Flux HLLC_Χ2(mesh& CL, mesh& CR, mesh& C, int method);//半点左侧右侧格点，坐标变换参考点

Flux HLLC_Υ(mesh& CD, mesh& CU, mesh& C, int method);
Flux HLLC_Υ2(mesh& CD, mesh& CU, mesh& C, int method);

Flux VanLeerA(mesh& C, double xix, double xiy, double xit, double J);
Flux VanLeerB(mesh& C, double xix, double xiy, double xit, double J);

void record();
double get_beta(mesh A, mesh B);//求出两个网格点与x轴的夹角
void reorder_neighbor();
double area(mesh A, mesh B, mesh C, mesh D);//求任意四点构成四边形面积 
void movemesh();
void findNeiborSec();
void update_Vm();
void clear_Vm();
void update_bound();
double get_theta(double x1, double y1, double x2, double y2);//求直线与x轴的夹角

double max(double a, double b);
double min(double a, double b);
double absmax(double a, double b);
double absmin(double a, double b);
mesh getCrossPoint(double theta, double a, double b, double r);
void polygonPoint(vector <mesh> &poly);

mesh getCrossPoint(Line L1, Line L2);
mesh getCrossPoint(mesh M, double a, double b, double r);//某点和圆心的连线与圆的交点
Line getLine(mesh A, mesh B);
Line getLine(double x1, double y1, double x2, double y2);
bool judgeFieldInOut(mesh& A, vector<mesh>& Poly);
bool judgeFieldInOut(double x, double y, vector <mesh> &poly);

Line getLine(double theta, mesh A);
double compute_res();//计算残差
void sortPoint();
void update_p3(mesh& p);
void update_p4_s(mesh& p);
void update_p4_u(mesh& p);
void update_bound();
void polymesh();
int findNearPoint(mesh A, vector <mesh> &poly);//find the closest point of given grid point
