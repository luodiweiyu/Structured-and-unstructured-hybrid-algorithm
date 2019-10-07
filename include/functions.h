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
Flux HLLC_��(mesh& CL, mesh& CR, mesh& C, int method);
Flux HLLC_��2(mesh& CL, mesh& CR, mesh& C, int method);//�������Ҳ��㣬����任�ο���

Flux HLLC_��(mesh& CD, mesh& CU, mesh& C, int method);
Flux HLLC_��2(mesh& CD, mesh& CU, mesh& C, int method);

Flux VanLeerA(mesh& C, double xix, double xiy, double xit, double J);
Flux VanLeerB(mesh& C, double xix, double xiy, double xit, double J);

void record();
double get_beta(mesh A, mesh B);//��������������x��ļн�
void reorder_neighbor();
double area(mesh A, mesh B, mesh C, mesh D);//�������ĵ㹹���ı������ 
void movemesh();
void findNeiborSec();
void update_Vm();
void clear_Vm();
void update_bound();
double get_theta(double x1, double y1, double x2, double y2);//��ֱ����x��ļн�

double max(double a, double b);
double min(double a, double b);
double absmax(double a, double b);
double absmin(double a, double b);
mesh getCrossPoint(double theta, double a, double b, double r);
void polygonPoint(vector <mesh> &poly);

mesh getCrossPoint(Line L1, Line L2);
mesh getCrossPoint(mesh M, double a, double b, double r);//ĳ���Բ�ĵ�������Բ�Ľ���
Line getLine(mesh A, mesh B);
Line getLine(double x1, double y1, double x2, double y2);
bool judgeFieldInOut(mesh& A, vector<mesh>& Poly);
bool judgeFieldInOut(double x, double y, vector <mesh> &poly);

Line getLine(double theta, mesh A);
double compute_res();//����в�
void sortPoint();
void update_p3(mesh& p);
void update_p4_s(mesh& p);
void update_p4_u(mesh& p);
void update_bound();
void polymesh();
int findNearPoint(mesh A, vector <mesh> &poly);//find the closest point of given grid point
