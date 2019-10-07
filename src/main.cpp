#include<iostream>
#include<fstream>
#include<cmath>
#include"include/const.h"
#include"include/functions.h"
#include"include/init.h"
#include"include/output.h"
#include"include/patition.h"
#include<stdlib.h>
#include<ctime>
#include<iomanip>
#include <vector>
#include<ctime>
#include<algorithm>
using namespace std;
vector <mesh> AP;//all points
vector<vector <mesh*>> A;
vector<mesh*> ps;//结构网格节点
vector<mesh*> pu;//非结构网格节点
vector<mesh*> bl;//左边界
vector<mesh*> br;//右边界
vector<mesh*> bu;//上边界
vector<mesh*> bd;//下边界
vector<mesh*> bb;//物体边界
vector <polygon_mesh> M0;
vector<vector <polygon_mesh>> M;
vector <mesh*> B;
vector <mesh*> C;
vector<vector <int>> ad;
vector<mesh> Ar;//记录上一时刻的物理量；

double dt;
double t_sim = 0;
int step = 0;
vector <polygon_mesh> PM;
vector<mesh> poly;
clock_t start;
clock_t finish;
double res;
int main()
{
	using namespace MeshPara;
	using ConstPara::t_end;
	//init_mesh1();
	init_mesh();
	polygonPoint(poly);
	init_polygon_mesh();
	getType();
	partition_Point();
	sortPoint();
	polymesh();
	//findConnectPoint();

	//out_neighbor();
	coordinate_trans();
	initFlow();
	ofstream fout;
	fout.open(to_string(Xnum) + "X" + to_string(Ynum) + ".dat", ios::app);
	fout << "variables=theta, delta_theta" << endl;
	fout << "zone T = \"dr = r +" << to_string(delta_r) + "r\"" << endl;
	vector <double> delta_theta;
	vector <double> delta_d;
	vector<double>theta;
	for (int i = 0; i < AP.size(); i++)
	{
		using namespace ConstPara;
		if (AP[i].type == "Body")
		{
			delta_d.push_back(distance(AP[i], *AP[i].neibor[0]));
			theta.push_back(get_theta(AP[i].x, AP[i].y, a, b));
			if (AP[i].x < a && AP[i].y > b)
				theta[theta.size() - 1] = pi + theta[theta.size() - 1];
			else if (AP[i].x < a && AP[i].y < b)
				theta[theta.size() - 1] = pi + theta[theta.size() - 1];
			else if (AP[i].x > a&& AP[i].y < b)
				theta[theta.size() - 1] = 2 * pi + theta[theta.size() - 1];
			else if (AP[i].x == a && AP[i].y > b)
				theta[theta.size() - 1] = pi / 2;
			else if (AP[i].x == a && AP[i].y < b)
				theta[theta.size() - 1] = 3 * pi / 2;
			else if (AP[i].x > a&& AP[i].y == b)
				theta[theta.size() - 1] = 0;
			else if (AP[i].x < a && AP[i].y == b)
				theta[theta.size() - 1] = pi;
		}
		sort(theta.begin(), theta.end());
		sort(delta_d.begin(), delta_d.end());
	}
	using ConstPara::pi;

	for (int i = 0; i < theta.size(); i++)
	{
		//fout << theta[i] << endl;
		if (i != theta.size() - 1)
			delta_theta.push_back(theta[i + 1] - theta[i]);
		else
			delta_theta.push_back(theta[0] + 2 * pi - theta[i]);
		fout << theta[i] * 180 / pi << "   " << delta_theta[i] * 180 / pi << endl;
	}
	fout.close();
	fout.clear();
	sort(delta_theta.begin(), delta_theta.end());
	//cout << delta_theta[0] * 180 / pi << "   " << delta_theta[delta_theta.size() - 1] * 180 / pi << endl;
	cout << delta_d[0] / dx << "dx   " << delta_d[delta_d.size() - 1] / dx << "dx" << endl;


	//out_M("mesh/" + methodType + "/step = " + to_string(step));
	start = clock();
	int i;
	while (t_sim < t_end)
	{
		record();
		get_dt();
		if (t_sim + dt > t_end)
			dt = t_end - t_sim;
		for (i = 0; i < ps.size(); i++)
		{
			if (ps[i]->sec_num == 0)
				continue;
			update_p4_s(*ps[i]);
		}
		//cout << "******************" << endl;

		for (i = 0; i < pu.size(); i++)
		{
			if (pu[i]->neibor.size() == 3)
				update_p3(*pu[i]);
			else
				update_p4_u(*pu[i]);
			AP[pu[i]->connectId] = *pu[i];//replace
		}

		update_bound();
		step++;
		t_sim = t_sim + dt;
		if (step % 100 == 0)
		{
			cout << t_sim << "  " << dt << "   " << step << endl;
			ofstream fout;
			fout.open("mesh/" + to_string(step) + ".dat");
			fout << "variables=x,y,rho,u,v,p" << endl;
			//fout << "solutiontime = "<<t_sim << endl;
			for (int i = 0; i < AP.size(); i++)
			{
				//if (AP[i].section != 0)
					fout << AP[i].x << "  " << AP[i].y << "  " << AP[i].section << "  " << AP[i].u << "  " << AP[i].v << "  " << AP[i].p << endl;
			}
			fout.close();
			fout.clear();
			fout.open("mesh/C/" + to_string(step) + ".dat");
			fout << "variables=x,y,rho,u,v,p" << endl;
			fout << "zone i = " << Xnum << " j = " << Ynum << " F = point" << endl;
			fout << "solutiontime = " << t_sim << endl;
			for (int i = 0; i < ps.size(); i++)
			{
				fout << ps[i]->x << "  " << ps[i]->y << "  " << ps[i]->rho << "  " << ps[i]->u << "  " << ps[i]->v << "  " << ps[i]->p << endl;
			}
			fout.close();
			fout.clear();
			out_M("mesh/step = " + to_string(step));

		}

		//if (step % 10 == 0)
		//{
		//	//if (abs(res - compute_res()) < 1e-10)
		//	//	break;
		//	//else
		//	res = compute_res();
		//	out_res();

		//	if (step % 1 == 0)
		//	{
		//		out_M("mesh/" + methodType + "/step = " + to_string(step));

		//		cout << "step = " << step << "   dt = " << dt << "   t = " << t_sim << "   CPU_t = " << (double)(double(clock()) - start) / CLOCKS_PER_SEC << "   res = " << res << endl;

		//		if (step % 100 == 0)
		//		{
		//			if (step % 100 == 0)
		//			{
		//				ofstream fout;
		//				fout.precision(15);
		//				int i, j;
		//				if (FlowType == "normal" || FlowType == "oblique")
		//				{
		//					fout.open("data/" + FlowType + methodType + "y=0.015step" + to_string(step) + ".dat");
		//					fout << "variables = x,rho" << endl;
		//					for (i = 0; i < A.size(); i++)
		//					{
		//						for (j = 0; j < A[i].size(); j++)
		//						{
		//							if (abs(A[i][j]->y - 0.015) < dx / 2)
		//								fout << A[i][j]->x << "   " << A[i][j]->rho << endl;
		//						}
		//					}
		//				}
		//				else if (FlowType == "intersection")
		//				{

		//					fout.open("data/" + FlowType + methodType + "_x=0.029 step" + to_string(step) + ".dat");
		//					fout << "variables = y,rho" << endl;
		//					for (i = 0; i < A.size(); i++)
		//					{
		//						for (j = 0; j < A[i].size(); j++)
		//						{
		//							if (abs(A[i][j]->x - 0.029) < dx / 2)
		//								fout << A[i][j]->y << "   " << A[i][j]->rho << endl;
		//						}
		//					}
		//					fout.close();
		//					fout.clear();
		//					fout.open("data/" + FlowType + methodType + "_x=0.0016 step" + to_string(step) + ".dat");
		//					fout << "variables = y,rho" << endl;
		//					for (i = 0; i < A.size(); i++)
		//					{
		//						for (j = 0; j < A[i].size(); j++)
		//						{
		//							if (abs(A[i][j]->x - 0.0016) < dx / 2)
		//								fout << A[i][j]->y << "   " << A[i][j]->rho << endl;
		//						}
		//					}

		//				}
		//			}
		//		}

		//	}
		//}


		//if (t_sim > 0.02)
		//	break;
		//if (step > 100000)
		//	break;
		//if (step > 1000 && res < 1e-10)
		//	break;
	}
	out_M("mesh/" + methodType + "/step = " + to_string(step));
	system("PAUSE");
}