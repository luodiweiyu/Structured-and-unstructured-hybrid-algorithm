#include<iostream>
#include<fstream>
#include<cmath>
#include"/Structured-and-unstructured-hybrid-algorithm/include/const.h"
#include"/Structured-and-unstructured-hybrid-algorithm/include/functions.h"
#include"/Structured-and-unstructured-hybrid-algorithm/include/init.h"
#include"/Structured-and-unstructured-hybrid-algorithm/include/output.h"
#include"/Structured-and-unstructured-hybrid-algorithm/include/patition.h"
#include<stdlib.h>
#include<ctime>
#include<iomanip>
#include <vector>
#include<ctime>
#include<algorithm>
using namespace std;
vector <mesh> A0;
vector<vector <mesh*>> A;
vector <mesh*> B;
vector <mesh*> C;
vector<vector <int>> ad;
vector<mesh> Ar;//记录上一时刻的物理量；

double dt;
double t_sim = 0;
int step = 0;
double xL, xR, yU, yD;
vector <polygon_mesh> M0;
vector<vector <polygon_mesh>> M;
clock_t start;
clock_t finish;
double res;
int main()
{
	srand((unsigned)time(NULL));
	using namespace MeshPara;
	using ConstPara::t_end;
	//init_mesh1();
	init_mesh();
	init_polygon_mesh();
	getType();
	reorderMesh();
	partition_Point();


	//findConnectPoint();

	//out_neighbor();
	reorder_neighbor();
	coordinate_trans();
	initFlow();
	ofstream fout;
	fout.open("mesh1.dat");
	fout << "variables=theta, delta_theta" << endl;
	vector <double> delta_theta;
	vector<double>theta;
	for (int i = 0; i < A0.size(); i++)
	{
		using namespace ConstPara;
		if (A0[i].type == "Cy")
		{
			theta.push_back(get_theta(A0[i].x, A0[i].y, a, b));
			if (A0[i].x < a && A0[i].y > b)
				theta[theta.size() - 1] = pi + theta[theta.size() - 1];
			else if (A0[i].x < a && A0[i].y < b)
				theta[theta.size() - 1] = pi + theta[theta.size() - 1];
			else if (A0[i].x > a && A0[i].y < b)
				theta[theta.size() - 1] = 2 * pi + theta[theta.size() - 1];
			else if (A0[i].x == a && A0[i].y > b)
				theta[theta.size() - 1] = pi / 2;
			else if (A0[i].x == a && A0[i].y < b)
				theta[theta.size() - 1] = 3 * pi / 2;
			else if (A0[i].x > a && A0[i].y == b)
				theta[theta.size() - 1] = 0;
			else if (A0[i].x < a && A0[i].y == b)
				theta[theta.size() - 1] = pi;

		}
		sort(theta.begin(), theta.end());
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
	cout << delta_theta[0] * 180 / pi << "   " << delta_theta[delta_theta.size() - 1] * 180 / pi << endl;


	//out_M("mesh/" + methodType + "/step = " + to_string(step));
	start = clock();
	while (t_sim < t_end)
	{
		record();

		get_dt();
		if (t_sim + dt > t_end)
			dt = t_end - t_sim;
		//if (t_sim < 1)
		//	dt = 2.5e-4;
		update_IN();
		update_bound_shockwave();
		step++;
		t_sim = t_sim + dt;
		if (step % 100 == 0)
		{
			cout << t_sim << "  " << dt << "   " << step << endl;
			ofstream fout;
			fout.open("mesh/" + to_string(step) + ".dat");
			fout << "variables=x,y,rho,u,v,p" << endl;
			//fout << "solutiontime = "<<t_sim << endl;
			for (int i = 0; i < A0.size(); i++)
			{
				if (A0[i].section != 0)
					fout << A0[i].x << "  " << A0[i].y << "  " << A0[i].rho << "  " << A0[i].u << "  " << A0[i].v << "  " << A0[i].p << endl;
			}
			fout.close();
			fout.clear();
			fout.open("mesh/C/" + to_string(step) + ".dat");
			fout << "variables=x,y,rho,u,v,p" << endl;
			fout << "zone i = 101 j=101 F=point" << endl;
			fout << "solutiontime = " << t_sim << endl;
			for (int i = 0; i < A[0].size(); i++)
			{
				fout << A[0][i]->x << "  " << A[0][i]->y << "  " << A[0][i]->rho << "  " << A[0][i]->u << "  " << A[0][i]->v << "  " << A[0][i]->p << endl;
			}
			fout.close();
			fout.clear();
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