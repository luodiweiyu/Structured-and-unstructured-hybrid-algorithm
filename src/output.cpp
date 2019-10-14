#include<iostream>
#include<fstream>
#include<iomanip>
#include<string>
#include<algorithm>
#include"include/const.h"
#include"include/shockwave.h"
#include"include/functions.h"
using namespace std;
using namespace MeshPara;
void out_mesh(string name)
{
	extern std::vector <mesh> AP;
	extern std::vector<std::vector<int>> ad;
	extern double t_sim;
	ofstream fout;
	fout.open(name + ".dat");
	int i;
	if (FlowType == "oblique" || FlowType == "intersection")
	{
		fout << "VARIABLES = \"X\", \"Y\", \"u\", \"v\", \"p\", \"rho\"" << endl;

		fout << "ZONE I =" << Xnum << ", J = " << Ynum << ", F = point" << endl;
		fout << "solutiontime = " << t_sim << endl;
		fout.scientific;
		for (i = 0; i < Pnum; i++)
		{
			//if(AP[i].section==1)
			fout << AP[i].x << "   " << AP[i].y << "   " << AP[i].u << "   " << AP[i].v << "   " << AP[i].p << "   " << AP[i].rho << endl;
		}

	}
	else if (FlowType == "Prandtl-Meyer")
	{
		int s = 0;
		for (i = 0; i < ad.size(); i++)
		{
			if (ad[i].size() == 3)
				s++;
		}
		fout << "VARIABLES = \"X\", \"Y\", \"u\", \"v\", \"p\", \"rho\",\"Ma\"" << endl;


		fout << "ZONE N =" << AP.size() << ", E = " << s << ", F = FEPOINT, ET = TRIANGLE" << endl;
		fout << "solutiontime = " << t_sim << endl;
		for (i = 0; i < AP.size(); i++)
		{
			fout << AP[i].x << "   " << AP[i].y << "   " << AP[i].u << "   " << AP[i].v << "   " << AP[i].p << "   " << AP[i].rho << "   " << get_Ma(AP[i].u, AP[i].v, AP[i].rho, AP[i].p) << endl;
		}
		for (i = 0; i < ad.size(); i++)
		{
			if (ad[i].size() == 3)
				fout << ad[i][0] + 1 << "   " << ad[i][1] + 1 << "   " << ad[i][2] + 1 << endl;
		}

		fout << "ZONE I =" << Xnum << ", J = " << Ynum << ", F = point" << endl;
		fout << "solutiontime = " << t_sim << endl;
		fout.scientific;
		for (i = 0; i < Pnum; i++)
		{
			fout << AP[i].x << "   " << AP[i].y << "   " << AP[i].u << "   " << AP[i].v << "   " << AP[i].p << "   " << AP[i].rho << "   " << get_Ma(AP[i].u, AP[i].v, AP[i].rho, AP[i].p) << endl;
		}
		fout << "ZONE I =" << Xnum << ", J = " << Ynum << ", F = point" << endl;
		fout << "solutiontime = " << t_sim << endl;
		fout.scientific;
		fout << AP[Xnum - 1].x << "   " << AP[Xnum - 1].y << "   " << AP[Xnum - 1].u << "   " << AP[Xnum - 1].v << "   " << AP[Xnum - 1].p << "   " << AP[Xnum - 1].rho << "   " << get_Ma(AP[Xnum - 1].u, AP[Xnum - 1].v, AP[Xnum - 1].rho, AP[Xnum - 1].p) << endl;
		for (i = Pnum; i < AP.size(); i++)
		{
			fout << AP[i].x << "   " << AP[i].y << "   " << AP[i].u << "   " << AP[i].v << "   " << AP[i].p << "   " << AP[i].rho << "   " << get_Ma(AP[i].u, AP[i].v, AP[i].rho, AP[i].p) << endl;
		}

	}
}

void out_Jacobin()
{
	extern std::vector <mesh> AP;
	ofstream out;
	out.open("J.dat");
	int i;
	for (i = 0; i < Pnum; i++)
	{
		for (int j = 0; j < AP[i].J.size(); j++)
			out << i << "  " << AP[i].J[j] << "   ";
		out << endl;
	}
}
void out_polygon_mesh(string name)
{
	extern std::vector <mesh> AP;
	extern double t_sim;
	extern vector <polygon_mesh> PM;
	ofstream fout;
	fout.open(name + ".dat");
	int i;
	//fout << "FILETYPE = GRID" << endl;
	//fout << "VARIABLES = \"X\", \"Y\"" << endl;
	fout << "VARIABLES =  \"X\", \"Y\"\"u\", \"v\", \"p\", \"rho\"" << endl;

	fout << "ZONE T=\"Test\"" << endl;
	fout << "ZONETYPE=FEPOLYGON" << endl;
	fout << "Nodes = " << Pnum << endl;
	fout << "Elements = " << PM.size() << endl;
	fout << "Faces = " << PM.size() * 6 << endl;
	fout << "NumConnectedBoundaryFaces=0 " << endl;
	fout << "TotalNumBoundaryConnections=0 " << endl;
	fout << "solutiontime = " << t_sim << endl;

	fout.scientific;
	for (i = 0; i < Pnum; i++)
	{
		fout << AP[i].x << endl;
	}
	fout << endl;

	for (i = 0; i < Pnum; i++)
	{
		fout << AP[i].y << endl;
	}
	fout << endl;
	for (i = 0; i < Pnum; i++)
	{
		fout << AP[i].u << endl;
	}
	fout << endl;

	for (i = 0; i < Pnum; i++)
	{
		fout << AP[i].v << endl;
	}
	fout << endl;

	for (i = 0; i < Pnum; i++)
	{
		fout << AP[i].p << endl;

	}
	fout << endl;

	for (i = 0; i < Pnum; i++)
	{
		fout << AP[i].rho << endl;
	}
	fout << endl;

	for (i = 0; i < PM.size(); i++)
	{
		for (int j = 0; j < PM[i].face_start.size(); j++)
		{
			fout << PM[i].face_start[j] + 1 << "   " << PM[i].face_end[j] + 1 << endl;
		}
	}
	fout << endl;
	for (i = 0; i < PM.size(); i++)
	{
		for (int j = 0; j < PM[i].face_start.size(); j++)
			fout << i + 1 << "  ";
		fout << endl;
	}
	fout << endl;
	for (i = 0; i < PM.size(); i++)
	{
		for (int j = 0; j < PM[i].face_start.size(); j++)
			fout << 0 << "  ";
		fout << endl;
	}

}
void out_M(std::string name)
{
	extern vector <mesh> AP;
	extern vector <polygon_mesh> PM;

	std::ofstream fout;
	extern double t_sim;
	using std::endl;
	fout.open(name + ".dat");
	fout << "VARIABLES =  \"X\", \"Y\", \"rho\", \"u\", \"v\", \"p\",\"um\",\"vm\",\"alpha1\",\"alpha2\",\"section\",\"secnum\"" << std::endl;

	int i, j, k;
	int face = 0;
	for (j = 0; j < PM.size(); j++)
	{
		face += PM[j].node.size();
	}

	fout << "ZONE T=\"Test" << 0 << "\"" << endl;
	fout << "ZONETYPE=FEPOLYGON" << endl;
	fout << "Nodes = " << AP.size() << endl;
	fout << "Elements = " << PM.size() << endl;
	fout << "Faces = " << face << endl;
	fout << "NumConnectedBoundaryFaces=0 " << endl;
	fout << "TotalNumBoundaryConnections=0 " << endl;
	fout << "solutiontime = " << t_sim << endl;

	fout.scientific;
	for (j = 0; j < AP.size(); j++)
	{
		fout << AP[j].x << "   ";
		if (j % 30 == 0)
			fout << endl;
	}
	fout << endl;
	for (j = 0; j < AP.size(); j++)
	{
		fout << AP[j].y << "   ";
		if (j % 30 == 0)
			fout << endl;
	}
	fout << endl;

	for (j = 0; j < AP.size(); j++)
	{
		fout << AP[j].rho << "   ";
		if (j % 30 == 0)
			fout << endl;
	}
	fout << endl;
	for (j = 0; j < AP.size(); j++)
	{
		fout << AP[j].u << "   ";
		if (j % 30 == 0)
			fout << endl;
	}
	fout << endl;
	for (j = 0; j < AP.size(); j++)
	{
		fout << AP[j].v << "   ";
		if (j % 30 == 0)
			fout << endl;
	}
	fout << endl;
	for (j = 0; j < AP.size(); j++)
	{
		fout << AP[j].p << "   ";
		if (j % 30 == 0)
			fout << endl;
	}
	fout << endl;
	for (j = 0; j < AP.size(); j++)
	{

		fout << AP[j].um << "   ";
		if (j % 30 == 0)
			fout << endl;
	}
	fout << endl;
	for (j = 0; j < AP.size(); j++)
	{
		fout << AP[j].vm << "   ";
		if (j % 30 == 0)
			fout << endl;
	}
	fout << endl;
	for (j = 0; j < AP.size(); j++)
	{
		fout << AP[j].alpha.f1 << "   ";
		if (j % 30 == 0)
			fout << endl;
	}
	fout << endl;
	for (j = 0; j < AP.size(); j++)
	{
		fout << AP[j].alpha.f2 << "   ";
		if (j % 30 == 0)
			fout << endl;
	}
	fout << endl;
	for (j = 0; j < AP.size(); j++)
	{
		fout << AP[j].section << "   ";
		if (j % 30 == 0)
			fout << endl;
	}
	fout << endl;
	for (j = 0; j < AP.size(); j++)
	{
		fout << AP[j].sec_num << "   ";
		if (j % 30 == 0)
			fout << endl;
	}
	fout << endl;
	for (j = 0; j < PM.size(); j++)
	{
		for (k = 0; k < PM[j].face_start.size(); k++)
		{
			fout << PM[j].face_start[k] + 1 << "   " << PM[j].face_end[k] + 1 << endl;
		}
		fout << endl;
	}
	fout << endl;
	for (j = 0; j < PM.size(); j++)
	{
		for (k = 0; k < PM[j].face_start.size(); k++)
			fout << j + 1 << "  ";
		fout << endl;
	}
	fout << endl;
	for (j = 0; j < PM.size(); j++)
	{
		for (k = 0; k < PM[j].face_start.size(); k++)
			fout << 0 << "  ";
		fout << endl;
	}

	//fout << "FILETYPE = GRID" << endl;
//fout << "VARIABLES = \"X\", \"Y\"" << endl;



}

void out_polygon_variables(string name)
{
	extern std::vector <mesh> AP;
	extern double t_sim;
	extern vector <polygon_mesh> PM;
	ofstream fout;
	fout.open(name + ".dat");
	int i;
	fout << "FILETYPE = SOLUTION" << endl;
	fout << "VARIABLES =  \"u\", \"v\", \"p\", \"rho\"" << endl;
	fout << "ZONE T=\"Test\"" << endl;
	fout << "ZONETYPE=FEPOLYGON" << endl;
	fout << "Nodes = " << Pnum << endl;
	fout << "Elements = " << PM.size() << endl;
	fout << "Faces = " << PM.size() * 6 << endl;
	fout << "NumConnectedBoundaryFaces=0 " << endl;
	fout << "TotalNumBoundaryConnections=0 " << endl;
	fout << "solutiontime = " << t_sim << endl;
	fout.precision(10);

	for (i = 0; i < Pnum; i++)
	{
		fout << AP[i].u << endl;
	}
	fout << endl;

	for (i = 0; i < Pnum; i++)
	{
		fout << AP[i].v << endl;
	}
	fout << endl;

	for (i = 0; i < Pnum; i++)
	{
		fout << AP[i].p << endl;

	}
	fout << endl;

	for (i = 0; i < Pnum; i++)
	{
		fout << AP[i].rho << endl;
	}
	fout << endl;
}
void out_res()
{
	extern double res;
	extern int step;
	extern double t_sim;
	ofstream fout;
	if (t_sim == 0)
	{
		fout.open("res.dat");
		fout << "Variables= t,res" << endl;
	}
	else if (step % 100 == 0)
	{
		fout.open("res.dat", ios::app);
		fout << t_sim << "   " << res << endl;
	}
	fout.close();
}
void outNeiborLines(vector<mesh*> m, string filename)
{
	int i, j;
	extern vector<mesh>AP;
	ofstream fout(filename);
	int E = 0;
	for (i = 0; i < m.size(); i++)
	{
		for (j = 0; j < m[i]->neibor.size(); j++)
			E++;
	}
	fout << "variables = x,y" << endl;
	fout << "ZONE N = " << AP.size() << ", E = " << E << ", F = FEPOINT, ET = TRIANGLE" << endl;
	for (i = 0; i < AP.size(); i++)
	{
		fout << AP[i].x << "  " << AP[i].y << endl;
	}
	for (i = 0; i < m.size(); i++)
	{
		for (j = 0; j < m[i]->neibor.size(); j++)
			fout << m[i]->id + 1 << "  " << m[i]->id + 1 << "  " << m[i]->neibor[j]->id + 1 << endl;
	}
}
//void outmesh_polygon(string name)
//{
//	extern std::vector <mesh> A;
//	extern double t_sim;
//	extern vector <polygon_mesh> M;
//	ofstream fout;
//	fout.open(name + ".dat");
//	int i;
//	fout << "VARIABLES = \"X\", \"Y\", \"u\", \"v\", \"p\", \"rho\"" << endl;
//	fout << "ZONE T=\"Test\"" << endl;
//	fout << "ZONETYPE=FEPOLYGON" << endl;
//	fout << "Nodes = " << Pnum << endl;
//	fout << "Elements = " << M.size() << endl;
//	fout << "Faces = " << M.size() * 6 << endl;
//	fout << "NumConnectedBoundaryFaces=0 " << endl;
//	fout << "TotalNumBoundaryConnections=0 " << endl;
//	fout << "solutiontime = " << t_sim << endl;
//	fout.scientific;
//	for (i = 0; i < Pnum; i++)
//	{
//		fout << A[i].x << " ";
//		if (i % 30 == 0)
//			fout << endl;
//	}
//	fout << endl;
//	fout << endl;
//
//	for (i = 0; i < Pnum; i++)
//	{
//		fout << A[i].y << " ";
//		if (i % 30 == 0)
//			fout << endl;
//	}
//	fout << endl;
//
//	for (i = 0; i < Pnum; i++)
//	{
//		fout << A[i].u << " ";
//		if (i % 30 == 0)
//			fout << endl;
//	}
//	fout << endl;
//
//	for (i = 0; i < Pnum; i++)
//	{
//		fout << A[i].v << " ";
//		if (i % 30 == 0)
//			fout << endl;
//	}
//	fout << endl;
//
//	for (i = 0; i < Pnum; i++)
//	{
//		fout << A[i].p << " ";
//		if (i % 30 == 0)
//			fout << endl;
//	}
//	fout << endl;
//
//	for (i = 0; i < Pnum; i++)
//	{
//		fout << A[i].rho << " ";
//		if (i % 30 == 0)
//			fout << endl;
//	}
//	fout << endl;
//
//	for (i = 0; i < M.size(); i++)
//	{
//		for (int j = 0; j < M[i].face_start.size(); j++)
//		{
//			fout << M[i].face_start[j] + 1 << "   " << M[i].face_end[j] + 1 << endl;
//		}
//	}
//	fout << endl;
//	for (i = 0; i < M.size(); i++)
//	{
//		for (int j = 0; j < M[i].face_start.size(); j++)
//			fout << i + 1 << "  ";
//		if (i % 30 == 0)
//			fout << endl;
//	}
//	fout << endl;
//	for (i = 0; i < M.size(); i++)
//	{
//		for (int j = 0; j < M[i].face_start.size(); j++)
//			fout << 0 << "  ";
//		if (i % 30 == 0)
//			fout << endl;
//
//	}
//
//}
void out_deltat_theta_d()
{
	ofstream fout;
	extern vector<mesh>AP;
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

}
void out_somepoint()
{
	extern vector<mesh*> ps;//结构网格节点
	extern vector<mesh*> pu;//非结构网格节点
	extern vector<mesh*> bl;//左边界
	extern vector<mesh*> br;//右边界
	extern vector<mesh*> bu;//上边界
	extern vector<mesh*> bd;//下边界
	extern vector<mesh*> bb;//物体边界
	extern vector<mesh> poly;
	extern vector<mesh>AP;
	outNeiborLines(ps, "ps.dat");
	outNeiborLines(pu, "pu.dat");
	ofstream fout;
	fout.open("poly.dat");
	for (int i = 0; i < poly.size(); i++)
		fout << poly[i].x << "  " << poly[i].y << endl;
	fout.open("pu_point.dat");
	fout << "variables = x, y" << endl;
	for (int i = 0; i < pu.size(); i++)
		fout << pu[i]->x << "  " << pu[i]->y << endl;
	fout.close();
	fout.open("ps_point.dat");
	fout << "variables = x, y" << endl;
	for (int i = 0; i < ps.size(); i++)
		if (ps[i]->section == 0)
			fout << ps[i]->x << "  " << ps[i]->y << endl;
	fout.close();
	fout.open("bb_point.dat");
	fout << "variables = x, y" << endl;
	for (int i = 0; i < bb.size(); i++)
		fout << bb[i]->x << "  " << bb[i]->y << endl;
	fout.close();
	fout.open("bb_neighbor.dat");
	int E = 0;
	for (int i = 0; i < bb.size(); i++)
	{
		for (int j = 0; j < bb[i]->neibor.size(); j++)
			if (bb[i]->neibor[j]->type == "Body")
				E++;
	}
	fout << "variables = x,y" << endl;
	fout << "ZONE N = " << AP.size() << ", E = " << E << ", F = FEPOINT, ET = TRIANGLE" << endl;
	for (int i = 0; i < AP.size(); i++)
		fout << AP[i].x << "  " << AP[i].y << endl;
	for (int i = 0; i < bb.size(); i++)
	{

		for (int j = 0; j < bb[i]->neibor.size(); j++)
			if (bb[i]->neibor[j]->type == "Body")
				fout << bb[i]->id << "  " << bb[i]->id << "  " << bb[i]->neibor[j]->id << endl;
	}
	fout.close();

}