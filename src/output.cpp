#include<fstream>
#include"include/const.h"
#include<iomanip>
#include<string>
#include"include/shockwave.h"
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
	extern vector <polygon_mesh> M0;
	ofstream fout;
	fout.open(name + ".dat");
	int i;
	//fout << "FILETYPE = GRID" << endl;
	//fout << "VARIABLES = \"X\", \"Y\"" << endl;
	fout << "VARIABLES =  \"X\", \"Y\"\"u\", \"v\", \"p\", \"rho\"" << endl;

	fout << "ZONE T=\"Test\"" << endl;
	fout << "ZONETYPE=FEPOLYGON" << endl;
	fout << "Nodes = " << Pnum << endl;
	fout << "Elements = " << M0.size() << endl;
	fout << "Faces = " << M0.size() * 6 << endl;
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

	for (i = 0; i < M0.size(); i++)
	{
		for (int j = 0; j < M0[i].face_start.size(); j++)
		{
			fout << M0[i].face_start[j] + 1 << "   " << M0[i].face_end[j] + 1 << endl;
		}
	}
	fout << endl;
	for (i = 0; i < M0.size(); i++)
	{
		for (int j = 0; j < M0[i].face_start.size(); j++)
			fout << i + 1 << "  ";
		fout << endl;
	}
	fout << endl;
	for (i = 0; i < M0.size(); i++)
	{
		for (int j = 0; j < M0[i].face_start.size(); j++)
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
	extern vector <polygon_mesh> M0;
	ofstream fout;
	fout.open(name + ".dat");
	int i;
	fout << "FILETYPE = SOLUTION" << endl;
	fout << "VARIABLES =  \"u\", \"v\", \"p\", \"rho\"" << endl;
	fout << "ZONE T=\"Test\"" << endl;
	fout << "ZONETYPE=FEPOLYGON" << endl;
	fout << "Nodes = " << Pnum << endl;
	fout << "Elements = " << M0.size() << endl;
	fout << "Faces = " << M0.size() * 6 << endl;
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
	extern vector<vector <mesh*>> A;
	extern double res;
	extern int step;
	extern double t_sim;
	ofstream fout;
	if (step == 0)
	{
		fout.open("res.dat");
		fout << "Variables= t,res" << endl;
	}
	else
	{
		fout.open("res.dat", ios::app);
		fout << t_sim << "   " << res << endl;
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