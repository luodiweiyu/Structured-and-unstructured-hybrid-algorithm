#include<fstream>
#include"/Structured-and-unstructured-hybrid-algorithm/include/const.h"
#include<iomanip>
#include<string>
#include"/Structured-and-unstructured-hybrid-algorithm/include/shockwave.h"
using namespace std;
using namespace MeshPara;
void out_mesh(string name)
{
	extern std::vector <mesh> A0;
	extern std::vector<std::vector<int>> ad;
	extern double t_sim;
	ofstream fout;
	fout.open(name + ".dat");
	int i;
	if (FlowType == "oblique" || FlowType == "intersection")
	{
		fout << "VARIABLES = \"X\", \"Y\", \"u\", \"v\", \"p\", \"老\"" << endl;

		fout << "ZONE I =" << Xnum << ", J = " << Ynum << ", F = point" << endl;
		fout << "solutiontime = " << t_sim << endl;
		fout.scientific;
		for (i = 0; i < Pnum; i++)
		{
			//if(A0[i].section==1)
			fout << A0[i].x << "   " << A0[i].y << "   " << A0[i].u << "   " << A0[i].v << "   " << A0[i].p << "   " << A0[i].老 << endl;
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
		fout << "VARIABLES = \"X\", \"Y\", \"u\", \"v\", \"p\", \"老\",\"Ma\"" << endl;


		fout << "ZONE N =" << A0.size() << ", E = " << s << ", F = FEPOINT, ET = TRIANGLE" << endl;
		fout << "solutiontime = " << t_sim << endl;
		for (i = 0; i < A0.size(); i++)
		{
			fout << A0[i].x << "   " << A0[i].y << "   " << A0[i].u << "   " << A0[i].v << "   " << A0[i].p << "   " << A0[i].老 << "   " << get_Ma(A0[i].u, A0[i].v, A0[i].老, A0[i].p) << endl;
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
			fout << A0[i].x << "   " << A0[i].y << "   " << A0[i].u << "   " << A0[i].v << "   " << A0[i].p << "   " << A0[i].老 << "   " << get_Ma(A0[i].u, A0[i].v, A0[i].老, A0[i].p) << endl;
		}
		fout << "ZONE I =" << Xnum << ", J = " << Ynum << ", F = point" << endl;
		fout << "solutiontime = " << t_sim << endl;
		fout.scientific;
		fout << A0[Xnum - 1].x << "   " << A0[Xnum - 1].y << "   " << A0[Xnum - 1].u << "   " << A0[Xnum - 1].v << "   " << A0[Xnum - 1].p << "   " << A0[Xnum - 1].老 << "   " << get_Ma(A0[Xnum - 1].u, A0[Xnum - 1].v, A0[Xnum - 1].老, A0[Xnum - 1].p) << endl;
		for (i = Pnum; i < A0.size(); i++)
		{
			fout << A0[i].x << "   " << A0[i].y << "   " << A0[i].u << "   " << A0[i].v << "   " << A0[i].p << "   " << A0[i].老 << "   " << get_Ma(A0[i].u, A0[i].v, A0[i].老, A0[i].p) << endl;
		}

	}
}

void out_Jacobin()
{
	extern std::vector <mesh> A0;
	ofstream out;
	out.open("J.dat");
	int i;
	for (i = 0; i < Pnum; i++)
	{
		for (int j = 0; j < A0[i].J.size(); j++)
			out << i << "  " << A0[i].J[j] << "   ";
		out << endl;
	}
}
void out_polygon_mesh(string name)
{
	extern std::vector <mesh> A0;
	extern double t_sim;
	extern vector <polygon_mesh> M0;
	ofstream fout;
	fout.open(name + ".dat");
	int i;
	//fout << "FILETYPE = GRID" << endl;
	//fout << "VARIABLES = \"X\", \"Y\"" << endl;
	fout << "VARIABLES =  \"X\", \"Y\"\"u\", \"v\", \"p\", \"老\"" << endl;

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
		fout << A0[i].x << endl;
	}
	fout << endl;

	for (i = 0; i < Pnum; i++)
	{
		fout << A0[i].y << endl;
	}
	fout << endl;
	for (i = 0; i < Pnum; i++)
	{
		fout << A0[i].u << endl;
	}
	fout << endl;

	for (i = 0; i < Pnum; i++)
	{
		fout << A0[i].v << endl;
	}
	fout << endl;

	for (i = 0; i < Pnum; i++)
	{
		fout << A0[i].p << endl;

	}
	fout << endl;

	for (i = 0; i < Pnum; i++)
	{
		fout << A0[i].老 << endl;
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
	extern vector<vector <mesh*>> A;
	extern vector<vector <polygon_mesh>> M;

	std::ofstream fout;
	extern double t_sim;
	using std::endl;
	fout.open(name + ".dat");
	fout << "VARIABLES =  \"X\", \"Y\", \"rho\", \"u\", \"v\", \"p\",\"um\",\"vm\",\"alpha1\",\"alpha2\",\"section\",\"secnum\"" << std::endl;

	int i, j, k;
	for (i = 0; i < M.size(); i++)
	{
		int face = 0;
		for (j = 0; j < M[i].size(); j++)
		{
			face += M[i][j].node.size();
		}

		fout << "ZONE T=\"Test" << i << "\"" << endl;
		fout << "ZONETYPE=FEPOLYGON" << endl;
		fout << "Nodes = " << A[i].size() << endl;
		fout << "Elements = " << M[i].size() << endl;
		fout << "Faces = " << face << endl;
		fout << "NumConnectedBoundaryFaces=0 " << endl;
		fout << "TotalNumBoundaryConnections=0 " << endl;
		fout << "solutiontime = " << t_sim << endl;

		fout.scientific;
		for (j = 0; j < A[i].size(); j++)
		{
			fout << A[i][j]->x << "   ";
			if (j % 30 == 0)
				fout << endl;
		}
		fout << endl;
		for (j = 0; j < A[i].size(); j++)
		{
			fout << A[i][j]->y << "   ";
			if (j % 30 == 0)
				fout << endl;
		}
		fout << endl;

		for (j = 0; j < A[i].size(); j++)
		{
			fout << A[i][j]->老 << "   ";
			if (j % 30 == 0)
				fout << endl;
		}
		fout << endl;
		for (j = 0; j < A[i].size(); j++)
		{
			fout << A[i][j]->u << "   ";
			if (j % 30 == 0)
				fout << endl;
		}
		fout << endl;
		for (j = 0; j < A[i].size(); j++)
		{
			fout << A[i][j]->v << "   ";
			if (j % 30 == 0)
				fout << endl;
		}
		fout << endl;
		for (j = 0; j < A[i].size(); j++)
		{
			fout << A[i][j]->p << "   ";
			if (j % 30 == 0)
				fout << endl;
		}
		fout << endl;
		for (j = 0; j < A[i].size(); j++)
		{

			fout << A[i][j]->um << "   ";
			if (j % 30 == 0)
				fout << endl;
		}
		fout << endl;
		for (j = 0; j < A[i].size(); j++)
		{
			fout << A[i][j]->vm << "   ";
			if (j % 30 == 0)
				fout << endl;
		}
		fout << endl;
		for (j = 0; j < A[i].size(); j++)
		{
			fout << A[i][j]->alpha.f1 << "   ";
			if (j % 30 == 0)
				fout << endl;
		}
		fout << endl;
		for (j = 0; j < A[i].size(); j++)
		{
			fout << A[i][j]->alpha.f2 << "   ";
			if (j % 30 == 0)
				fout << endl;
		}
		fout << endl;
		for (j = 0; j < A[i].size(); j++)
		{
			fout << A[i][j]->section << "   ";
			if (j % 30 == 0)
				fout << endl;
		}
		fout << endl;
		for (j = 0; j < A[i].size(); j++)
		{
			fout << A[i][j]->sec_num<< "   ";
			if (j % 30 == 0)
				fout << endl;
		}
		fout << endl;
		for (j = 0; j < M[i].size(); j++)
		{
			for (k = 0; k < M[i][j].face_start.size(); k++)
			{
				fout << M[i][j].face_start[k] + 1 << "   " << M[i][j].face_end[k] + 1 << endl;
			}
			fout << endl;
		}
		fout << endl;
		for (j = 0; j < M[i].size(); j++)
		{
			for (k = 0; k < M[i][j].face_start.size(); k++)
				fout << j + 1 << "  ";
			fout << endl;
		}
		fout << endl;
		for (j = 0; j < M[i].size(); j++)
		{
			for (k = 0; k < M[i][j].face_start.size(); k++)
				fout << 0 << "  ";
			fout << endl;
		}

	}
	//fout << "FILETYPE = GRID" << endl;
	//fout << "VARIABLES = \"X\", \"Y\"" << endl;



}

void out_polygon_variables(string name)
{
	extern std::vector <mesh> A0;
	extern double t_sim;
	extern vector <polygon_mesh> M0;
	ofstream fout;
	fout.open(name + ".dat");
	int i;
	fout << "FILETYPE = SOLUTION" << endl;
	fout << "VARIABLES =  \"u\", \"v\", \"p\", \"老\"" << endl;
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
		fout << A0[i].u << endl;
	}
	fout << endl;

	for (i = 0; i < Pnum; i++)
	{
		fout << A0[i].v << endl;
	}
	fout << endl;

	for (i = 0; i < Pnum; i++)
	{
		fout << A0[i].p << endl;

	}
	fout << endl;

	for (i = 0; i < Pnum; i++)
	{
		fout << A0[i].老 << endl;
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
//		fout << A[i].老 << " ";
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