//分区文件，用于激波装配法，将原网格分成若干区域
#include<iostream>
#include"/Structured-and-unstructured-hybrid-algorithm/include/const.h"
#include"/Structured-and-unstructured-hybrid-algorithm/include/functions.h"
#include<vector>
#include<string>
#include<fstream>
using std::vector;
extern vector <mesh> A0;
extern vector <polygon_mesh> M0;
using namespace ShockwaveCross;
using MeshPara::Pnum;

extern vector<vector <mesh*>> A;
extern vector<vector <polygon_mesh>> M;

Line L12;
Line L13;
int findAd(mesh A, vector <mesh> An)
{
	int s = -1;
	for (int i = 0; i < An.size(); i++)
	{
		if (A.x == An[i].x && A.y == An[i].y)
			s = i;
		else
			continue;
	}
	return s;
}

void partition_Point()//对已有的网格点进行分区
{
	using namespace MeshPara;
	if (methodType == "C")
	{
		vector<mesh*> t;
		int i;
		if (FlowType != "cylinder")
		{
			A.push_back(t);
			for (i = 0; i < A0.size(); i++)
				A[0].push_back(&A0[i]);
		}
		if (FlowType == "cylinder")
		{
			A.push_back(t);
			A.push_back(t);
			A.push_back(t);
			for (i = 0; i < A0.size(); i++)
				A[0].push_back(&A0[i]), A[0][i]->id = i;
			double x, y, xl, yl, xr, yr, xu, yu, xd, yd;
			int n1, n2, n3, n4, n, size;
			for (i = 0; i < A[0].size(); i++)
			{
				x = A[0][i]->x;
				y = A[0][i]->y;
				xl = x - dx;
				yl = y;
				xr = x + dx;
				yr = y;
				xu = x;
				yu = y + dy;
				xd = x;
				yd = y - dy;
				n1 = n2 = n3 = n4 = 0;
				if ((xl - 0.02)*(xl - 0.02) + (yl - 0.01)*(yl - 0.01) > (0.55 / 75)*(0.55 / 75))
					n1 = 1;
				if ((xr - 0.02)*(xr - 0.02) + (yr - 0.01)*(yr - 0.01) > (0.55 / 75)*(0.55 / 75))
					n2 = 1;
				if ((xu - 0.02)*(xu - 0.02) + (yu - 0.01)*(yu - 0.01) > (0.55 / 75)*(0.55 / 75))
					n3 = 1;
				if ((xd - 0.02)*(xd - 0.02) + (yd - 0.01)*(yd - 0.01) > (0.55 / 75)*(0.55 / 75))
					n4 = 1;
				n = n1 + n2 + n3 + n4;
				if ((x - 0.02)*(x - 0.02) + (y - 0.01)*(y - 0.01) <= (0.55 / 75)*(0.55 / 75))
					A[0][i]->section = 0, A[0][i]->sec_num = 0;
				if ((x - 0.02)*(x - 0.02) + (y - 0.01)*(y - 0.01) > (0.55 / 75)*(0.55 / 75))
				{
					if (n == 4)
					{
						A[0][i]->section = 1, A[0][i]->sec_num = 1;
						if (A[0][i]->type == "IN")
						{
							A[0][i]->neibor.push_back(A[0][i + 1]);
							A[0][i]->neibor.push_back(A[0][i + Xnum]);
							A[0][i]->neibor.push_back(A[0][i - 1]);
							A[0][i]->neibor.push_back(A[0][i - Xnum]);
						}
					}
					else if (n == 3)
					{
						A[1].push_back(&A0[i]);
						A[0][i]->section = -1, A[0][i]->sec_num = 0;
						size = A[1].size() - 1;
						A[1][size]->neiborsec = 0;
						A[1][size]->neiborsec_ad = i;
						A0.push_back(getCrossPoint(*A[1][size], 0.02, 0.01, 0.55 / 75));
						A[1][size]->id = i;
						A0[A0.size() - 1].id = A0.size() - 1;
						A0[A0.size() - 1].type = "Cy";
						A0[A0.size() - 1].neibor.push_back(A[1][size]);
						if (n1 == 0)
						{
							A[1][size]->neibor.push_back(&A0[i + 1]);
							A[1][size]->neibor.push_back(&A0[i + Xnum]);
							A[1][size]->neibor.push_back(&A0[A0.size() - 1]);
							A[1][size]->neibor.push_back(&A0[i - Xnum]);
						}
						else if (n2 == 0)
						{ 
							A[1][size]->neibor.push_back(&A0[A0.size() - 1]);
							A[1][size]->neibor.push_back(&A0[i + Xnum]);
							A[1][size]->neibor.push_back(&A0[i - 1]);
							A[1][size]->neibor.push_back(&A0[i - Xnum]);
						}
						else if (n3 == 0)
						{
							A[1][size]->neibor.push_back(&A0[i + 1]);
							A[1][size]->neibor.push_back(&A0[A0.size() - 1]);
							A[1][size]->neibor.push_back(&A0[i - 1]);
							A[1][size]->neibor.push_back(&A0[i - Xnum]);
						}
						else if (n4 == 0)
						{
							A[1][size]->neibor.push_back(&A0[i + 1]);
							A[1][size]->neibor.push_back(&A0[i + Xnum]);
							A[1][size]->neibor.push_back(&A0[i - 1]);
							A[1][size]->neibor.push_back(&A0[A0.size() - 1]);
						}
					}
					else if (n == 2)
					{
						A0[i].section = -1, A0[i].sec_num = 0;
						A[2].push_back(&A0[i]);
						size = A[2].size() - 1;
						A[2][size]->neiborsec = 0;
						A[2][size]->neiborsec_ad = i;
						A0.push_back(getCrossPoint(*A[2][size], 0.02, 0.01, 0.55 / 75));
						A[2][size]->id = i;
						A0[A0.size() - 1].id = A0.size() - 1;
						A0[A0.size() - 1].type = "Cy";
						A0[A0.size() - 1].neibor.push_back(A[2][size]);
						if (A[0][i]->type == "IN")
						{
							if (n1 == 0 && n3 == 0)
							{
								A[2][size]->neibor.push_back(&A0[i + 1]);
								A[2][size]->neibor.push_back(&A0[A0.size() - 1]);
								A[2][size]->neibor.push_back(&A0[i - Xnum]);
							}
							if (n1 == 0 && n4 == 0)
							{
								A[2][size]->neibor.push_back(&A0[i + 1]);
								A[2][size]->neibor.push_back(&A0[i + Xnum]);
								A[2][size]->neibor.push_back(&A0[A0.size() - 1]);
							}
							if (n2 == 0 && n3 == 0)
							{
								A[2][size]->neibor.push_back(&A0[A0.size() - 1]);
								A[2][size]->neibor.push_back(&A0[i - 1]);
								A[2][size]->neibor.push_back(&A0[i - Xnum]);
							}
							if (n2 == 0 && n4 == 0)
							{
								A[2][size]->neibor.push_back(&A0[i + Xnum]);
								A[2][size]->neibor.push_back(&A0[i - 1]);
								A[2][size]->neibor.push_back(&A0[A0.size() - 1]);
							}

						}
					}
				}

			}
		}
	}


}

