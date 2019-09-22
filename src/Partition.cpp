//Partition file for computing , divid the original mesh into several regions
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

void partition_Point()//Partition existing grid points
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
			double a, b, r, r1;
			a = 0.02;
			b = 0.01;
			r = 0.55 / 75;
			//the 2D cylinder (x-a)^2+(y-b)^2=r^2
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
				r1 = r + sqrt(dx * dx + dy * dy);
				//r1 decides the unstructural grid region 
				//r1 is larger than r
				//the larger of r1, the larger of the unstructural grid region
				if ((xl - a) * (xl - a) + (yl - b) * (yl - b) > r1 * r1)
					n1 = 1;
				if ((xr - a) * (xr - a) + (yr - b) * (yr - b) > r1 * r1)
					n2 = 1;
				if ((xu - a) * (xu - a) + (yu - b) * (yu - b) > r1 * r1)
					n3 = 1;
				if ((xd - a) * (xd - a) + (yd - b) * (yd - b) > r1 * r1)
					n4 = 1;
				n = n1 + n2 + n3 + n4;
				if ((x - a) * (x - a) + (y - b) * (y - b) <= r1 * r1)
					A[0][i]->section = 0, A[0][i]->sec_num = 0;
				if ((x - a) * (x - a) + (y - b) * (y - b) > r1 * r1)
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
						A0.push_back(getCrossPoint(*A[1][size], a, b, r));
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
						A0.push_back(getCrossPoint(*A[2][size], a, b, r));
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

