//Partition file for computing , divid the original mesh into several regions
#include<iostream>
#include"include/const.h"
#include"include/functions.h"
#include<vector>
#include<string>
#include<fstream>
using std::vector;
extern vector <mesh> AP;
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
		if (FlowType == "cylinder")
		{
			extern vector <mesh> poly;
			double x, y, xl, yl, xr, yr, xu, yu, xd, yd;
			double r1;
			//the 2D cylinder (x-a)^2+(y-b)^2=r^2
			int n1, n2, n3, n4, n, sizeA, size;
			sizeA = AP.size();
			for (i = 0; i < sizeA; i++)
			{
				AP[i].id = i;
				if (AP[i].type == "IN")
				{
					AP[i].neibor.push_back(&AP[i + 1]);
					AP[i].neibor.push_back(&AP[i + Xnum]);
					AP[i].neibor.push_back(&AP[i - 1]);
					AP[i].neibor.push_back(&AP[i - Xnum]);
				}
				x = AP[i].x;
				y = AP[i].y;
				xl = x - dx;
				yl = y;
				xr = x + dx;
				yr = y;
				xu = x;
				yu = y + dy;
				xd = x;
				yd = y - dy;
				n1 = n2 = n3 = n4 = 0;
				if (i == 492)
					i = 492;
				//r1 = r + delta_r * dx;
				//r1 decides the unstructural grid region 
				//r1 is larger than r
				//the larger of r1, the larger of the unstructural grid region
				//if ((xl - a) * (xl - a) + (yl - b) * (yl - b) > r1* r1)
				//	n1 = 1;
				//if ((xr - a) * (xr - a) + (yr - b) * (yr - b) > r1* r1)
				//	n2 = 1;
				//if ((xu - a) * (xu - a) + (yu - b) * (yu - b) > r1* r1)
				//	n3 = 1;
				//if ((xd - a) * (xd - a) + (yd - b) * (yd - b) > r1* r1)
				//	n4 = 1;
				//if ((x - a) * (x - a) + (y - b) * (y - b) <= r1 * r1)


				//if (judgeFieldInOut(AP[i], poly))
				//	AP[i].section = 0, AP[i].sec_num = 0;
				//else if (!judgeFieldInOut(AP[i], poly))

				//{
					//if (!judgeFieldInOut(xl, yl, poly))
					//	n1 = 1;
					//if (!judgeFieldInOut(xr, yr, poly))
					//	n2 = 1;
					//if (!judgeFieldInOut(xu, yu, poly))
					//	n3 = 1;
					//if (!judgeFieldInOut(xd, yd, poly))
					//	n4 = 1;
				if (judgeFieldInOut(AP[i], poly))
					AP[i].section = 0, AP[i].sec_num = 0;
				else
				{

					if (!judgeFieldInOut(xl, yl, poly))
						n1 = 1;
					if (!judgeFieldInOut(xr, yr, poly))
						n2 = 1;
					if (!judgeFieldInOut(xu, yu, poly))
						n3 = 1;
					if (!judgeFieldInOut(xd, yd, poly))
						n4 = 1;

					n = n1 + n2 + n3 + n4;
					if (n == 4)
						AP[i].section = 1, AP[i].sec_num = 1;
					else
					{
						double r1 = r;
						int near_id = findNearPoint(AP[i], poly);
						AP[i].section = -1, AP[i].sec_num = 0;
						AP.push_back(AP[i]);
						size = AP.size() - 1;
						AP[size].id = size;
						AP[size].connectId = i;
						AP[i].connectId = size;
						vector<mesh*> v;
						AP[size].neibor.swap(v);
						AP.push_back(poly[near_id]);

						//if (AP[i].x < a)
						//{
						//	AP.push_back(getCrossPoint(AP[AP.size() - 1], a, b, r1));//create a new point on the body
						//}
						//else if (AP[i].y > b)
						//{
						//	double x0 = a, y0 = b + r1;
						//	double k = tan(4.6 * ConstPara::pi / 180);
						//	Line L1, L2;
						//	L1.A = k, L1.B = -1, L1.C = y0 - k * x0;
						//	L2.A = -1 / k, L2.B = -1, L2.C = AP[i].y - (-1 / k) * AP[i].x;

						//	AP.push_back(getCrossPoint(L1, L2));
						//}
						//else if (AP[i].y < b)
						//{
						//	double x0 = a, y0 = b - r1;
						//	double k = tan(-4.6 * ConstPara::pi / 180);
						//	Line L1, L2;
						//	L1.A = k, L1.B = -1, L1.C = y0 - k * x0;
						//	L2.A = -1 / k, L2.B = -1, L2.C = AP[i].y - (-1 / k) * AP[i].x;
						//	AP.push_back(getCrossPoint(L1, L2));
						//}
						AP[AP.size() - 1].id = AP.size() - 1;
						AP[AP.size() - 1].type = "Body";
						AP[AP.size() - 1].neibor.push_back(&AP[size]);
						AP[AP.size() - 1].section = -2;
						if (n == 3)
						{
							if (n1 == 0)
							{
								AP[size].neibor.push_back(&AP[i + 1]);
								AP[size].neibor.push_back(&AP[i + Xnum]);
								AP[size].neibor.push_back(&AP[AP.size() - 1]);
								AP[size].neibor.push_back(&AP[i - Xnum]);
							}
							else if (n2 == 0)
							{
								AP[size].neibor.push_back(&AP[AP.size() - 1]);
								AP[size].neibor.push_back(&AP[i + Xnum]);
								AP[size].neibor.push_back(&AP[i - 1]);
								AP[size].neibor.push_back(&AP[i - Xnum]);
							}
							else if (n3 == 0)
							{
								AP[size].neibor.push_back(&AP[i + 1]);
								AP[size].neibor.push_back(&AP[AP.size() - 1]);
								AP[size].neibor.push_back(&AP[i - 1]);
								AP[size].neibor.push_back(&AP[i - Xnum]);
							}
							else if (n4 == 0)
							{
								AP[size].neibor.push_back(&AP[i + 1]);
								AP[size].neibor.push_back(&AP[i + Xnum]);
								AP[size].neibor.push_back(&AP[i - 1]);
								AP[size].neibor.push_back(&AP[AP.size() - 1]);
							}
						}
						else if (n == 2)
						{
							if (n1 == 0 && n3 == 0)
							{
								AP[size].neibor.push_back(&AP[i + 1]);
								AP[size].neibor.push_back(&AP[AP.size() - 1]);
								AP[size].neibor.push_back(&AP[i - Xnum]);
							}
							if (n1 == 0 && n4 == 0)
							{
								AP[size].neibor.push_back(&AP[i + 1]);
								AP[size].neibor.push_back(&AP[i + Xnum]);
								AP[size].neibor.push_back(&AP[AP.size() - 1]);
							}
							if (n2 == 0 && n3 == 0)
							{
								AP[size].neibor.push_back(&AP[AP.size() - 1]);
								AP[size].neibor.push_back(&AP[i - 1]);
								AP[size].neibor.push_back(&AP[i - Xnum]);
							}
							if (n2 == 0 && n4 == 0)
							{
								AP[size].neibor.push_back(&AP[i + Xnum]);
								AP[size].neibor.push_back(&AP[i - 1]);
								AP[size].neibor.push_back(&AP[AP.size() - 1]);
							}
						}
					}
				}
			}
		}
	}


}

