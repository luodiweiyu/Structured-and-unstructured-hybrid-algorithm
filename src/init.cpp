#include"include/const.h"
#include"include/functions.h"
#include"include/shockwave.h"
#include<iostream>
#include<ctime>
#include<stdlib.h>
#include<vector>
#include"include/Prandtl-Meyer.h"
#include"include/init.h"

using std::vector;
using namespace MeshPara;
void init_mesh()
{
	extern vector <mesh> AP;
	extern vector<vector<int>> ad;
	mesh t;
	vector<int> a;
	int i, j;
	i = 0;

	if (FlowType == "oblique" || FlowType == "intersection" || FlowType == "normal" || FlowType == "couette" || FlowType == "cylinder")
	{
		if (meshType == 10 || meshType == 11 || meshType == 12)
		{
			while (i < Pnum)
			{
				for (j = 0; j < Xnum; j++)
				{
					AP.push_back(t);
					if (j == 0)
					{
						if (i == 0)
							AP[i].x = dx / 2, AP[i].y = 0;
						else if (AP[i - Xnum].x == 0)
							AP[i].x = dx / 2, AP[i].y = AP[i - Xnum].y + dy;
						else
							AP[i].x = 0, AP[i].y = AP[i - Xnum].y + dy;
					}
					else if (AP[i - j].x == 0)
					{
						if (j % 2 == 0)
							AP[i].x = AP[i - 1].x + dx, AP[i].y = AP[i - 1].y;
						else
							AP[i].x = AP[i - 1].x + 2 * dx, AP[i].y = AP[i - 1].y;
					}
					else
					{
						if (j % 2 == 0)
							AP[i].x = AP[i - 1].x + 2 * dx, AP[i].y = AP[i - 1].y;
						else
							AP[i].x = AP[i - 1].x + dx, AP[i].y = AP[i - 1].y;
					}
					AP[i].id = i;
					i++;
				}
			}
		}
		else
		{
			while (i < Pnum)
			{
				for (j = 0; j < Xnum; j++)
				{
					AP.push_back(t);

					if (j == 0)
					{
						if (i == 0)
							AP[i].x = 0, AP[i].y = 0;
						else
							AP[i].x = 0, AP[i].y = AP[i - Xnum].y + dy;
					}
					else
						AP[i].x = AP[i - 1].x + dx, AP[i].y = AP[i - 1].y;
					i++;
				}
			}
		}
	}
	if (FlowType == "Prandtl-Meyer")
	{
		using namespace Prandtl_Meyer;
		using ConstPara::gama;
		double LineLength = 0.003;
		//左边结构网格生成
		double dx;
		double dl = LineLength / Ynum;
		double ex = cos(mu1);
		double ey = sin(mu1);
		double xend, yend;
		while (i < Pnum)
		{
			if (i == 0)
			{
				xend = LineLength;
				yend = 0;
			}
			else
			{
				xend = dl * ex + AP[i - 1].x;
				yend = dl * ey + AP[i - 1].y;
			}
			dx = xend / (Xnum - 1);
			for (j = 0; j < Xnum; j++)
			{
				AP.push_back(t);
				if (j == 0)
				{
					if (i == 0)
						AP[i].x = 0, AP[i].y = 0;
					else
						AP[i].x = 0, AP[i].y = yend;
				}
				else
				{
					AP[i].x = AP[i - 1].x + dx, AP[i].y = AP[i - 1].y;
				}
				i++;
			}
		}
		//左边网格数据
		int k = 0;
		int s = 0;
		while (k < (Xnum - 1) * (Ynum - 1))
		{
			for (j = 0; j < Xnum - 1; j++)
			{
				ad.push_back(a);
				ad[k].push_back(k + s);
				ad[k].push_back(k + 1 + s);
				ad[k].push_back(k + 1 + Xnum + s);
				ad[k].push_back(k + Xnum + s);
				if (j == Xnum - 2)
					s++;
				k++;
			}
		}
		//右边结构网格生成
		ex = cos(mu2 - delta2);
		ey = sin(mu2 - delta2);
		double xstart, ystart;
		double dl2 = 1.5 * LineLength * sin(mu2) / Ynum;
		s = 0;
		while (i < 2 * Pnum - 2)
		{
			if (i == Pnum)
			{
				xstart = LineLength;
				ystart = 0;
				xend = LineLength + 2 * LineLength * cos(delta2);
				yend = -2 * LineLength * sin(delta2);
			}
			else
			{
				if (s == 1)
				{
					xstart = 1.3 * dl * ex + AP[Xnum - 1].x;
					ystart = 1.3 * dl * ey + AP[Xnum - 1].y;
				}
				else
				{
					xstart = 1.3 * dl * ex + AP[i - Xnum].x;
					ystart = 1.3 * dl * ey + AP[i - Xnum].y;
				}
				xend = dl2 * sin(delta2) + AP[i - 1].x;
				yend = dl2 * cos(delta2) + AP[i - 1].y;
			}
			double d = sqrt((xend - xstart) * (xend - xstart) + (yend - ystart) * (yend - ystart)) / (Xnum - 1);
			double ex1 = (xend - xstart) / sqrt((xend - xstart) * (xend - xstart) + (yend - ystart) * (yend - ystart));
			double ey1 = (yend - ystart) / sqrt((xend - xstart) * (xend - xstart) + (yend - ystart) * (yend - ystart));
			for (j = 0; j < Xnum; j++)
			{
				if (j == 0 && i == Pnum)
					continue;
				AP.push_back(t);
				if (j == 0)
					AP[i].x = xstart, AP[i].y = ystart;
				else
				{
					if (i == Pnum)
						AP[i].x = AP[Xnum - 1].x + d * ex1, AP[i].y = AP[Xnum - 1].y + d * ey1;
					else
						AP[i].x = AP[i - 1].x + d * ex1, AP[i].y = AP[i - 1].y + d * ey1;
				}
				i++;
				if (i >= 2 * Pnum - 1)
					break;
			}
			s++;
		}
		//右边网格数据
		k = ad.size();
		s = Xnum + Ynum - 1;
		while (k < 2 * (Xnum - 1) * (Ynum - 1))
		{
			for (j = 0; j < Xnum - 1; j++)
			{
				ad.push_back(a);
				if (k == (Xnum - 1) * (Ynum - 1))
					ad[k].push_back(Xnum - 1);
				else
					ad[k].push_back(k + s - 1);
				ad[k].push_back(k + s);
				ad[k].push_back(k + Xnum + s);
				ad[k].push_back(k + Xnum + s - 1);
				if (j == Xnum - 2)
					s++;
				k++;
			}
		}
		//中间夹角非结构网格
		for (j = 0; j < Ynum; j++)
		{
			if (j == 0)
				continue;
			if (j == 1)
			{
				ad.push_back(a);
				ad[ad.size() - 1].push_back(2 * Xnum - 1);
				ad[ad.size() - 1].push_back(Xnum - 1);
				ad[ad.size() - 1].push_back(Xnum - 1 + Pnum);
			}
			else
			{
				//本层的起始点和结束点
				xstart = AP[(j + 1) * Xnum - 1].x;
				ystart = AP[(j + 1) * Xnum - 1].y;
				xend = AP[j * Xnum - 1 + Pnum].x;
				yend = AP[j * Xnum - 1 + Pnum].y;
				//起始点到结束点的等分长度，以及单位方向向量的xy值
				double d = sqrt((xend - xstart) * (xend - xstart) + (yend - ystart) * (yend - ystart)) / j;
				double ex1 = (xend - xstart) / sqrt((xend - xstart) * (xend - xstart) + (yend - ystart) * (yend - ystart));
				double ey1 = (yend - ystart) / sqrt((xend - xstart) * (xend - xstart) + (yend - ystart) * (yend - ystart));
				int size = AP.size();
				for (k = 0; k < j - 1; k++)
				{
					AP.push_back(t);
					AP[AP.size() - 1].x = xstart + (k + 1) * d * ex1;
					AP[AP.size() - 1].y = ystart + (k + 1) * d * ey1;
					if (j == 2)
					{
						ad.push_back(a);
						ad[ad.size() - 1].push_back((j + 1) * Xnum - 1);
						ad[ad.size() - 1].push_back(j * Xnum - 1);
						ad[ad.size() - 1].push_back(AP.size() - 1);
						ad.push_back(a);
						ad[ad.size() - 1].push_back(AP.size() - 1);
						ad[ad.size() - 1].push_back(j * Xnum - 1);
						ad[ad.size() - 1].push_back((j - 1) * Xnum - 1 + Pnum);
						ad.push_back(a);
						ad[ad.size() - 1].push_back(AP.size() - 1);
						ad[ad.size() - 1].push_back((j - 1) * Xnum - 1 + Pnum);
						ad[ad.size() - 1].push_back(j * Xnum - 1 + Pnum);
					}
					else if (k == 0)
					{
						ad.push_back(a);
						ad[ad.size() - 1].push_back((j + 1) * Xnum - 1);
						ad[ad.size() - 1].push_back(j * Xnum - 1);
						ad[ad.size() - 1].push_back(AP.size() - 1);
						ad.push_back(a);
						ad[ad.size() - 1].push_back(AP.size() - 1);
						ad[ad.size() - 1].push_back(j * Xnum - 1);
						ad[ad.size() - 1].push_back(size - j + 2);
					}
					else if (k == j - 2)
					{
						ad.push_back(a);
						ad[ad.size() - 1].push_back(AP.size() - 1);
						ad[ad.size() - 1].push_back((j - 1) * Xnum - 1 + Pnum);
						ad[ad.size() - 1].push_back(j * Xnum - 1 + Pnum);
						ad.push_back(a);
						ad[ad.size() - 1].push_back(AP.size() - 1);
						ad[ad.size() - 1].push_back(AP.size() - 2);
						ad[ad.size() - 1].push_back(size - j + 2 + k - 1);
						ad.push_back(a);
						ad[ad.size() - 1].push_back(AP.size() - 1);
						ad[ad.size() - 1].push_back(size - j + 2 + k - 1);
						ad[ad.size() - 1].push_back((j - 1) * Xnum - 1 + Pnum);

					}
					else
					{
						ad.push_back(a);
						ad[ad.size() - 1].push_back(AP.size() - 1);
						ad[ad.size() - 1].push_back(AP.size() - 2);
						ad[ad.size() - 1].push_back(size - j + 2 + k - 1);
						ad.push_back(a);
						ad[ad.size() - 1].push_back(AP.size() - 1);
						ad[ad.size() - 1].push_back(size - j + 2 + k - 1);
						ad[ad.size() - 1].push_back(size - j + 2 + k);
					}
				}
			}

		}

	}
}

void reorderMesh()//将边界变为规则
{
	extern vector <mesh> AP;
	extern vector<vector<int>> ad;
	using ConstPara::pi;
	int i = 0;
	if ((FlowType == "oblique" || FlowType == "intersection" || FlowType == "normal" || FlowType == "couette") && meshType != 20)
	{
		for (i = 0; i < AP.size(); i++)
		{
			if (abs(AP[i].x - xL) <= 1e-10 || abs(AP[i].x - xR) <= 1e-10)
				continue;
			if (abs(AP[i].x - (xL + dx * cos(pi / 3))) <= 1e-10)
				AP[i].x = xL, AP[i].type = "L";
			else if (abs(AP[i].x - (xR - dx * cos(pi / 3))) <= 1e-10)
				AP[i].x = xR, AP[i].type = "R";
			else if (abs(AP[i].y - (yD + dx * sin(pi / 3))) <= 1e-10)
				AP[i].y = yD, AP[i].type = "D";
			else if (abs(AP[i].y - (yU - dx * sin(pi / 3))) <= 1e-10)
				AP[i].y = yU, AP[i].type = "U";
		}
	}
}
void remesh_bound()//将边界x或y坐标移动到和内部点相同，改善边界条件
{
	extern vector<vector <mesh*>> A;
	int i, j, k;
	mesh *m;
	if (FlowType == "oblique" || FlowType == "intersection" || FlowType == "normal")
	{
		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j]->type == "D" || A[i][j]->type == "U")
				{
					for (k = 0; k < A[i][j]->neibor.size(); k++)
					{
						m = A[i][j]->neibor[k];
						if (abs(A[i][j]->y - m->y) < 1e-10)
							continue;
						else
							A[i][j]->x = m->x;
					}
				}
			}
		}
	}
}
void init_mesh1()
{
	extern vector <mesh> AP;
	int i, j;
	i = 0;
	mesh t;
	while (i < Pnum)
	{
		for (j = 0; j < Xnum; j++)
		{
			AP.push_back(t);

			if (j == 0)
			{
				if (i == 0)
					AP[i].x = 0, AP[i].y = 0;
				else if (AP[i - Xnum].x == 0)
					AP[i].x = 0, AP[i].y = AP[i - Xnum].y + dy;
				else
					AP[i].x = 0, AP[i].y = AP[i - Xnum].y + dy;
			}
			else if (AP[i - j].x == 0)
			{
				if (j % 2 == 0)
					AP[i].x = AP[i - 1].x + dx, AP[i].y = AP[i - 1].y;
				else
					AP[i].x = AP[i - 1].x + dx, AP[i].y = AP[i - 1].y;
			}
			else
			{
				if (j % 2 == 0)
					AP[i].x = AP[i - 1].x + dx, AP[i].y = AP[i - 1].y;
				else
					AP[i].x = AP[i - 1].x + dx, AP[i].y = AP[i - 1].y;
			}
			AP[i].x0 = AP[i].x;
			AP[i].y0 = AP[i].y;
			i++;
		}
	}
}
void init_polygon_mesh()
{

	extern vector <mesh> AP;
	extern vector <polygon_mesh> PM;
	polygon_mesh t;
	int i, j;
	i = 0;
	int k = 0;
	if (meshType == 10 || meshType == 11 || meshType == 12)
	{
		while (i < Pnum)
		{

			if (i + 1 + 2 * Xnum > Pnum - 1)
				break;

			k++;
			for (int j = 0; j < Xnum; j++)
			{
				if (abs(abs(AP[i].x - AP[i + 1].x) - dx) > 1e-10)
				{
					i++;
					continue;
				}
				if ((i + k - 1) % 2 != 0)
				{
					i++;
					continue;
				}
				if (i + 1 + 2 * Xnum > Pnum - 1)
					continue;
				PM.push_back(t);
				int s = int(PM.size()) - 1;
				PM[s].node.push_back(i);
				PM[s].node.push_back(i + 1);
				PM[s].node.push_back(i + 1 + Xnum);
				PM[s].node.push_back(i + 1 + 2 * Xnum);
				PM[s].node.push_back(i + 2 * Xnum);
				PM[s].node.push_back(i + Xnum);

				PM[s].face_start.push_back(i);
				PM[s].face_start.push_back(i + 1);
				PM[s].face_start.push_back(i + 1 + Xnum);
				PM[s].face_start.push_back(i + 1 + 2 * Xnum);
				PM[s].face_start.push_back(i + 2 * Xnum);
				PM[s].face_start.push_back(i + Xnum);

				PM[s].face_end.push_back(i + 1);
				PM[s].face_end.push_back(i + 1 + Xnum);
				PM[s].face_end.push_back(i + 1 + 2 * Xnum);
				PM[s].face_end.push_back(i + 2 * Xnum);
				PM[s].face_end.push_back(i + Xnum);
				PM[s].face_end.push_back(i);
				i++;
			}
		}
	}
	else
	{
		int s;
		while (i < Pnum)
		{
			if (abs(AP[i].x - xR) < 1e-10 || abs(AP[i].y - yU) < 1e-10)
			{
				i++;
				continue;
			}
			PM.push_back(t);
			s = int(PM.size()) - 1;
			PM[s].node.push_back(i);
			PM[s].node.push_back(i + 1);
			PM[s].node.push_back(i + 1 + Xnum);
			PM[s].node.push_back(i + Xnum);

			PM[s].face_start.push_back(i);
			PM[s].face_start.push_back(i + 1);
			PM[s].face_start.push_back(i + 1 + Xnum);
			PM[s].face_start.push_back(i + Xnum);

			PM[s].face_end.push_back(i + 1);
			PM[s].face_end.push_back(i + 1 + Xnum);
			PM[s].face_end.push_back(i + Xnum);
			PM[s].face_end.push_back(i);
			i++;
		}

	}
}
void getType()
{
	extern vector <mesh> AP;
	if (FlowType == "oblique" || FlowType == "intersection" || FlowType == "normal" || FlowType == "couette" || FlowType == "cylinder")
	{
		int i;
		float DELTA = 1e-10;
		for (i = 0; i < Pnum; i++)
		{
			if (abs(AP[i].x - xL) < DELTA)
			{
				AP[i].type = "L";
			}
			else if (abs(AP[i].x - xR) < DELTA)
			{
				AP[i].type = "R";

			}
			else if (abs(AP[i].y - yD) < DELTA)
			{
				AP[i].type = "D";
			}
			else if (abs(AP[i].y - yU) < DELTA)
			{
				AP[i].type = "U";
			}
			else
			{
				AP[i].type = "IN";
			}
		}

	}
	if (FlowType == "Prandtl-Meyer")
	{
		int i;
		for (i = 0; i < Pnum; i++)
		{
			if (AP[i].x == 0)
				AP[i].type = "L";
			else if (AP[i].y == 0)
				AP[i].type = "DL";//下边界的左半边
			else
				AP[i].type = "IN";

		}
		for (i = Pnum; i < 2 * Pnum - 1; i++)
		{
			if (i < Pnum + Xnum - 1)
				AP[i].type = "DR";//下边界的右半边
			else if ((i - (Pnum + Xnum - 2)) % Xnum == 0)
				AP[i].type = "R";
			else
				AP[i].type = "IN";
		}
		for (i = 2 * Pnum - 1; i < AP.size(); i++)
		{
			AP[i].type = "IN";
		}
		AP[Xnum - 1].type = "IN";
	}
}

void initFlow()
{
	if (methodType == "F")
	{
		int i, j;
		extern vector<vector <mesh*>> A;
		if (FlowType == "normal")
		{
			using namespace Normal;
			mesh A1, A2;
			A1.rho = rho1, A1.p = p1, A1.u = u1, A1.v = v1;
			get_down(A1, A2, ConstPara::pi / 2);
			rho2 = A2.rho, p2 = A2.p, u2 = A2.u, v2 = A2.v;
			u2 += 1;
			u1 += 1;
			for (i = 0; i < A.size(); i++)
			{
				for (j = 0; j < A[i].size(); j++)
				{
					if (i == 0)
					{
						A[i][j]->rho = rho1;
						A[i][j]->u = u1;
						A[i][j]->v = v1;
						A[i][j]->p = p1;
					}
					else if (i == 1)
					{
						A[i][j]->rho = rho2;
						A[i][j]->u = u2;
						A[i][j]->v = v2;
						A[i][j]->p = p2;
					}

				}

			}
		}
		else if (FlowType == "oblique")
		{
			using namespace Oblique;
			mesh A1, A2;
			A1.rho = rho1, A1.p = p1, A1.u = u1, A1.v = v1;
			get_down(A1, A2, beta);
			rho2 = A2.rho, p2 = A2.p, u2 = A2.u, v2 = A2.v;

			for (i = 0; i < A.size(); i++)
			{
				for (j = 0; j < A[i].size(); j++)
				{
					if (i == 0)
					{
						A[i][j]->rho = rho1;
						A[i][j]->u = u1;
						A[i][j]->v = v1;
						A[i][j]->p = p1;
					}
					else if (i == 1)
					{
						A[i][j]->rho = rho2;
						A[i][j]->u = u2;
						A[i][j]->v = v2;
						A[i][j]->p = p2;
					}

				}

			}
		}
		else if (FlowType == "intersection")
		{
			using namespace ShockwaveCross;
			using namespace ConstPara;
			double beta40 = 45 * pi / 180, beta41 = 50 * pi / 180, beta42 = 60 * pi / 180;
			double p41, p40;
			double p51, p50;
			double delta41, delta40;
			double delta51, delta50;
			double beta51, beta50;
			double fbeta41, fbeta40;
			double un20, un21, Mu20, Mu21;
			double un30, un31, Mu30, Mu31;
			mesh A1, A2, A3;
			A1.rho = rho1, A1.p = p1, A1.u = u1, A1.v = v1;
			get_down(A1, A2, beta2/*-0.05*pi/180*/);
			get_down(A1, A3, beta3);

			rho2 = A2.rho, p2 = A2.p, u2 = A2.u, v2 = A2.v;
			rho3 = A3.rho, p3 = A3.p, u3 = A3.u, v3 = A3.v;

			//A2.rho = rho2, A2.p = p2, A2.u = u2, A2.v = v2;
			//A3.rho = rho3, A3.p = p3, A3.u = u3, A3.v = v3;
			mesh A40, A41, A50, A51;
			while (abs(beta42 - beta41) > 1e-15)
			{
				while (beta42 > -beta2 || beta42 < beta2)
				{
					if (abs(beta42) > 1e5)
					{
						beta42 = beta42 / 1e5;
						continue;
					}
					if (beta42 > -beta2)
						beta42 = beta42 + beta2;
					else if (beta42 < beta2)
						beta42 = beta42 - beta2;
					//if (beta42 > pi / 2)
					//	beta42 = beta42 - pi / 2;
					//	else if (beta42 < -pi / 2)
					//	beta42 = beta42 + pi / 2;
				}
				if (abs(beta42) < 1e-7)
					beta42 = 30 * pi / 180;

				beta40 = (beta41 + beta42) / 2;
				beta41 = beta42;

				get_down(A2, A40, beta40);
				A50.p = A40.p;
				beta50 = get_beta(A3, A50.p, -1);
				get_down(A3, A50, beta50);

				get_down(A2, A41, beta41);
				A51.p = A41.p;
				beta51 = get_beta(A3, A51.p, -1);
				get_down(A3, A51, beta51);

				delta40 = get_delta(A40.u, A40.v);
				delta41 = get_delta(A41.u, A41.v);
				delta50 = get_delta(A50.u, A50.v);
				delta51 = get_delta(A51.u, A51.v);
				fbeta41 = delta41 - delta51;
				fbeta40 = delta40 - delta50;
				if (abs(fbeta41 - fbeta40) <= 1e-20)
				{
					beta42 = beta41;
					break;
				}
				beta42 = beta41 - fbeta41 / (fbeta41 - fbeta40) * (beta41 - beta40);
			}
			beta4 = beta42;
			get_down(A2, A41, beta4);
			A51.p = A41.p;
			beta5 = get_beta(A3, A51.p, -1);
			get_down(A3, A51, beta5);
			delta4 = get_delta(A41.u, A41.v);
			delta5 = get_delta(A51.u, A51.v);

			double un2 = u2 * sin(beta4) - v2 * cos(beta4);
			double ut2 = u2 * cos(beta4) + v2 * sin(beta4);
			double Mu2 = get_Ma(un2, 0, rho2, p2);
			rho4 = rho2 * get_rho2rho1(Mu2);
			p4 = p2 * get_p2p1(Mu2);
			double Md4 = sqrt((Mu2 * Mu2 + 2 / (gama - 1)) / (2 * gama * Mu2 * Mu2 / (gama - 1) - 1));
			double c4 = sqrt(gama * p4 / rho4);
			double un4 = Md4 * c4;
			double ut4 = ut2;
			u4 = un4 * sin(beta4) + ut4 * cos(beta4);
			v4 = -un4 * cos(beta4) + ut4 * sin(beta4);


			double un3 = -u3 * sin(beta5) + v3 * cos(beta5);
			double ut3 = u3 * cos(beta5) + v3 * sin(beta5);
			double Mu3 = get_Ma(un3, 0, rho3, p3);
			rho5 = rho3 * get_rho2rho1(Mu3);
			p5 = p3 * get_p2p1(Mu3);
			double Md5 = sqrt((Mu3 * Mu3 + 2 / (gama - 1)) / (2 * gama * Mu3 * Mu3 / (gama - 1) - 1));
			double c5 = sqrt(gama * p5 / rho5);
			double un5 = Md5 * c5;
			double ut5 = ut3;
			u5 = un5 * sin(-beta5) + ut5 * cos(-beta5);
			v5 = un5 * cos(-beta5) - ut5 * sin(-beta5);
			std::cout.precision(20);
			std::cout << "p2= " << p2 << "   p3= " << p3 << std::endl;
			std::cout << "u2= " << u2 << "   u3= " << u3 << std::endl;
			std::cout << "v2= " << v2 << "   v3= " << v3 << std::endl;
			std::cout << "rho2= " << rho2 << "   rho3= " << rho3 << std::endl;
			std::cout << "beta4= " << beta4 * 180 / pi << "   beta5= " << beta5 * 180 / pi << std::endl;
			std::cout << "delta4= " << delta4 * 180 / pi << "   delta5= " << delta5 * 180 / pi << std::endl;
			std::cout << "p4= " << p4 << "   p5= " << p5 << std::endl;
			std::cout << "u4= " << u4 << "   u5= " << u5 << std::endl;
			std::cout << "v4= " << v4 << "   v5= " << v5 << std::endl;
			std::cout << "rho4= " << rho4 << "   rho5= " << rho5 << std::endl;
			std::cout << std::endl;
			for (i = 0; i < A.size(); i++)
			{
				for (j = 0; j < A[i].size(); j++)
				{
					A[i][j]->um = A[i][j]->vm = 0;
					if (i == 0)
					{
						A[i][j]->rho = rho1;
						A[i][j]->u = u1;
						A[i][j]->v = v1;
						A[i][j]->p = p1;
					}
					else if (i == 1)
					{
						A[i][j]->rho = rho2;
						A[i][j]->u = u2;
						A[i][j]->v = v2;
						A[i][j]->p = p2;
					}
					else if (i == 2)
					{
						//if (A[i][j]->type == "L")
						//{
						A[i][j]->rho = rho3;
						A[i][j]->u = u3;
						A[i][j]->v = v3;
						A[i][j]->p = p3;

						//}
						//else
						//{
						//	A[i][j]->rho = rho3;
						//	A[i][j]->u = u3 - 5;
						//	A[i][j]->v = v3 -0.5;
						//	A[i][j]->p = p3;
						//}
					}
					else if (i == 3)
					{
						A[i][j]->rho = rho4;
						A[i][j]->u = u4;
						A[i][j]->v = v4;
						A[i][j]->p = p4;
					}
					else if (i == 4)
					{
						A[i][j]->rho = rho5;
						A[i][j]->u = u5;
						A[i][j]->v = v5;
						A[i][j]->p = p5;
					}
					else
					{
						std::cout << "something wrong in intersection!" << std::endl;
					}

				}

			}

		}
	}
	else if (methodType == "C")
	{
		if (FlowType == "normal")
			init_flow_normal();
		else if (FlowType == "cylinder")
			init_flow_cylinder();//圆柱绕流
		else if (FlowType == "oblique")
			init_flow_shockwave();
		else if (FlowType == "intersection")
			init_flow_shockwaveCross();
		else if (FlowType == "couette")
		{
			using namespace Couette;
			mesh A1, A2;
			A1.rho = rho1, A1.p = p1, A1.u = u1, A1.v = v1;
			int i, j;
			extern vector<vector <mesh*>> A;
			for (i = 0; i < A.size(); i++)
			{
				for (j = 0; j < A[i].size(); j++)
				{
					A[i][j]->rho = rho1;
					A[i][j]->u = 0;
					A[i][j]->v = v1;
					A[i][j]->p = p1;
				}

			}

		}

	}
}
void init_flow_normal()
{
	extern vector<vector <mesh*>> A;
	using namespace Normal;
	int i, j;
	mesh A1, A2;
	A1.rho = rho1, A1.p = p1, A1.u = u1, A1.v = v1;
	get_down(A1, A2, ConstPara::pi / 2);
	rho2 = A2.rho, p2 = A2.p, u2 = A2.u, v2 = A2.v;
	u2 += 1;
	u1 += 1;
	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			if (A[i][j]->x < 0.002)
			{
				A[i][j]->rho = rho1;
				A[i][j]->u = u1;
				A[i][j]->v = v1;
				A[i][j]->p = p1;
			}
			else
			{
				A[i][j]->rho = rho2;
				A[i][j]->u = u2;
				A[i][j]->v = v2;
				A[i][j]->p = p2;
			}
		}
	}

}
void init_flow_uniform()//均匀流
{
	extern vector<vector <mesh*>> A;
	int i;
	using namespace Init;
	for (i = 0; i < A[0].size(); i++)
	{
		A[0][i]->rho = rho0;
		A[0][i]->u = u0;
		A[0][i]->v = v0;
		A[0][i]->p = p0;
	}
}
void init_flow_shockwave()//斜激波
{
	extern vector<vector <mesh*>> A;
	using namespace Oblique;
	mesh A1, A2;
	A1.rho = rho1, A1.p = p1, A1.u = u1, A1.v = v1;
	get_down(A1, A2, beta);
	rho2 = A2.rho, p2 = A2.p, u2 = A2.u, v2 = A2.v;
	using std::cout;
	using std::endl;
	cout.precision(10);
	cout << "rho1=" << rho1 << ", u1= " << u1 << ", v1=" << v1 << ",  p1=" << p1 << endl;
	cout << "rho2=" << rho2 << ", u2= " << u2 << ", v2=" << v2 << ",  p2=" << p2 << endl;

	int i, j;
	for (i = 0; i < A.size(); i++)
	{
		for (j = 0; j < A[i].size(); j++)
		{
			if (A[i][j]->x < 0.0075)
			{
				A[i][j]->rho = rho1;
				A[i][j]->u = u1;
				A[i][j]->v = v1;
				A[i][j]->p = p1;
			}
			else
			{
				A[i][j]->rho = rho2;
				A[i][j]->u = u2;
				A[i][j]->v = v2;
				A[i][j]->p = p2;
			}
		}
		//if (A[0][i].type == "IN")
		//{
		//	A[0][i].rho = rho1;
		//	A[0][i].u = u1;
		//	A[0][i].v = v1;
		//	A[0][i].p = p1;
		//}
		//else if (A[0][i].type == "L")
		//{
		//	A[0][i].rho = rho1;
		//	A[0][i].u = u1;
		//	A[0][i].v = v1;
		//	A[0][i].p = p1;
		//}
		//else if (A[0][i].type == "U")
		//{
		//	if (A[0][i].x < A[0][startpoint].x)
		//	{
		//		A[0][i].rho = rho1;
		//		A[0][i].u = u1;
		//		A[0][i].v = v1;
		//		A[0][i].p = p1;
		//	}
		//	else
		//	{
		//		A[0][i].rho = rho2;
		//		A[0][i].u = u2;
		//		A[0][i].v = v2;
		//		A[0][i].p = p2;
		//	}
		//}
		//else if (A[0][i].type == "D" || A[0][i].type == "R")
		//{
		//	if (A[0][i].x < A[0][startpoint].x + dy * Ynum / tan(beta))
		//	{
		//		A[0][i].rho = rho1;
		//		A[0][i].u = u1;
		//		A[0][i].v = v1;
		//		A[0][i].p = p1;
		//	}
		//	else
		//	{
		//		A[0][i].rho = rho2;
		//		A[0][i].u = u2;
		//		A[0][i].v = v2;
		//		A[0][i].p = p2;
		//	}
		//}
		//else
		//	std::cout << "something wrong!" << std::endl;
	}
}
void init_flow_shockwaveCross()//同侧激波相交
{
	extern vector<vector <mesh*>> A;
	extern vector <mesh> AP;
	double Ma1;
	using namespace ShockwaveCross;
	using namespace Init;
	int i;
	using namespace ConstPara;
	Line L12;
	Line L13;
	mesh C;

	mesh A1, A2, A3;
	A1.rho = rho1, A1.p = p1, A1.u = u1, A1.v = v1;
	get_down(A1, A2, beta2);
	get_down(A1, A3, beta3);
	rho2 = A2.rho, p2 = A2.p, u2 = A2.u, v2 = A2.v;
	rho3 = A3.rho, p3 = A3.p, u3 = A3.u, v3 = A3.v;
	double theta12 = -30 / 180.0 * ConstPara::pi;
	double theta13 = 30 / 180.0 * ConstPara::pi;
	mesh A12;
	mesh A13;
	vector<mesh> t;
	A12.x = -dx / 2, A12.y = AP[Pnum - 1].y - 2 * dy * 5 * Ynum / 43;
	A13.x = -dx / 2, A13.y = 2 * dy * 5 * Ynum / 43;
	L12 = getLine(theta12, A12);
	L13 = getLine(theta13, A13);
	C = getCrossPoint(L12, L13);
	L12 = getLine(beta2, C);
	L13 = getLine(beta3, C);
	for (i = 0; i < Pnum; i++)
	{
		if (A[0][i]->type == "L")
		{
			if (A[0][i]->x * L12.A + A[0][i]->y * L12.B + L12.C < 0)
			{
				A[0][i]->rho = rho2;
				A[0][i]->u = u2;
				A[0][i]->v = v2;
				A[0][i]->p = p2;
			}
			else if (A[0][i]->x * L13.A + A[0][i]->y * L13.B + L13.C > 0)
			{
				A[0][i]->rho = rho3;
				A[0][i]->u = u3;
				A[0][i]->v = v3;
				A[0][i]->p = p3;
			}
			else
			{
				A[0][i]->rho = rho1;
				A[0][i]->u = u1;
				A[0][i]->v = v1;
				A[0][i]->p = p1;
			}
		}
		else
		{
			A[0][i]->rho = rho0;
			A[0][i]->u = u0;
			A[0][i]->v = v0;
			A[0][i]->p = p0;
		}
	}
}
void init_flow_cylinder()//圆柱绕流
{
	extern vector<vector <mesh*>> A;
	extern vector <mesh> AP;
	double Ma1;
	using namespace Init;
	int i;
	using namespace ConstPara;
	for (i = 0; i < AP.size(); i++)
	{
		//if (AP[i].section == 0)
		//	continue;
		AP[i].rho = rho0;
		AP[i].u = u0;
		AP[i].v = v0;
		AP[i].p = p0;
	}
}

//void init_U()
//{
//	extern vector<vector <mesh*>> A;
//	extern vector<vector<vector <double>>> U;
//	vector<vector <double>> u0;
//	vector <double> u1;
//	int i, j;
//	for (i = 0; i < A.size(); i++)
//	{
//
//		U.push_back(u0);
//		for (j = 0; j < A[i].size(); j++)
//		{
//			U[i].push_back(u1);
//			U[i][j].push_back(A[i][j]->rho);
//			U[i][j].push_back(A[i][j]->rho * A[i][j]->u);
//			U[i][j].push_back(A[i][j]->rho * A[i][j]->v);
//			U[i][j].push_back(0.5 * A[i][j]->rho * (A[i][j]->u * A[i][j]->u + A[i][j]->v * A[i][j]->v) + A[i][j]->p / (ConstPara::gama - 1));
//		}
//	}
//}
void findConnectPoint()
{
	extern vector<vector <mesh*>> A;
	int i, j, k;
	if (FlowType == "oblique" || FlowType == "normal")
	{
		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j]->type == "SHOCK")
					continue;
				for (k = 0; k < A[i].size(); k++)
				{
					if (A[i][k]->type != "SHOCK")
						continue;
					if (abs(A[i][j]->y - A[i][k]->y) <= dy + 1e-10)
						A[i][j]->moveConnct.push_back(k);
				}

			}
		}
	}
	else if (FlowType == "intersection")
	{

		for (i = 0; i < A.size(); i++)
		{
			for (j = 0; j < A[i].size(); j++)
			{
				if (A[i][j]->type == "SHOCK" || A[i][j]->type == "DISCON" || A[i][j]->type == "CENTER")
					continue;
				for (k = 0; k < A[i].size(); k++)
				{
					if (A[i][k]->type == "SHOCK" || A[i][k]->type == "DISCON")
					{
						if (abs(A[i][k]->x - A[i][j]->x) <= dx / 2 + 1e-10)
							A[i][j]->moveConnct.push_back(k);
					}

				}

			}
		}



	}
}
//void init_Ar()
//{
//	extern vector <mesh> AP;
//	extern vector <mesh> Ar;
//	for (int i = 0; i < Pnum; i++)
//	{
//		Ar[i].x = AP[i].x;
//		Ar[i].y = AP[i].y;
//		Ar[i].rho = AP[i].rho;
//		Ar[i].u = AP[i].u;
//		Ar[i].v = AP[i].v;
//		Ar[i].p = AP[i].p;
//		Ar[i].neibor = AP[i].neibor;
//	}
//}