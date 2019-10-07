#include"include/const.h"
#include"include/functions.h"
#include<vector>
using std::vector;
void coordinate_trans()
{
	extern vector<mesh> AP;
	int i, j, k;
	int n1, n2, n3, n4;
	extern double t_sim;
	extern double dt;
	extern  std::vector <mesh> Ar;

	for (i = 0; i < AP.size(); i++)
	{

		if (AP[i].neibor.size() < 3)//小于3的点无法做坐标变换
			continue;
		if (t_sim == 0)
		{
			if (AP[i].neibor.size() == 3)
			{
				using MeshPara::method;
				for (k = 0; k < 12; k++)
				{
					n1 = method[k][0];
					n2 = method[k][1];
					n3 = method[k][2];
					n4 = method[k][3];
					AP[i].xxi.push_back(0.5 * (AP[i].neibor[n1]->x - AP[i].neibor[n3]->x));
					AP[i].yxi.push_back(0.5 * (AP[i].neibor[n1]->y - AP[i].neibor[n3]->y));
					AP[i].xeta.push_back(0.5 * (AP[i].neibor[n2]->x - AP[i].neibor[n4]->x));
					AP[i].yeta.push_back(0.5 * (AP[i].neibor[n2]->y - AP[i].neibor[n4]->y));
					AP[i].J.push_back(1 / (AP[i].xxi[k] * AP[i].yeta[k] - AP[i].xeta[k] * AP[i].yxi[k]));
					AP[i].xix.push_back(AP[i].yeta[k] * AP[i].J[k]);
					AP[i].xiy.push_back(-AP[i].xeta[k] * AP[i].J[k]);
					AP[i].etax.push_back(-AP[i].yxi[k] * AP[i].J[k]);
					AP[i].etay.push_back(AP[i].xxi[k] * AP[i].J[k]);
					AP[i].xtau.push_back(0);
					AP[i].ytau.push_back(0);
					AP[i].xit.push_back(0);
					AP[i].etat.push_back(0);
				}
			}
			else if (AP[i].neibor.size() >= 4)
			{
				AP[i].xxi.push_back(0.5 * (AP[i].neibor[0]->x - AP[i].neibor[2]->x));
				AP[i].yxi.push_back(0.5 * (AP[i].neibor[0]->y - AP[i].neibor[2]->y));
				AP[i].xeta.push_back(0.5 * (AP[i].neibor[1]->x - AP[i].neibor[3]->x));
				AP[i].yeta.push_back(0.5 * (AP[i].neibor[1]->y - AP[i].neibor[3]->y));
				AP[i].J.push_back(1 / (AP[i].xxi[0] * AP[i].yeta[0] - AP[i].xeta[0] * AP[i].yxi[0]));
				AP[i].xix.push_back(AP[i].yeta[0] * AP[i].J[0]);
				AP[i].xiy.push_back(-AP[i].xeta[0] * AP[i].J[0]);
				AP[i].etax.push_back(-AP[i].yxi[0] * AP[i].J[0]);
				AP[i].etay.push_back(AP[i].xxi[0] * AP[i].J[0]);
				AP[i].xtau.push_back(0);
				AP[i].ytau.push_back(0);
				AP[i].xit.push_back(0);
				AP[i].etat.push_back(0);
			}
		}
		else//t_sim不等于零时还需要坐标变换，应为动网格！
		{
			if (AP[i].neibor.size() >= 4)
			{
				AP[i].xxi[0] = 0.5 * (AP[i].neibor[0]->x - AP[i].neibor[2]->x);
				AP[i].yxi[0] = 0.5 * (AP[i].neibor[0]->y - AP[i].neibor[2]->y);
				AP[i].xeta[0] = 0.5 * (AP[i].neibor[1]->x - AP[i].neibor[3]->x);
				AP[i].yeta[0] = 0.5 * (AP[i].neibor[1]->y - AP[i].neibor[3]->y);
				AP[i].J[0] = 1 / (AP[i].xxi[0] * AP[i].yeta[0] - AP[i].xeta[0] * AP[i].yxi[0]);
				AP[i].xix[0] = AP[i].yeta[0] * AP[i].J[0];
				AP[i].xiy[0] = -AP[i].xeta[0] * AP[i].J[0];
				AP[i].etax[0] = -AP[i].yxi[0] * AP[i].J[0];
				AP[i].etay[0] = AP[i].xxi[0] * AP[i].J[0];
				AP[i].xtau[0] = (AP[i].x - Ar[AP[i].id].x) / dt;
				AP[i].ytau[0] = (AP[i].y - Ar[AP[i].id].y) / dt;
				AP[i].xit[0] = -AP[i].xtau[0] * AP[i].xix[0] - AP[i].ytau[0] * AP[i].xiy[0];
				AP[i].etat[0] = -AP[i].xtau[0] * AP[i].etax[0] - AP[i].ytau[0] * AP[i].etay[0];
			}
			if (AP[i].neibor.size() == 3)
			{
				AP[i].xxi[0] = 0.5 * (AP[i].neibor[0]->x - AP[i].neibor[2]->x);
				AP[i].yxi[0] = 0.5 * (AP[i].neibor[0]->y - AP[i].neibor[2]->y);
				AP[i].xeta[0] = 0.5 * (AP[i].neibor[1]->x - AP[i].neibor[2]->x);
				AP[i].yeta[0] = 0.5 * (AP[i].neibor[1]->y - AP[i].neibor[2]->y);
				AP[i].J[0] = 1 / (AP[i].xxi[0] * AP[i].yeta[0] - AP[i].xeta[0] * AP[i].yxi[0]);
				AP[i].xix[0] = AP[i].yeta[0] * AP[i].J[0];
				AP[i].xiy[0] = -AP[i].xeta[0] * AP[i].J[0];
				AP[i].etax[0] = -AP[i].yxi[0] * AP[i].J[0];
				AP[i].etay[0] = AP[i].xxi[0] * AP[i].J[0];
				AP[i].xtau[0] = (AP[i].x - Ar[AP[i].id].x) / dt;
				AP[i].ytau[0] = (AP[i].y - Ar[AP[i].id].y) / dt;
				AP[i].xit[0] = -AP[i].xtau[0] * AP[i].xix[0] - AP[i].ytau[0] * AP[i].xiy[0];
				AP[i].etat[0] = -AP[i].xtau[0] * AP[i].etax[0] - AP[i].ytau[0] * AP[i].etay[0];

				AP[i].xxi[1] = 0.5 * (AP[i].neibor[1]->x - AP[i].neibor[0]->x);
				AP[i].yxi[1] = 0.5 * (AP[i].neibor[1]->y - AP[i].neibor[0]->y);
				AP[i].xeta[1] = 0.5 * (AP[i].neibor[2]->x - AP[i].neibor[0]->x);
				AP[i].yeta[1] = 0.5 * (AP[i].neibor[2]->y - AP[i].neibor[0]->y);
				AP[i].J[1] = 1 / (AP[i].xxi[1] * AP[i].yeta[1] - AP[i].xeta[1] * AP[i].yxi[1]);
				AP[i].xix[1] = AP[i].yeta[1] * AP[i].J[1];
				AP[i].xiy[1] = -AP[i].xeta[1] * AP[i].J[1];
				AP[i].etax[1] = -AP[i].yxi[1] * AP[i].J[1];
				AP[i].etay[1] = AP[i].xxi[1] * AP[i].J[1];
				AP[i].xtau[1] = (AP[i].x - Ar[AP[i].id].x) / dt;
				AP[i].ytau[1] = (AP[i].y - Ar[AP[i].id].y) / dt;
				AP[i].xit[1] = -AP[i].xtau[1] * AP[i].xix[1] - AP[i].ytau[1] * AP[i].xiy[1];
				AP[i].etat[1] = -AP[i].xtau[1] * AP[i].etax[1] - AP[i].ytau[1] * AP[i].etay[1];

				AP[i].xxi[2] = 0.5 * (AP[i].neibor[2]->x - AP[i].neibor[1]->x);
				AP[i].yxi[2] = 0.5 * (AP[i].neibor[2]->y - AP[i].neibor[1]->y);
				AP[i].xeta[2] = 0.5 * (AP[i].neibor[0]->x - AP[i].neibor[1]->x);
				AP[i].yeta[2] = 0.5 * (AP[i].neibor[0]->y - AP[i].neibor[1]->y);
				AP[i].J[2] = 1 / (AP[i].xxi[2] * AP[i].yeta[2] - AP[i].xeta[2] * AP[i].yxi[2]);
				AP[i].xix[2] = AP[i].yeta[2] * AP[i].J[2];
				AP[i].xiy[2] = -AP[i].xeta[2] * AP[i].J[2];
				AP[i].etax[2] = -AP[i].yxi[2] * AP[i].J[2];
				AP[i].etay[2] = AP[i].xxi[2] * AP[i].J[2];
				AP[i].xtau[2] = (AP[i].x - Ar[AP[i].id].x) / dt;
				AP[i].ytau[2] = (AP[i].y - Ar[AP[i].id].y) / dt;
				AP[i].xit[2] = -AP[i].xtau[2] * AP[i].xix[2] - AP[i].ytau[2] * AP[i].xiy[2];
				AP[i].etat[2] = -AP[i].xtau[2] * AP[i].etax[2] - AP[i].ytau[2] * AP[i].etay[2];

			}

		}
	}
}

