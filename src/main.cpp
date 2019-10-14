#include<iostream>
#include<fstream>
#include<cmath>
#include<omp.h>
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
vector<mesh*> ps;//�ṹ����ڵ�
vector<mesh*> pu;//�ǽṹ����ڵ�
vector<mesh*> bl;//��߽�
vector<mesh*> br;//�ұ߽�
vector<mesh*> bu;//�ϱ߽�
vector<mesh*> bd;//�±߽�
vector<mesh*> bb;//����߽�
vector<vector <int>> ad;
vector<mesh> Ar;//��¼��һʱ�̵���������

double dt;
double t_sim = 0;
int step = 0;
vector <polygon_mesh> PM;//������������������������޹�
vector<mesh> poly;//������棬��ɺܶ�ɢ�㣬�����жϽṹ����ڵ㴦����������
double res;
int main()
{
	using namespace MeshPara;
	using ConstPara::t_end;
	//init_mesh1();
	init_mesh();
	polygonPoint(poly);
	//init_polygon_mesh();
	getType();
	partition_Point();
	sortPoint();
	polymesh();
	out_M("mesh/step = " + to_string(step));
	//out_neighbor();
	coordinate_trans();
	initFlow();
	int i;
	while (t_sim < t_end)
	{
		record();
		get_dt();
		for (i = 0; i < ps.size(); i++)
			update_p4_s(*ps[i]);
		for (i = 0; i < pu.size(); i++)
		{
			if (pu[i]->neibor.size() == 3)
				update_p3(*pu[i]);
			else
				update_p4_u(*pu[i]);
			AP[pu[i]->connectId] = *pu[i];//replace
		}
		update_bound();
		if (++step % 100 == 0)
		{
			if (abs(res - compute_res()) < 1e-20)
				break;
			else
				res = compute_res();
			cout << "step = " << step << "  t_sim = " << t_sim << "  dt = " << dt << "  res = " << res << endl;
			out_M("mesh/step = " + to_string(step));
		}
		out_res();
		t_sim = t_sim + dt;
	}
	out_M("mesh/step = " + to_string(step));
	system("PAUSE");
}