#pragma once
#include<cmath>
#include<string>
#include<vector>
using std::vector;
using std::string;
namespace ConstPara//����
{
	const double CFL = 0.3;
	const double t_end = 300;
	const double pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628;
	const double gama = 1.4;
	const double random = 0.495;
}
namespace Init//��ʼ����
{
	extern double rho0;
	extern double v0;
	extern double p0;
	extern double u0;
}
namespace Normal//��ʼ����
{
	extern double rho1;
	extern double v1;
	extern double p1;
	extern double u1;
	extern double rho2;
	extern double v2;
	extern double p2;
	extern double u2;
}
namespace Couette//��ʼ����
{
	extern double rho1;
	extern double v1;
	extern double p1;
	extern double u1;
}

namespace MeshPara
{
	const string FlowType = "cylinder";
	//Ŀǰ�У����Ͳ�Prandtl-Meyer;������normal��б����oblique �������ཻintersection,������couette��Բ������cylinder�� ����other
	const int meshType = 20;
	//10������������11�Ŷ�������������12ש������
	//20�����νṹ����
	const string methodType = "C";
	//���㷽����������׽C-Capturing������װ��F-Fitting
	const int method[12][4] = {
{1,2,0,0},//����3
{0,1,2,2},
{2,0,1,1},

{0,1,1,2},//����2
{1,2,2,0},
{2,0,0,1},


{0,0,1,2},//����1
{1,1,2,0},
{2,2,0,1},

{0,1,2,0},//����4
{2,0,1,2},
{1,2,0,1},
	};

	const int Xnum = 101;//x�����ĸ���20,40,150
	const int Ynum = 101;//y�����ĸ���43,87,215
	const int Pnum = Xnum * Ynum;//һ����ĸ���
	const double dx = 0.0002;
	const double dy = dx;
}
struct Flux
{
	double f1;
	double f2;
	double f3;
	double f4;
};

struct mesh
{
	double x;
	double y;
	int id;
	double rho;//�ܶ�
	double u;//ˮƽ�ٶ�
	double v;//��ֱ�ٶ�
	double p;//ѹǿ
	double x0;//������ĳ�ʼ����
	double y0;
	double um;//����˶��ٶ�
	double vm;
	int section = 1;//����
	int sec_num;//���������Ĳ���1��0
	int step = 0;//�����õ��ںβ�����µ�
	int neiborsec = -1;//���ڷ����������ֵΪ�������ʾ�õ㲻�Ǽ����㣬��֮Ϊ���ڷ������õļ�����
	int neiborsec_ad = -1;//���ڷ����ĵ�ַ
	vector <double> xix;
	vector <double> xiy;
	vector <double> xit;
	vector <double> etax;
	vector <double> etay;
	vector <double> etat;
	vector <double> xxi;
	vector <double> xeta;
	vector <double> xtau;
	vector <double> yxi;
	vector <double> yeta;
	vector <double> ytau;
	vector <double> J;//�ſɱ�����ʽ
	vector <mesh*>neibor;//��¼�ø������ڸ��λ����Ϣ
	vector <int>moveConnct;//�˶������㣬�������˶����moveConnect���й�
	string type;//������ͣ���Ϊ�������ұ߽��Լ��ڲ���U��D,L,R,IN,������SHOCK���Ӵ���ϵ�DISCON��contact discontinuity��,�����ཻ��CENTER
	string change = "Y";//��Χ�����������Ƿ�ı䣬û�С�N�����ı䡰Y�������ڼ������
	Flux alpha = { 0 };
	Flux beta = { 0 };


};
struct Coor//��ʾ����
{
	double x;
	double y;
};
struct Line//��ʾֱ��Ax+Bx+C=0
{
	double A;
	double B;
	double C;
};
struct polygon_mesh
{
	//��ֻ��¼�����ţ�����������������
	vector <int> node;
	vector <int>face_start;
	vector <int>face_end;
	int section;
};
namespace Oblique//б����
{
	extern double beta;
	extern double delta;
	extern int startpoint;
	extern double  rho1, p1, Ma1, u1, v1;
	extern double  rho2, p2, Ma2, u2, v2;
}
namespace Prandtl_Meyer//��������Ү������
{
	extern double  rho0, p0, theta0, lambda0, delta0, mu0, Ma0, u0, v0;
	extern double  rho1, p1, theta1, lambda1, delta1, mu1, Ma1, u1, v1;
	extern double  rho2, p2, theta2, lambda2, delta2, mu2, Ma2, u2, v2;
}
namespace ShockwaveCross//�����ཻ
{

	extern double  rho1, p1, beta1, delta1, Ma1, u1, v1;
	extern double  rho2, p2, beta2, delta2, Ma2, u2, v2;
	extern double  rho3, p3, beta3, delta3, Ma3, u3, v3;
	extern double  rho4, p4, beta4, delta4, Ma4, u4, v4;
	extern double  rho5, p5, beta5, delta5, Ma5, u5, v5;
}

