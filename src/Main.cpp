#include<iostream>
#include<string>
#include<time.h>

#include"LevelSetSolver.h"

using namespace std;

// �ı�ģ��������ٶ�/�����ʼ�����Ϳ���У�鲻ͬ��������ֵ��ʽ
int main()
{
	// ģ�����
	//	n		������
	//	s		�ռ���ɢ��ʽ
	//			1	QUICK
	//	t		ʱ���ƽ���ʽ
	//			1	forward Euler
	//			2	TVD-RK2
	//			3	TVD-RK3
	int n = 256, s = 1, t = 3, r = 2;

	// ģ��ʱ�䣨�̶����cfl�����㣩
	//	tmax	��ֹʱ��
	//	step_r	���³�ʼ���������
	double umax = 2.0 * pi;
	double cfl = 0.8;
	double dt = cfl * 1.0 / n / umax, tmax = 0.8;
	int step_r = floor(0.1 / dt);

	LevelSetSolver ls(n, dt, tmax, s, t, r, step_r);


	// �ٶ�����
	//	1	������
	//	2	tearing flow (Rider & Kothe 1995)
	//	3	��ת��
	ls.InitVelocity(2);

	// ��ʼ��������
	//	1	Һ�Σ���ҵ1��
	//	2	ȱ��Բ�̣���ҵ2��
	//	3	Һ�Σ�Rider & Kothe 1995��
	ls.InitPhi(3);

	clock_t start, end;
	start = clock();
	ls.Calculation();
	end = clock();
	cout << "time = " << double(end - start) / CLOCKS_PER_SEC << "s" << endl;

	string fname = "n" + to_string(n) + "s" + to_string(s) + to_string(t) + to_string(r) \
		+ "t" + to_string(int(tmax * 10.0)) + ".dat";
	ls.WriteAll(fname);
}