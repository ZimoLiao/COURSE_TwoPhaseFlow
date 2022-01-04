#include<iostream>
#include<string>
#include<time.h>

#include"LevelSetSolver.h"

using namespace std;

// 改变模拟参数和速度/界面初始化类型可以校验不同算例、数值格式
int main()
{
	// 模拟参数
	//	n		网格量
	//	s		空间离散格式
	//			1	QUICK
	//	t		时间推进格式
	//			1	forward Euler
	//			2	TVD-RK2
	//			3	TVD-RK3
	int n = 256, s = 1, t = 3, r = 2;

	// 模拟时间（固定最大cfl数计算）
	//	tmax	终止时间
	//	step_r	重新初始化间隔步数
	double umax = 2.0 * pi;
	double cfl = 0.8;
	double dt = cfl * 1.0 / n / umax, tmax = 0.8;
	int step_r = floor(0.1 / dt);

	LevelSetSolver ls(n, dt, tmax, s, t, r, step_r);


	// 速度设置
	//	1	剪切流
	//	2	tearing flow (Rider & Kothe 1995)
	//	3	旋转流
	ls.InitVelocity(2);

	// 初始界面设置
	//	1	液滴（作业1）
	//	2	缺角圆盘（作业2）
	//	3	液滴（Rider & Kothe 1995）
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