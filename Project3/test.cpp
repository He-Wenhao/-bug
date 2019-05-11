#include<algorithm>
#include<iostream>
#include<fstream>
#include<iostream>
#include<complex>
#include"const.h"
constexpr double E = 2.71828182845904523536028747;
constexpr double PI = 3.14159265358979323846264338;
using namespace std;
using com = complex<double>;
//x轴离散成2*N0个格子
//constexpr int N0 = 16384;
//x间隔为delta_x
//constexpr double delta_x = 0.1;
//归一化
template<int n>
void normalize(com* x,double delta_x) {
	double norm=0;
	for (int i = 0; i < n; i++) {
		norm += pow(abs(x[i]),2)*delta_x;
	}
	norm = sqrt(norm);
	for (int i = 0; i < n; i++) {
		x[i] = x[i] / norm;
	}
}
//特征矢+特征值
struct eigen {
	com* vec;
	double val;
};
//打印向量
template<int n>
void print_sol(com x[n]) {
	for (int i = 0; i < n; i++) { cout << i + 1 << "   " << x[i] << endl; }
}

//解对称三对角线性方程组
template<int n>
com* solve_eq(com c_[n][2], com b_[n]) {
	
	//complex c[n][2];
	com(*c)[2] = new com[n][2];
	com* b = new com[n];
	for (int i = 0; i < n; i++) {
		c[i][0] = c_[i][0];
		c[i][1] = c_[i][1];
		b[i] = b_[i];
	}
	//计算L矩阵元
	for (int i = 1; i <= n; i++) {
		int r = max(1, i - 1);
		for (int j = r; j <= i; j++) {
			for (int k = r; k <= j - 1; k++) {
				//由于数组下标从0开始,故需要-1
				c[i - 1][j - i + 2 - 1] = c[i - 1][j - i + 2 - 1] - c[i - 1][k - i + 2 - 1] * c[j - 1][k - j + 2 - 1] / c[k - 1][2 - 1];
			}
		}
	}
	//计算b'
	for (int i = 1; i <= n; i++) {
		int r = max(1, i - 1);
		for (int j = r; j <= i - 1; j++) {
			b[i - 1] = b[i - 1] - c[i - 1][j - i + 2 - 1] * b[j - 1] / c[j - 1][2 - 1];
		}
	}
	//上三角回代
	for (int i = n; i >= 1; i--) {
		int t = min(n, i + 1);
		for (int j = i + 1; j <= t; j++) {
			b[i - 1] = b[i - 1] - c[j - 1][i - j + 2 - 1] * b[j - 1];
		}
		b[i - 1] = b[i - 1] / c[i - 1][2 - 1];
	}
	delete[] c;
	return b;
}

//反幂法求本征问题
template<int N>
com* inver_exp(com C[2*N-1][2], com u_[2*N-1]) {
	//初始化u,lambda
	constexpr int n = 2 * N - 1;
	com* u=new com[n];
	for (int i = 0; i < n; i++) {
		u[i] = u_[i];
	}
	com lambda = 0;
	com* v=new com[n];
	for (int k = 0; k < 10000; k++) {
		delete[] v;
		v = nullptr;
		v = solve_eq<n>(C, u);
		//求v绝对值最大的分量max
		com max = 0;
		for (int i = 0; i < n; i++) {
			if (abs(max) < abs(v[i])) {
				max = v[i];
			}
		}
		//储存前一次本征值
		com lambda_before = lambda;
		lambda = 1. / max;
		//u=v*lambda
		for (int i = 0; i < n; i++) {
			u[i] = v[i]*lambda;
		}
		//检查精度是否合适
		if (abs(lambda_before - lambda)/abs(lambda) < 1e-10) {
			break;
		}
	}
	delete[] v;
	//return eigen{ u,abs(lambda) };
	return u;
}


//乘高斯吸收函数
template<int n>
void absorb(com A[n],double delta_x) {
	int N0 = (n + 1) / 2;
	double x0 = 0.75*N0*delta_x;
	double x = -N0 * delta_x;
	for (int i = 0; i < n; i++) {
		x += delta_x;
		if (abs(x) > x0) {
			com temp = pow(E, -pow((abs(x) - x0) / 0.2, 2));

			com temp2 = (A)[i];


			(A)[i] = temp2 * temp;
		}

	}
}


//波函数传播delta_t
//其中I(10^16为单位),omega(原子单位),N是激光场参数
template<int N0>
void Udelta_t(com psai[2*N0-1], double t, double I, double omega, int N, double delta_t,double delta_x) {
	constexpr int n = 2 * N0 - 1;
	//构造哈密顿量H,装在C中
	com(*C)[2] = new com[n][2];
	C[0][0] = 0;
	for (int i = 1; i < n; i++) {//构造非对角元
		C[i][0] = 0;
	}
	for (int i = 0; i < n; i++) {//构造对角元
		C[i][1] = 1;
	}


	//计算psai*(1-0.5*i*H*delta_t)
	com *psai_temp = new com[n];//临时储存psai
	for (int i = 0; i < n; i++) {
		psai_temp[i] = psai[i];
	}
	
	delete[] psai;
	psai = nullptr;
	psai = solve_eq<n>(C, psai_temp);//开始计算
	/*
	com* temp_solve = solve_eq<n>(C, psai_temp);
	for (int i = 0; i < n; i++) {
		psai[i] = temp_solve[i];
	}
	delete[] temp_solve;
	*/
	delete[] C;
	delete[] psai_temp;



}


//生成t=0的波函数
template<int N0>
com* generate_t0(double delta_x) {
	//构造三对角型的系数矩阵A,装在C中
	com (*C)[2]=new com[2*N0-1][2];
	constexpr int n = 2 * N0 - 1;
	C[0][0] = 0;
	double x = -N0 * delta_x;
	for (int i = 1; i < n; i++) {//构造非对角元
		C[i][0] = -0.5 / pow(delta_x, 2);
	}
	for (int i = 0; i < n; i++) {//构造对角元
		x += delta_x;
		C[i][1] = 1 / pow(delta_x, 2) - 1 / sqrt(2 + pow(x, 2)) + 0.48;
	}
	//u 为初始向量
	com* u = new com[n];
	for (int i = 0; i < n; i++) {
		u[i] = 1 / sqrt(2 * N0*delta_x);
	}

	com* psai0 = inver_exp<N0>(C, u);
	normalize<n>(psai0,delta_x);
	delete[] u,C;
	return psai0;
}


//求加速度和t的函数
com* generate_at() {
	double I = 0.01;
	double omega300 = 45.5633525316 / 400;
	constexpr int N0 = 6000;//取2N0-1个坐标点
	double delta_x = 0.1;
	int Ndata = 400;//取400个时间点
	int N = 4;
	constexpr int n = 2 * N0 - 1;
	double tf = 2 * N*PI / omega300;
	double sqrtI = sqrt(I / 3.5094448314);
	//对应的delta_t
	double delta_t = tf / Ndata;
	com* at = new com[Ndata];
	com* psai = generate_t0<N0>(delta_x);
	double t = 0;
	for (int i = 0; i < Ndata; i++) {
		Udelta_t<N0>(psai, t, I, omega300, N, delta_t, delta_x);//波函数传播
		absorb<n>(psai,delta_x);//乘吸收函数
		//求t时刻at
		double x = -N0 * delta_x;
		for (int j = 0; j < 2 * N0 - 1; j++) {
			x += delta_x;
			//受力算符 -dV/dx+E
			com F = -x / pow(x*x + 2, 1.5) + sqrtI * pow(sin(omega300*t / 2. / N), 2)*sin(omega300*t);
			at[i] += conj(psai[j])*psai[j] *F* delta_x;
		}
		t += delta_t;
		if (i % 100 == 0) {//!!!!!!!!!!!!!
			cout << i << endl;
		}
	}
	delete[] psai;
	return at;
}

void testbug1() {

	generate_at();
}

void testbug2() {
	constexpr int N0 = 6000;//取2N0-1个坐标点
	double delta_x = 0.1;
	constexpr int n = 2 * N0 - 1;
	int N = 4;
	double I = 0.01;
	double omega300 = 45.5633525316 / 400;
	double tf = 2 * N*PI / omega300;
	int Ndata = 400;//取400个时间点
	double delta_t = tf / Ndata;



	com* psai = generate_t0<N0>(delta_x);
	Udelta_t<N0>(psai, 1, I, omega300, N, delta_t, delta_x);//波函数传播
	absorb<n>(psai, delta_x);//乘吸收函数
}

int main() {
	testbug2();
	system("pause");
}


