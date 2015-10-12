#include "math.h"

// test data
//static const int N = 20;
//static const int N1 = 10;
//static const int q = 3;
//static const int qq = 5;
//static const double hx = 1.0 / N1;
//static const double hy = 1.0 / N;
//static const int m = 1;
//--------------------------------------

// real data
static const int N = 1200;
static const int N1 = 800;
static const int q = 101;
static const int qq = 20;
static const double hx = 1.0 / 100;
static const double hy = 1.0 / 200;
static const int m = 25000;
//--------------------------------------
static const int M1 = N1 + 1;
static const int M = N + 1;
static const int M2 = M1 * M;
static const double tau = 0.0005;
static const double tg = 2;
static const int cntr = N / 2;
static const int w = q;
//M1 - ���������� ����� �� ��� �
//M - ���������� ����� �� ��� y
//(2*q-1) - ���������� ����� � ��������� �����
//(qq,cntr) - ����� ���� ������� �����
//tg = hx/hy
//m - ���������� ����� �� �������

static const double Gamma = 1.4;
static const double Re = 10000;
static const double Mah = 4;
static const double Pr = 0.72;
static const double omega = 0.8;
static const double epsilon = 0.0000000001;
static double A[2 * M2][12];
static double D[2 * M2];
static double f[2 * M2];
static double Sigma_k[M2];
static double Sigma_k1[M2];
static double u_k[M2];
static double u_k1[M2];
static double v_k[M2];
static double v_k1[M2];
static double B[2 * M2];
static double u2[M2];
static double v2[M2];
static double e_k[M2];
static double e_k1[M2];
static double e2[M2];
static double T[M2];
static double Sigma_kk[M2];
static double u_kk[M2];
static double v_kk[M2];
static double e_kk[M2];
static double SigmaX_k[M2];
static double uX_k[M2];
static double vY_k[M2];
static double eR_k[M2];

// � �������� � _k1 �������� �������� ������� �� d-�� ���� �� �������
// � �������� � _k �������� �������� ������� � ���������� �������� �� ������������
// � �������� � _kk �������� �������� ������� c (d-1) ���� �� �������
// ������� � "2" ����������� � ��������� ������ �������
// � �������� � X_k �������� �������� �������, ����������� ������� ����������

int s_m, s_e, s_itr, s_end, itr = 5, itr_tr = itr;


inline double Mu(double e_k)
{
	return pow(Gamma * (Gamma - 1) * Mah * Mah * e_k * e_k, omega);
}

inline double P(double sigma_k, double e_k)
{
	return (Gamma - 1) * sigma_k * sigma_k * e_k * e_k;
}