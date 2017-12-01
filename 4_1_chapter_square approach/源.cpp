#include <iostream>
using namespace std;

void scanf_data(double ** &a, double *&b, double &x, int &n);
void calc_b(double temp_a[],double *b, const int n);
double integral_f(double temp_a[], int in, const int n);

void solution(double **a, double *b, int n);
void LU_decom(double **&a, double **&l, double **&u, int n);
void Calc_X(double **&l, double **&u, double *&b, int n);


int main()
{
	double **a = NULL;
	double *b = NULL, x = 0, temp_a[2] = {0,1};
	int n = 0;
	scanf_data(a, b,x, n);

	calc_b(temp_a,b,n);

	solution(a, b, n);

	getchar();
	getchar();
}

void calc_b(double temp_a[],double *b, const int n)
{

	for (int  i = 0; i < n; i++)
	{
		b[i] = integral_f(temp_a, i, n);
		cout << "\tb" << i << "= " << b[i];
	}
	cout << endl;
}
double integral_f(double temp_a[], int in, const int n)
{
	double piecewise = 5000;
	double len = temp_a[1] - temp_a[0];
	double sum = 0;
	double piece_len = len / piecewise;
	for (int i = 0; i < piecewise; i++)
	{
		double x = i * piece_len;
		sum += piece_len*sqrt(1 + x*x)*pow(x, in);
	}
	return sum;
}
void scanf_data(double ** &a,double *&b, double &x, int &n)
{
	cout << "********请输入数值的个数:**********\n";
	cin >> n;

	a = new double *[n];
	b = new double[n];

//	cout << "\n\n********请先输入" << n << "个x,然后输入" << n << "个f(x):**********\n";
	for (size_t i = 0; i < n; i++)
	{
		a[i] = new double[n];
		for (size_t j = 0; j < n; j++)
		{
			a[i][j] = 1.0 / (i+j+1);
			cout << a[i][j] << " ";
		}
		cout << endl;
	}

	cout << "\n\n********需要计算的x的值**********\n";
	cin >> x;

}

void solution(double **a, double *b, int n)
{
	double **L = NULL, **U = NULL;
	//LU分解
	LU_decom(a, L, U, n);
	//根据LU计算x的值
	Calc_X(L, U, b, n);
}

void LU_decom(double **&a, double **&l, double **&u, int n)
{
	//初始化空间l u
	l = new double *[n];
	u = new double *[n];
	for (int i = 0; i < n; i++)
	{
		l[i] = new double[n];
		u[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			l[i][j] = 0;
			u[i][j] = 0;
			if (i == j)
				l[i][j] = 1;
		}
	}

	//计算l u矩阵的结果
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (j == 0 && i != 0)
				l[i][j] = a[i][j] / u[0][0];
			else if (i == 0)
				u[i][j] = a[i][j];
			else {
				int min_temp = i > j ? j : i;

				if (i > j)
				{
					double sum = 0;
					for (int k = 0; k < min_temp; k++)
					{
						sum += l[i][k] * u[k][j];
					}
					l[i][j] = (a[i][j] - sum) / u[min_temp][min_temp];
				}
				else {
					double sum = 0;
					for (int k = 0; k < min_temp; k++)
					{
						sum += l[i][k] * u[k][j];
					}
					u[i][j] = (a[i][j] - sum);
				}
			}
			//cout << l[i][j] << " ";
		}
		//cout << endl;
	}

}

void Calc_X(double **&l, double **&u, double *&b, int n)
{
	double *y = NULL, *x = NULL;
	x = new double[n];
	y = new double[n];
	//计算y
	for (int i = 0; i < n; i++)
	{

		double sum_x = 0;
		for (int j = 0; j < i; j++)
		{
			sum_x += l[i][j] * y[j];
		}
		y[i] = b[i] - sum_x;
		//cout << y[i] << " ";
	}

	//计算x

	for (int i = n - 1; i >= 0; i--)
	{

		double sum_y = 0;
		for (int j = n - 1; j > i; j--)
		{
			sum_y += u[i][j] * x[j];
		}
		x[i] = (y[i] - sum_y) / u[i][i];
	}

	cout << "\n\n********a的值为**********：\n";
	for (int i = 0; i <n; i++)
	{
		cout << "\ta" << "[" << i << "]=" << x[i];
	}
}
