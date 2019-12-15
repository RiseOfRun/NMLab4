#include <functional>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;
typedef function<double(vector<double>)> Func;

class functions
{
public:
	static double F1(vector<double> x)
	{
		return 0;
	}

	static double F2(vector<double> x)
	{
		return 0;
	}

	static double dF11(vector<double> x)
	{

	}
	static double dF12(vector<double> x)
	{

	}

	static double dF21(vector<double> x)
	{

	}
	static double dF22(vector<double> x)
	{

	}

	vector <function<double(vector<double>)>> F{F1,F2};
	vector <vector<Func>> dF = vector<vector<Func>>{
	vector<Func>{dF11, dF12},
	vector<Func> {dF21, dF22}
	};

};

class SNE
{
public:
	int m, n;
	functions F;
	double e1, eps, norm0, normF ,Bk;
	int maxIter,iterations;
	vector<double> x0, xk,dxk,f0, f;
	vector<vector<double>> A;
	SNE(ifstream &in, ifstream &pr)
	{
		in >> n >> m >> maxIter >> e1 >> eps;
		x0 = vector<double>(n);
		for (size_t i = 0; i < n; i++)
		{
			pr >> x0[i];
		}

		iterations = 0;
		xk = x0;
		f0 = CalculateF(x0);
		norm0 = Norm(f0);
	}

	void Sort()
	{
		int count = m - n + 1;
		vector<Func> sortedFuncs;
		vector<vector<Func>> sortedDfs;
		vector<double> values;
		for (size_t i = 0; i < m; i++)
		{
			double tmp = F.F[i](xk);
			bool pushed = false;
			for (size_t j = 0; j < values.size; j++)
			{
				if (abs(tmp) > abs(values[j]))
				{
					values.emplace(values.begin() + j, tmp);
					sortedFuncs.emplace(sortedFuncs.begin() + j, F.F[i]);
					sortedDfs.emplace(sortedDfs.begin() + j, F.dF[i]);
					pushed = true;
					break;
				}
			}
			if (!pushed)
			{
				values.push_back(tmp);
				sortedFuncs.push_back(F.F[i]);
				sortedDfs.push_back(F.dF[i]);
			}
		}
		f = values;
		F.dF = sortedDfs;
		F.F = sortedFuncs;
	}

	void JacobiV2()
	{
		Sort();
		A = vector<vector<double>>(n);
		for (size_t i = 0; i < n; i++)
		{
			A[i] = vector<double>(n);
			for (size_t j = 0; j < n; j++)
			{
				A[i][j] = F.dF[i][j](xk);
			}
		}
	}

	void JacobiV4()
	{
		A = vector<vector<double>>(m);
		for (size_t i = 0; i < m; i++)
		{
			A[i] = vector<double>(n);
			for (size_t j = 0; j < n; j++)
			{
				A[i][j] = F.dF[i][j](xk);
			}
		}
		f = CalculateF(xk);
		f = MultMTF(xk);
		A = Symetric();
	}

	vector<vector<double>> Symetric()
	{
		vector<vector<double>> result(n);
		for (size_t i = 0; i < n; i++)
		{
			result[i] = vector<double>(m);
			for (size_t j = 0; j < n; j++)
			{
				double sum = 0;
				for (size_t k = 0; k < m; k++)
				{
					result[i][j] += A[k][i] * A[k][j];
				}
				
			}
		}
		return result;
	}

	vector<double> MultMTF(vector<double> &x)
	{
		vector<double> result(m);
		for (size_t j = 0; j < n; j++)
		{
			double sum = 0;
			for (size_t i = 0; i < m; i++)
			{
				result[i] -= A[i][j] * x[j];
			}
		}
		return result;
	}

	vector<double> Gauss(vector<double> b, int n)
	{
		for (size_t i = 0; i < n; i++)
		{
			double max = abs(A[i][i]);
			int l = i;
			for (int j = i; j < n; j++)
			{
				if (abs(A[j][i]) > max)
				{
					l = j;
					max = abs(A[j][i]);
				}
			}

			std::swap(A[i], A[l]);
			std::swap(b[i], b[l]);

			for (size_t j = i + 1; j < n; j++)
			{
				double mult = A[j][i] / A[i][i];
				if (mult ==0)
				{
					throw exception("BadMatrix");
				}
				for (size_t k = i; k < n; k++)
				{
					if (k == i)
					{
						A[j][k] = 0;
					}
					else A[j][k] -= mult * A[i][k];
				}
				b[j] -= b[i] * mult;
			}
		}

		vector<double>x(n);
		for (int i = n - 1; i >= 0; i--)
		{
			x[i] = b[i] / A[i][i];
			for (int j = 0; j < i; j++)
			{
				b[j] -= x[i] * A[j][i];
			}
		}
		return x;
	}

	vector<double> CalculateF(vector<double>& x)
	{
		vector<double> result(m);
		for (size_t i = 0; i < m; i++)
		{
			result[i] = F.F[i](x);
		}
		return result;
	}

	double Norm(vector<double>& x)
	{
		double result = 0;
		for (size_t i = 0; i < x.size; i++)
		{
			result += x[i] * x[i];
		}
		return result;
	}

	vector<double> Shift(double B)
	{
		vector<double> result(n);
		for (size_t i = 0; i < n; i++)
		{
			result[i] = xk[i] + B * dxk[i];
		}
		return result;
	}

	double FindB()
	{
		double B = 1;
		vector<double> shifted = Shift(B);
		vector<double> Fv = CalculateF(shifted);
		double NormV = Norm(Fv);
		while (NormV>norm0)
		{
			B = B / 2;
			if (B<e1)
			{
				return B;
			}
			shifted = Shift(B);
			Fv = CalculateF(shifted);
			NormV = Norm(Fv);
		}
		normF = NormV;
		return B;
	}

	bool Step()
	{
		if (normF/norm0<=eps*eps||iterations>=maxIter||Bk<e1)
		{
			return false;
		}

		dxk = Gauss(f, n);
		Bk = FindB();
		if (Bk < e1)
		{
			return false;
		}
		dxk = Shift(Bk);
		return true;
	}

	~SNE();
private:

};

int main()
{
	ifstream in,pr;
	in.open("input.txt");
	pr.open("pr.txt");
	SNE Sys(in, pr);
	bool flag = true;
	while (flag)
	{
		Sys.JacobiV2();
		try
		{
			flag = Sys.Step();
		}
		catch (const std::exception& ex)
		{
			string messege(ex.what);
			if (messege == "BadMatrix")
			{
				cout << messege << endl;
				return -1;
			}
		}
	}
}
