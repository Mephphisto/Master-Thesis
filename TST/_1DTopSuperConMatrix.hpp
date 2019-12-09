#pragma once
#include <Eigen/Dense>
class _1DTopSuperConMatrix
{
public:
	_1DTopSuperConMatrix(size_t size_in, double t_in, double mu_in, double Delta_in)
	{
		this->t = -t_in;
		this->mu = -mu_in;
		this->Delta = Delta_in;
		Gsize = size_in / 2;
		m = Eigen::MatrixXcd(size_in,size_in);
		m.setZero();
		Build_A();
		//Build_B();
	}
	Eigen::MatrixXcd get()
	{
		return m;
	}
	double trace_A()
	{
		return 2 * Gsize * mu;
	}
private:
	Eigen::MatrixXcd m;
	double t, mu, Delta;
	size_t Gsize;
	void Build_A()
	{
		for (size_t i = 0; i < Gsize - 1; i++)
		{
			m(i, i) = mu;
			m(i + Gsize, i + Gsize) = mu;
			m(i, i+1) = t;
			m(i+1, i) = t;
			m(i + Gsize, i + Gsize + 1) = t;
			m(i + Gsize + 1, i + Gsize ) = t;
			m(i, i + 1 + Gsize) = Delta;
			m(i + 1, i + Gsize) = Delta;
			m(i + Gsize, i + 1) = Delta;
			m(i + 1 + Gsize, i) = Delta;
		}
		m(Gsize-1, Gsize-1) = mu;
		m(2*Gsize - 1, 2*Gsize - 1) = mu;
	}
};

