#include "ENO_Reconstruction.h"

double ENO_Reconstruction::m_step = 0;
int ENO_Reconstruction::m_OFP = 0;
std::vector<double> ENO_Reconstruction::m_primitive;
std::vector<double> ENO_Reconstruction::m_cell_X;
std::vector<double> ENO_Reconstruction::m_cell_Q;

ENO_Reconstruction::ENO_Reconstruction(int i)
{
	m_base = i;
	m_primitive.clear();
}

ENO_Reconstruction::~ENO_Reconstruction()
{

}

void ENO_Reconstruction::update()
{
	m_primitive.clear();
	m_primitive.push_back(m_step*m_cell_Q[0]);
	for (int i = 1; i < m_cell_Q.size(); ++i)
	{
		m_primitive.push_back(m_step*m_cell_Q[i] + m_primitive[i - 1]);
	}
}

double ENO_Reconstruction::reconstruct(double x)
{
	double L = 0.0;
	double No;
	double De;
	for (int j = m_stencil; j <= m_stencil + m_OFP; ++j)
	{
		No = 0.0;
		De = 1.0;
		for (int i = m_stencil; i <= m_stencil + m_OFP; ++i)
		{
			if (j != i)
			{
				De *= (m_cell_X[j] - m_cell_X[i]);
			}
		}
		double temp;
		for (int i = m_stencil; i <= m_stencil + m_OFP; ++i)
		{
			if (i != j)
			{
				temp = 1.0;
				for (int k = m_stencil; k <= m_stencil + m_OFP; ++k)
				{
					if ((j != k) && (i != k))
					{
						temp *= (x - (m_cell_X[k] + 0.5*m_step));
					}
				}
				No += temp;
			}			
		}

		L += m_primitive[j] * No / De;
	}

	return L;
}

double ENO_Reconstruction::reconstruct3(double X)
{
	double L = 0.0;
	double x[4] = { m_cell_X[m_stencil] + 0.5*m_step, m_cell_X[m_stencil + 1] + 0.5*m_step, m_cell_X[m_stencil + 2] + 0.5*m_step, m_cell_X[m_stencil + 3] + 0.5*m_step };
	double y[4] = { m_primitive[m_stencil], m_primitive[m_stencil + 1], m_primitive[m_stencil + 2], m_primitive[m_stencil + 3] };
	L += 0.5*(3.0*y[1] - y[0] - 3.0*y[2] + y[3])*X*X;
	L += ((x[1] + x[2] + x[3])*y[0] / 3.0 - (x[0] + x[2] + x[3])*y[1] + (x[0] + x[1] + x[3])*y[2] - (x[0] + x[1] + x[2])*y[3] / 3.0)*X;
	L += 0.5*(-((x[1] * x[2] + x[1] * x[3] + x[2] * x[3])*y[0] / 3.0) + ((x[0] * x[2] + x[0] * x[3] + x[2] * x[3])*y[1]) - ((x[0] * x[1] + x[0] * x[3] + x[1] * x[3])*y[2]) + ((x[0] * x[1] + x[0] * x[2] + x[1] * x[2])*y[3] / 3.0));
	L /= m_step*m_step*m_step;

	return L;
}

double ENO_Reconstruction::divDiff(std::vector<double> q, int i, int j)
{
	if (i == j) return q[i];
	return ((divDiff(q, i + 1, j) - divDiff(q, i, j - 1)) / (m_cell_X[j] - m_cell_X[i]));
}

int ENO_Reconstruction::stencil_begin(std::vector<double> q, int j, int m)
{
	int i = j - 1;
	for (int n = 1; n < m; ++n)
	{
		if (i <= 0) return 0;
		else if (i == 1) return 1;
		else if ((i + n + 1) >= (q.size() - 1))  return (q.size() - m - 1);
		else if (abs(divDiff(q, i - 1, i + n)) < abs(divDiff(q, i, i + n + 1))) i = i - 1;
	}

	return i;
}

double ENO_Reconstruction::cell_reconst(std::vector<double> q, int im, int m)
{
	std::vector<double> W(m + 1);
	for (int j = 0; j <= m; ++j)
	{
		W[j] = 0.0;
		for (int i = im; i < (im + j); ++i)
		{
			W[j] += 0.5*(m_cell_X[i + 1] - m_cell_X[i - 1])*q[i];
		}
	}

	return 0.0;
}
