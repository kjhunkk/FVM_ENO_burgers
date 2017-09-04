#include "ENO_Flux.h"
#include <iostream>

std::vector<ENO_Reconstruction> ENO_Flux::m_Rw;

ENO_Flux::ENO_Flux(bool direct, const double DX, const double DT, const double a, std::vector<double> x)
	:Flux(direct, DX, DT)
{
	m_cell_X = x;
	m_speed = a;
}

ENO_Flux::~ENO_Flux()
{

}

double ENO_Flux::flux(const int index)
{
	double flux = 0.0;
	if (m_Rw[index].reconstruct(m_cell_X[index]) >= 0.0)
	{
		if (m_direction == forward) flux = 0.5*pow(m_Rw[index].reconstruct(m_cell_X[index] + 0.5*m_dx), 2.0);
		else flux = 0.5*pow(m_Rw[index - 1].reconstruct(m_cell_X[index] - 0.5*m_dx), 2.0);
	}
	else
	{
		if (m_direction == forward) flux = 0.5*pow(m_Rw[index + 1].reconstruct(m_cell_X[index] + 0.5*m_dx), 2.0);
		else flux = 0.5*pow(m_Rw[index].reconstruct(m_cell_X[index] - 0.5*m_dx), 2.0);
	}

	return flux;
}