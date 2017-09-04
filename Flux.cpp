#include "Flux.h"

std::vector<double> Flux::m_cell_Q;

Flux::Flux(bool direct, const double DX, const double DT)
{
	m_direction = direct;
	m_dx = DX;
	m_dt = DT;
	m_cell_Q.resize(0);
}

Flux::~Flux()
{

}

double Flux::flux(const int index)
{
	double flux = 0.0;
	std::vector<double> shock_speed = cal_shock_speed();
	double adv_p = shock_speed[index - m_direction] > 0.0 ? shock_speed[index - m_direction] : 0.0;
	double adv_m = shock_speed[index + 1 - m_direction] < 0.0 ? shock_speed[index + 1 - m_direction] : 0.0;

	flux = adv_m*m_cell_Q[index + m_direction] + adv_p*m_cell_Q[index - 1 + m_direction];

	return flux;
}

std::vector<double> Flux::cal_shock_speed() const
{
	std::vector<double> shock_speed(m_cell_Q.size());

	for (int i = 1; i < m_cell_Q.size() - 1; ++i)
	{
		shock_speed[i] = m_cell_Q[i + m_direction] - m_cell_Q[i - 1 + m_direction];
	}

	return shock_speed;
}
