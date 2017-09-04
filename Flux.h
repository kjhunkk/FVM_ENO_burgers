#pragma once

#include <vector>

#define s_f shock_speed_forward
#define s_b shock_speed_backward

#define forward 1
#define backward 0

class Flux
{
public:
	// constructor / p.m. direction of flux, dx, dt
	Flux(bool, const double, const double);

	// destructor
	virtual ~Flux();

public:
	// functions
	// set cell quantities for flux calculation / p.m. cell quantity
	static inline void setCellQ(std::vector<double> Q) { m_cell_Q = Q; };

	inline void setDt(double dt) { m_dt = dt; };

	inline bool getDirection() const { return m_direction; };

	// calculate flux / p.m. index / e.g. #0 : backward, #1 : forward
	virtual double flux(const int);

protected:
	// variables
	bool m_direction; // #0 : backward, #1: forward
	double m_speed;
	double m_dx;
	double m_dt;
	std::vector<double> m_cell_X;
	static std::vector<double> m_cell_Q;

protected:
	// functions
	// calculate shock speed / p.m. directional
	std::vector<double> cal_shock_speed() const;
};