#pragma once

#include <memory>
#include <vector>

#include "ENO_Flux.h"
#include "ENO_Reconstruction.h"
#include "ControlArea.h"
#include "Condition.h"

#define forward 1
#define backward 0

typedef int Type;

class RHS
{
public:
	// constructor / p.m. Condition, ControlArea / limiter type = none
	RHS(std::shared_ptr<Condition>, std::shared_ptr<ControlArea>);

	// destructor
	~RHS();

public:
	// functions
	// calculate RHS with flux / p.m. cell quantity, dt
	std::vector<double> cal_RHS(std::shared_ptr<ControlArea>, double);

protected:
	// variables
	Type m_bc;
	double m_dx;
	double m_dt;
	int m_NFX;
	int m_nSideGhost;
	std::vector<double> m_Q;
	std::vector<double> m_Q1;
	std::vector<double> m_Q2;
	std::vector<double> m_L; // L(Qn)
	std::shared_ptr<ENO_Flux> m_flux_p;
	std::shared_ptr<ENO_Flux> m_flux_m;
	std::vector<ENO_Reconstruction> m_Rx;

protected:
	// functions
	// RHS solver / p.m. Q / r.t. L
	std::vector<double> solver(const std::vector<double>);
};