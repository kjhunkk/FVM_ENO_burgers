#pragma once

#include <vector>

#include "Flux.h"
#include "ENO_Reconstruction.h"

class ENO_Flux : public Flux
{
public:
	// constructor / p.m. direction of flux, dx, dt, advection speed, x coordinate
	ENO_Flux(bool, const double, const double, const double, std::vector<double>);
	~ENO_Flux();

public:
	// functions
	// set reconstruction object / p.m. reconstruction
	static void setReconstruction(std::vector<ENO_Reconstruction> Rw) { m_Rw = Rw; };

	// calculate flux / p.m. index / e.g. #0 : backward, #1 : forward
	double flux(const int);

protected:
	// variables
	static std::vector<ENO_Reconstruction> m_Rw;
};