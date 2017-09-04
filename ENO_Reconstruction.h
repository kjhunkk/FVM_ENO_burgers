#pragma once

#include <vector>

class ENO_Reconstruction
{
public:
	// constructor / p.m. base point, grid size
	ENO_Reconstruction(int);
	~ENO_Reconstruction();

public:
	// functions
	// set grid size / p.m. grid size
	static void setDX(double dx) { m_step = dx; };

	// set order of polynomial / p.m. order of polynomial
	static void setOFP(int n) { m_OFP = n; };

	// set cell X coordinate / p.m. x coordinate
	static void setX(std::vector<double> x) { m_cell_X = x; };

	// set cell averaged quantity / p.m. cell-averaged value
	static void setQ(std::vector<double> q) { m_cell_Q = q; };

	// set stencil begin
	void set_stencil() { m_stencil = stencil_begin(m_primitive, m_base, m_OFP); };

	// update primitive value
	static void update();

	// reconstruction. should have been updated before use / p.m. x coordinate / r.t. reconsturction value
	double reconstruct(double);

	// reconstruction by degree of 3. should have been updated befor use / p.m. x coordinate / r.t. reconstruction value
	double reconstruct3(double);

protected:
	// variables
	int m_stencil;
	int m_base;

	// static variables
	static double m_step;
	static int m_OFP;
	static std::vector<double> m_primitive;
	static std::vector<double> m_cell_X;
	static std::vector<double> m_cell_Q;

	// functions
	// divided difference / p.m. quantity, i, i + k / r.t. divided difference
	double divDiff(std::vector<double>, int, int);

	// calculate adaptive stencil / p.m. quantity, base, order of polynomial / r.t. stencil begin
	int stencil_begin(std::vector<double>, int, int);

	// calculate reconstruction for cell j / p.m. quantity, stencil begin, order of polynomial / r.t. reconstruction for cell j
	double cell_reconst(std::vector<double>, int, int);
};