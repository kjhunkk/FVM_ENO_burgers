#pragma once
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

typedef int direction;
typedef int Type;

class ControlArea
{
public:
	// constructor / p.m. number of cell / d.f. number of ghost cell of one side = 2
	ControlArea(int);

	// constructor / p.m. number of cell, number of ghost cell of one side
	ControlArea(int, int);

	// constructor / p.m. number of cell, number of ghost cell of one side, boundary value begin, boundary value end
	ControlArea(int, int, double, double);

	// destructror
	~ControlArea();

public:
	// functions
	inline void setCellQ(std::vector<double> Q) { m_cell_Q = Q; };

	inline double getBV_begin() const { return m_boundary_value_begin; };

	inline double getBV_end() const { return m_boundary_value_end; };

	inline double getStep() const { return m_step; };

	inline int getSideGhost() const { return m_nSideGhost; };

	inline std::vector<double> getCell_Q() const { return m_cell_Q; };

	inline std::vector<double> getCell_Q_previous() const { return m_cell_Q_previous; };

	inline std::vector<double> getCell_X() const { return m_cell_X; };

	// apply boundary condition / p.m. boundary condition type
	void BC(Type);

	// apply boundary condition / p.m. boundary condition type, cell quantity
	void BC(Type, std::vector<double>&);

	// update cell variables / e.g. previous Q
	void update();

	// export solution variables / p.m. file index, initial condition, flux type, order of polynomial
	bool output(std::string, Type, Type, int) const;

	// initialize cell quantity with averaging & generate cell coordinate / p.m. area, initial condition
	void initialize(double, Type);

protected:
	// variables
	double m_boundary_value_begin;
	double m_boundary_value_end;
	double m_step;
	int m_nSideGhost;
	int m_NFX;

	// solution variables
	std::vector<double> m_cell_Q;
	std::vector<double> m_cell_Q_previous;
	std::vector<double> m_cell_X;

protected:
	// functions
	// assign initial conditions / p.m. grid size / gauss quadrature order 5
	void initial_condition_0(double);
	void initial_condition_1(double);
	void initial_condition_2(double);

	// initial condition functions / p.m. x coordinate
	double function_0(double);
	double function_1(double);
	double function_2(double);
};