#pragma once
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

typedef int Type;

class Condition
{
public:
	Condition();
	~Condition();

public:
	inline double getAdvSpeed() const { return m_speed; };

	inline double getGridSize() const { return m_dx; };

	inline double getTimeStep() const { return m_dt; };

	inline double getCFL() const { return m_CFL; };

	inline double getTargetTime() const { return m_target_time; };

	inline double getArea() const{ return m_area; };
	
	inline int getOFP() const { return m_OFP; };

	inline int getNFX() const { return m_NFX; };

	inline Type getFlux() const { return m_flux; };

	inline Type getInitial() const { return m_initial; };

	inline Type getBC() const { return m_BC; };

	// read input file
	bool input();

	// update conditions / e.g. dt
	void update(std::vector<double>);

protected:
	// variables
	double m_speed;
	double m_dx;
	double m_dt;
	double m_CFL;
	double m_target_time;
	double m_area;
	int m_NFX;
	int m_OFP;
	Type m_flux;
	Type m_initial;
	Type m_BC;

protected:
	// functions
	// read input file
	bool readInput();

	// read option file
	bool readOption(std::string, std::vector<std::string>&);

	// print options & input a condition
	bool inputType(Type&, std::string, std::vector<std::string>);

	// print conditions
	void printConditions() const;
};