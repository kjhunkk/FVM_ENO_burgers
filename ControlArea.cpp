#include "ControlArea.h"

ControlArea::ControlArea(int N)
{
	// number of cells
	m_NFX = N;

	// default number of ghost cell(one side)
	m_nSideGhost = 2;

	// generate cell quantity array / size = NFX + ghost cells(2*nSideGhost)
	m_cell_Q.resize(m_NFX + 2 * m_nSideGhost);
	m_cell_Q_previous.resize(m_NFX + 2 * m_nSideGhost);
	m_cell_X.resize(m_NFX + 2 * m_nSideGhost);
}

ControlArea::ControlArea(int N, int ghost)
{
	// number of cells
	m_NFX = N;

	// default number of ghost cell(one side)
	m_nSideGhost = ghost;

	// generate cell quantity array / size = NFX + ghost cells(2*nSideGhost)
	m_cell_Q.resize(m_NFX + 2*m_nSideGhost);
	m_cell_Q_previous.resize(m_NFX + 2*m_nSideGhost);
	m_cell_X.resize(m_NFX + 2 * m_nSideGhost);
}

ControlArea::ControlArea(int N, int ghost, double begin, double end)
{
	// number of cells
	m_NFX = N;

	// default number of ghost cell(one side)
	m_nSideGhost = ghost;

	// generate cell quantity array / size = NFX + ghost cells(4)
	m_cell_Q.resize(m_NFX + 2*m_nSideGhost);
	m_cell_Q_previous.resize(m_NFX + 2*m_nSideGhost);
	m_cell_X.resize(m_NFX + 2*m_nSideGhost);

	// save boundary values
	m_boundary_value_begin = begin;
	m_boundary_value_end = end;
}

ControlArea::~ControlArea()
{

}

void ControlArea::BC(Type boundary_type)
{
	switch (boundary_type)
	{
		// constant boundary condition
	case 0:
		// ghost cells & boundary cells
		for (int i = 0; i < m_nSideGhost; ++i)
		{
			m_cell_Q[i] = m_boundary_value_begin;
			m_cell_Q[m_cell_Q.size() - 1 - i] = m_boundary_value_end;
		}
		break;

		// periodic boundary condition
	case 1:
		// ghost cells & boundary cells
		for (int i = 0; i < m_nSideGhost; ++i)
		{
			m_cell_Q[i] = m_cell_Q[m_cell_Q.size() - 2 * m_nSideGhost + i];
			m_cell_Q[m_cell_Q.size() - m_nSideGhost + i] = m_cell_Q[m_nSideGhost + i];
		}
		break;
	default:
		std::cout << "BC error\n";
	}

	m_cell_Q_previous = m_cell_Q;
}

void ControlArea::BC(Type boundary_type, std::vector<double>& q)
{
	switch (boundary_type)
	{
		// constant boundary condition
	case 0:
		// ghost cells & boundary cells
		for (int i = 0; i <= m_nSideGhost; ++i)
		{
			q[i] = m_boundary_value_begin;
			q[q.size() - 1 - i] = m_boundary_value_end;
		}
		break;

		// periodic boundary condition
	case 1:
		// ghost cells & boundary cells
		for (int i = 0; i < m_nSideGhost; ++i)
		{
			q[i] = q[q.size() - 2 * m_nSideGhost + i];
			q[q.size() - m_nSideGhost + i] = q[m_nSideGhost + i];
		}
		break;
	default:
		std::cout << "BC error\n";
	}
}

void ControlArea::update()
{
	m_cell_Q_previous = m_cell_Q;
}

bool ControlArea::output(std::string index, Type initial, Type flux, int OFP) const
{
	std::ofstream file;
	std::string FileName = "./output/";
	FileName += "ENO_Burgers_";
	std::string Scheme = "";
	switch (initial)
	{
	case 0:
		FileName += "shock_";
		break;
	case 1:
		FileName += "expansion_";
		break;
	case 2:
		FileName += "sine_";
		break;
	default:
		FileName += "TypeError_";
	}
	switch (flux)
	{
	case 0:
		FileName += "Upwind_";
		Scheme += "Upwind";
		break;
	case 1:
		FileName += "Lax Friedrich_";
		Scheme += "Lax Friedrich";
		break;
	case 2:
		FileName += "Lax Wendroff_";
		Scheme += "Lax Wendroff";
		break;
	case 3:
		FileName += "Beam Warming_";
		Scheme += "Beam Warming";
		break;
	default:
		FileName += "TypeError_";
	}
	FileName += std::to_string(OFP);
	FileName += "_";
	FileName += index;
	FileName += ".plt";
	file.open(FileName, std::ios::trunc);
	if (file.is_open())
	{
		std::cout << "output file open\n";
		file << "Variables = X, Q\n";
		file << "Zone t = \"" << Scheme << "\", i=" << m_cell_X.size() - 2*m_nSideGhost << ", f=point\n";
		for (int i = m_nSideGhost; i < m_cell_X.size() - m_nSideGhost; ++i)
		{
			file << m_cell_X[i] << "\t" << m_cell_Q[i] << "\n";
		}
		file.close();
		return 1;
	}
	else
	{
		std::cout << "output file error\n";
		return 0;
	}
}

void ControlArea::initialize(double area, Type initial)
{
	// calculate x step size
	m_step = area / double(m_cell_X.size() - 2 * m_nSideGhost);

	// assign cell coordinate
	for (int i = 0; i < m_cell_X.size(); ++i)
	{
		m_cell_X[i] = -(0.5*area) + (i - m_nSideGhost + 0.5)*m_step;
	}

	// assgin initial condition
	switch (initial)
	{
		// shock
	case 0:
		initial_condition_0(m_step);
		break;

		// expansion
	case 1:
		initial_condition_1(m_step);
		break;

		// sine
	case 2:
		initial_condition_2(m_step);
		break;

	default:
		std::cout << "initial condition error\n";
	}

	// save boundary values
	m_boundary_value_begin = m_cell_Q[m_nSideGhost];
	m_boundary_value_end = m_cell_Q[m_cell_X.size() - 1 - m_nSideGhost];
}

void ControlArea::initial_condition_0(double h)
{
	for (int i = 0; i < m_cell_Q.size(); ++i)
	{
		m_cell_Q[i] = 0.5*(0.2369268851*function_0(m_cell_X[i] - 0.5*h*0.9061798459) + 0.4786286705*function_0(m_cell_X[i] - 0.5*h*0.5384693101)
			+ 0.5688888889*function_0(m_cell_X[i]) + 0.4786286705*function_0(m_cell_X[i] + 0.5*h*0.5384693101) + 0.2369268851*function_0(m_cell_X[i] + 0.5*h*0.9061798459));
	}
	m_cell_Q_previous = m_cell_Q;
}

void ControlArea::initial_condition_1(double h)
{
	for (int i = 0; i < m_cell_Q.size(); ++i)
	{
		m_cell_Q[i] = 0.5*(0.2369268851*function_1(m_cell_X[i] - 0.5*h*0.9061798459) + 0.4786286705*function_1(m_cell_X[i] - 0.5*h*0.5384693101)
			+ 0.5688888889*function_1(m_cell_X[i]) + 0.4786286705*function_1(m_cell_X[i] + 0.5*h*0.5384693101) + 0.2369268851*function_1(m_cell_X[i] + 0.5*h*0.9061798459));
	}
	m_cell_Q_previous = m_cell_Q;
}

void ControlArea::initial_condition_2(double h)
{
	for (int i = 0; i < m_cell_Q.size(); ++i)
	{
		m_cell_Q[i] = 0.5*(0.2369268851*function_2(m_cell_X[i] - 0.5*h*0.9061798459) + 0.4786286705*function_2(m_cell_X[i] - 0.5*h*0.5384693101)
			+ 0.5688888889*function_2(m_cell_X[i]) + 0.4786286705*function_2(m_cell_X[i] + 0.5*h*0.5384693101) + 0.2369268851*function_2(m_cell_X[i] + 0.5*h*0.9061798459));
	}
	m_cell_Q_previous = m_cell_Q;
}

double ControlArea::function_0(double x)
{
	double out = 0.0;

	if (x < 0.0) out = 2.0;
	else out = 0.0;

	return out;
}

double ControlArea::function_1(double x)
{
	double out = 0.0;

	if (x < 0.0) out = -1.0;
	else out = 1.0;
	
	return out;
}

double ControlArea::function_2(double x)
{
	double out = 1.0 + 0.5*sin(2*M_PI*x);

	return out;
}