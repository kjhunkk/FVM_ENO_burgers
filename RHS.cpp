#include "RHS.h"

RHS::RHS(std::shared_ptr<Condition> condition, std::shared_ptr<ControlArea> controlArea)
{
	m_dx = controlArea->getStep();
	m_nSideGhost = controlArea->getSideGhost();
	m_flux_p = std::make_shared<ENO_Flux>(forward, m_dx, condition->getTimeStep(), condition->getAdvSpeed(), controlArea->getCell_X());
	m_flux_m = std::make_shared<ENO_Flux>(backward, m_dx, condition->getTimeStep(), condition->getAdvSpeed(), controlArea->getCell_X());
	m_bc = condition->getBC();
	m_NFX = condition->getNFX() + 2 * m_nSideGhost;

	// initilize reconstruction object
	m_Rx.clear();

	for (int i = 0; i < m_NFX; ++i)
	{
		m_Rx.push_back(ENO_Reconstruction(i));
	}
	ENO_Reconstruction::setDX(m_dx); // set grid size
	ENO_Reconstruction::setOFP(condition->getBC()); // set order of polynomial
	ENO_Reconstruction::setX(controlArea->getCell_X()); // set X coordinate
}

RHS::~RHS()
{

}

std::vector<double> RHS::cal_RHS(std::shared_ptr<ControlArea> controlArea, double dt)
{
	std::vector<double> Q = controlArea->getCell_Q();
	int NFX = Q.size();
	m_dt = dt;
	m_flux_p->setDt(dt);
	m_flux_m->setDt(dt);

	// Temporary control area
	std::shared_ptr<ControlArea> field = std::make_shared<ControlArea>(*controlArea);

	// step 1(Q1)
	// Apply boundary condition
	field->BC(m_bc);

	// Calculate RHS
	m_L = solver(field->getCell_Q());

	m_Q1.clear();
	m_Q1 = Q;
	for (int i = m_nSideGhost; i < NFX - m_nSideGhost; ++i)
	{
		m_Q1[i] = Q[i] + dt*m_L[i];
	}
	field->setCellQ(m_Q1);

	// step 2(Q2)
	// Apply boundary condition
	field->BC(m_bc);

	// Calculate RHS
	m_L = solver(field->getCell_Q());

	m_Q2.clear();
	m_Q2 = Q;
	for (int i = m_nSideGhost; i < NFX - m_nSideGhost; ++i)
	{
		m_Q2[i] = 0.75*Q[i] + 0.25*(m_Q1[i] + dt*m_L[i]);
	}
	field->setCellQ(m_Q2);

	// step 3(Q next)
	// Apply boundary conditoin
	field->BC(m_bc);

	// Calculate RHS
	m_L = solver(field->getCell_Q());

	m_Q.clear();
	m_Q = Q;
	for (int i = m_nSideGhost; i < NFX - m_nSideGhost; ++i)
	{
		m_Q[i] = (Q[i] + 2.0*m_Q2[i] + 2.0*dt*m_L[i]) / 3.0;
	}

	return m_Q;
}

std::vector<double> RHS::solver(const std::vector<double> Q)
{
	ENO_Reconstruction::setQ(Q); // set average cell quantity
	ENO_Reconstruction::update(); // update primitive value

	for (int i = 0; i < m_Rx.size(); ++i)
	{
		m_Rx[i].set_stencil(); // update stencil
	}

	ENO_Flux::setReconstruction(m_Rx);


	std::vector<double> solution;
	solution.resize(m_NFX);
	solution = Q;

	// solution update
	for (int i = m_nSideGhost; i < (m_NFX - m_nSideGhost); ++i)
	{
		solution[i] = -((m_flux_p->flux(i) - m_flux_m->flux(i)) / m_dx);
	}

	return solution;
}