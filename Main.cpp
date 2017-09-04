#include <memory>
#include <iostream>

#include "DataType.h"
#include "RHS.h"
#include "Condition.h"
#include "ControlArea.h"
#include "Flux.h"
#include "ENO_Reconstruction.h"
#include "ENO_Flux.h"
#include "OrderTest.h"

#define condition Advection_Condition
#define controlArea Advection_ControlArea

int main()
{
	// initialization
	std::shared_ptr<Condition> Advection_Condition = std::make_shared<Condition>();

	if (condition->input() == 0) return 0;

	// Order test object
	std::shared_ptr<OrderTest> orderTest = std::make_shared<OrderTest>();
	
	// default ghost cell = 2
	std::shared_ptr<ControlArea> Advection_ControlArea = std::make_shared<ControlArea>(condition->getNFX(), 4);

	// inimializing cell quantities / e.g. averaging
	controlArea->initialize(condition->getArea(), condition->getInitial());

	double dx = controlArea->getStep();
	std::cout << "grid size = " << dx << "\n";

	controlArea->output("initial", condition->getInitial(), condition->getFlux(), condition->getOFP());

	// Set exact solution to order test
	orderTest->setExact(controlArea->getCell_Q());

	int NFX = condition->getNFX();
	int nSideGhost = controlArea->getSideGhost(); // number of ghost cell
	int OFP = condition->getOFP(); // order of polynomials

	NFX += 2 * nSideGhost; // number of grids

	// RHS calculater object
	std::shared_ptr<RHS> RHS_solver = std::make_shared<RHS>(condition, controlArea);

	// Reconstruction objects
	std::vector<ENO_Reconstruction> Rx;
	Rx.clear();

	for (int i = 0; i < NFX; ++i)
	{
		Rx.push_back(ENO_Reconstruction(i));
	}
	ENO_Reconstruction::setDX(dx); // set grid size
	ENO_Reconstruction::setOFP(OFP); // set order of polynomial
	ENO_Reconstruction::setX(controlArea->getCell_X()); // set X coordinate
	
	// time assign
	double time;
	double dt = condition->getTimeStep();
	double TargetTime = condition->getTargetTime();

	// advection flux structure
	std::shared_ptr<ENO_Flux> flux_p = std::make_shared<ENO_Flux>(forward, dx, dt, condition->getAdvSpeed(), controlArea->getCell_X());
	std::shared_ptr<ENO_Flux> flux_m = std::make_shared<ENO_Flux>(backward, dx, dt, condition->getAdvSpeed(), controlArea->getCell_X());
	
	// solution parameter
	std::vector<double> Solution;
	std::vector<double> prev_Solution;
	std::vector<double> shockSpeed;
	prev_Solution.resize(NFX);
	Solution.resize(NFX);
	shockSpeed.resize(NFX - 1);
	
	// time marching
	for (time = 0.0; time < TargetTime - dt;)
	{
		controlArea->BC(condition->getBC()); // boundary condition
		prev_Solution = controlArea->getCell_Q(); // previous solution variable

		// Calculate time step
		shockSpeed.clear();
		for (int_t i = 0; i < NFX - 1; ++i)
			shockSpeed.push_back(0.5*(prev_Solution[i] + prev_Solution[i + 1]));
		condition->update(shockSpeed);
		dt = condition->getTimeStep();
		
		ENO_Reconstruction::setQ(controlArea->getCell_Q()); // set average cell quantity
		ENO_Reconstruction::update(); // update primitive value

		for (int i = 0; i < Rx.size(); ++i)
		{
			Rx[i].set_stencil(); // update stencil
		}

		ENO_Flux::setReconstruction(Rx); // set reconstruction in flux object
		
		Solution = RHS_solver->cal_RHS(controlArea, dt);

		controlArea->setCellQ(Solution);

		time += dt;
		std::cout << "current time = " << time << "\n";
	}
	if (time < TargetTime) // last step
	{
		dt = TargetTime - time;
		flux_p->setDt(dt);
		flux_m->setDt(dt);

		controlArea->BC(condition->getBC()); // boundary condition
		prev_Solution = controlArea->getCell_Q(); // previous solution variable
		
		ENO_Reconstruction::setQ(controlArea->getCell_Q()); // set average cell quantity
		ENO_Reconstruction::update(); // update primitive value

		for (int i = 0; i < Rx.size(); ++i)
		{
			Rx[i].set_stencil(); // update stencil
		}

		ENO_Flux::setReconstruction(Rx); // set reconstruction in flux object

		Solution = RHS_solver->cal_RHS(controlArea, dt);

		controlArea->setCellQ(Solution);

		time += dt;
		std::cout << "Last time step = " << dt << "\n";
	}

	std::cout << "Finished in <" << time << ">(sec)\n";

	// Exporting result
	controlArea->output("result", condition->getInitial(), condition->getFlux(), condition->getOFP());

	// Calculate error
	double L1 = orderTest->L1error(controlArea->getCell_Q());
	double L2 = orderTest->L2error(controlArea->getCell_Q());
	double Linf = orderTest->Linf_error(controlArea->getCell_Q());

	std::cout << "L1 error    = " << L1 << "\n";
	std::cout << "L2 error    = " << L2 << "\n";
	std::cout << "L inf error = " << Linf << "\n";
	
	return 0;
}