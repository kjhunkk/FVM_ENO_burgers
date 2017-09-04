#include "Condition.h"

Condition::Condition()
{
	m_speed = 0.0;
	m_dx = 0.0;
	m_dt = 0.0;
	m_CFL = 0.0;
	m_target_time = 0.0;
	m_area = 0.0;
	m_OFP = 0;
	m_NFX = 0;
	m_initial = 0;
	m_BC = 0;
}

Condition::~Condition()
{

}

bool Condition::input()
{
	// read input file
	if (readInput() == 0) return 0;

	// read option list files
	std::vector<std::string> fluxOption;
	if (readOption("flux", fluxOption) == 0) return 0;
	std::vector<std::string> initialOption;
	if (readOption("initial", initialOption) == 0) return 0;
	std::vector<std::string> BCOption;
	if (readOption("BC", BCOption) == 0) return 0;

	// input options
	inputType(m_flux, "flux condition", fluxOption);
	inputType(m_initial, "initial condition", initialOption);
	inputType(m_BC, "Boundary condition", BCOption);
	
	m_NFX = int(m_area / m_dx);
	m_dt = m_CFL*m_dx / abs(m_speed);
	printConditions();

	return 1;
}

void Condition::update(std::vector<double> Shockspeed)
{
	// calculate maximum shock speed
	double maxS = *std::max_element(Shockspeed.begin(), Shockspeed.end());
	double minS = *std::min_element(Shockspeed.begin(), Shockspeed.end());
	
	maxS = maxS > abs(minS) ? maxS : abs(minS);

	// update conditions
	m_dt = m_CFL*m_dx / maxS;
}

bool Condition::readInput()
{
	std::ifstream file;
	std::string buff;
	file.open("./input.inp");

	if (file.is_open())
	{
		std::cout << "----------input file open----------\n";
		file >> buff >> buff >> buff >> buff >> buff >> buff >> m_speed >> m_target_time >> m_CFL >> m_dx >> m_area >> m_OFP;
		file.close();
		return 1;
	}
	else
	{
		std::cout << "----------cannot find input file----------\n";
		file.close();
		return 0;
	}
}

bool Condition::readOption(std::string name, std::vector<std::string>& option)
{
	std::ifstream file;
	std::string buff;
	std::string filename = "./";
	filename += name;
	filename += ".inp";
	file.open(filename);

	if (file.is_open())
	{
		std::cout << "----------option file open----------\n";
		file >> buff;
		while (buff != "End")
		{
			file >> buff >> buff;
			option.push_back(buff);
			if (option.size() > 100) break;
		}
		option.pop_back();
		return 1;
	}
	else
	{
		std::cout << "----------cannot find option file----------\n";
		std::cout << "file name = " << filename << "\n";
		file.close();
		return 0;
	}
}

bool Condition::inputType(Type& type, std::string name, std::vector<std::string> option)
{
	std::cout << "<input " << name << " type>\n";
	std::cout << "--------------------\n";
	std::cout << "-------Options------\n";
	for (int i = 0; i < option.size(); ++i)
	{
		std::cout << "#" << i << " : " << option[i] << "\n";
	}
	std::cout << name << " option = ";
	std::cin >> type;
	if ((type >= 0) && (type < option.size()))
	{
		std::cout << "--------------------\n";
		return 1;
	}
	else
	{
		std::cout << "!!Option error!!\n";
		return 0;
	}
}

void Condition::printConditions() const
{
	std::cout << "----------Conditions----------\n";
	std::cout << "$ Grid Size           = " << m_dx << "\n";
	std::cout << "$ CFL number          = " << m_CFL << "\n";
	std::cout << "$ Time step           = " << m_dt << "\n";
	std::cout << "$ Target Time         = " << m_target_time << "\n";
	std::cout << "$ Area                = " << m_area << "\n";
	std::cout << "$ Order of Polynomial = " << m_OFP << "\n";
	std::cout << "$ Number of Grids     = " << m_NFX << "\n";
	std::cout << "$ flux type           = " << m_flux << "\n";
	std::cout << "$ initial condition   = " << m_initial << "\n";
	std::cout << "$ boundary condition  = " << m_BC << "\n";
}
