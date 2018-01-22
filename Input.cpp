/* The following program runs Tully Surface Hopping as described by Tully 
   (J. Chem. Phys. 93 (2), 1061 (1990)) for the  four test cases listed in 
   that publication. This was done as part of my reading group work for the
   Van Voorhis Group at MIT. */
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <string>
#include "Input.h"

void InputObj::Init()
{
	std::cout << "Enter case (A / B / C / D):" << std::endl;
	std::cin >> Case;
	std::cout << "Enter mass:" << std::endl;
	std::cin >> Mass;
	std::cout << "Enter time step:" << std::endl;
	std::cin >> TimeStep;
	std::cout << "Enter initial momentum: " << std::endl;
	std::cin >> Momentum;
	std::cout << "Enter initial position: " << std::endl;
	std::cin >> Position;
	std::cout << "Enter number of steps in the dynamics:" << std::endl;
	std::cin >> NumSteps;
	std::cout << "Enter number of walkers:" << std::endl;
	std::cin >> NumWalkers;
}