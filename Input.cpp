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

int Init()
{
	std::string Case;
	std::cout << "Enter case (A / B / C / D):" << std::endl;
	std::cin >> Case;

	return 0;
}