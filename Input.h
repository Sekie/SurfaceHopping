#include <Eigen/Dense>
#include <vector>
#include <string>

class InputObj
{
public:
	int NumSurfaces = 2;
	int NumSteps;
	int NumWalkers;
	double TimeStep;
	double Momentum;
	double Position;
	double Mass;
	std::string Case;
	void Init();
};