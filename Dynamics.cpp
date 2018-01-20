#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <string>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

// This function sets up the potential matrix. The first two elements of the vector are the first and second diabatic surfaces and
// the third element is the electronic coupling between these two surfaces.
Eigen::MatrixXd Potential(double x, std::string Case)
{
	double tmpDouble;
	Eigen::MatrixXd V;
	if (Case == "A")
	{
		V = Eigen::MatrixXd::Zero(2, 2);
		double A = 0.01;
		double B = 1.6;
		double C = 0.005;
		double D = 1.0;

		// V11
		if (x > 0)
		{
			tmpDouble = A * (1 - std::exp(-1 * B * x));
		}
		else
		{
			tmpDouble = -1 * A * (1 - std::exp(B * x));
		}
		V(0, 0) = tmpDouble;
		
		// V22
		tmpDouble = -1 * tmpDouble;
		V(1, 1) = tmpDouble;

		// V12
		tmpDouble = C * std::exp(-1 * D * x * x);
		V(0, 1) = tmpDouble;
		V(1, 0) = tmpDouble;
	}
	return V;
}

// This function calculates -dV/dx, the force at the given point x. The first element of the vector is the force on the lower surface
// and the second element is the force on the higher surface.
std::vector<double> CalcForce(double x, double TimeStep, std::string Case)
{
	std::vector<double> F;
	double dx = 0.01;
	// We will do dV/dx = V(x+dx) - V(x-dx) / 2dx
	Eigen::MatrixXd VPlusDX = Potential(x + dx, Case);
	Eigen::MatrixXd VMinusDX = Potential(x - dx, Case);
	
	double tmpDouble;
	for (int i = 0; i < VPlusDX.rows(); i++)
	{
		tmpDouble = -1 * (VPlusDX(i, i) - VMinusDX(i, i)) / (2 * dx);
		F.push_back(tmpDouble);
	}

	return F;
}

void MDStep(double &x, double &p, double TimeStep, int Surface, std::string Case)
{

}