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

#define MASS 1.0

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
std::vector<double> CalcForce(double x, std::string Case)
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

/* Implementation of the Gill's method of 4th order Runge-Kutta for two variables,
   solving Newton's Equation of motion 
		F = ma
		
		d  [ y2 ] = [ dv/dt ] = [ F / m ] = [ F(y1) / m ]
		dt [ y1 ] = [ dx/dt ] = [ p / m ] = [ y2 / m    ]
*/
void RKUpdate(double &x, double &p, double TimeStep, int Surface, std::string Case)
{
	// First, we define the constants in 4th order RK.
	double Mass = 1.0;

	double h = 0.5;
	double RKa = (sqrt(2.0) - 1.0) / 2.0;
	double RKb = (2.0 - sqrt(2.0)) / 2.0;
	double RKc = -1 * sqrt(2.0) / 2.0;
	double RKd = (2 + sqrt(2.0)) / 2.0;

	double k1, k2, k3, k4, j1, j2, j3, j4;

	k1 = h * p / MASS;
	std::vector<double> F = CalcForce(x, Case);
	j1 = h * F[Surface] / MASS;

	k2 = h * (p + 1 / 2 * j1) / MASS;
	F = CalcForce(x + k1 / 2, Case);
	j2 = h * F[Surface] / MASS;

	k3 = h * (p + RKa * j1 + RKb * j2) / MASS;
	F = CalcForce(x + RKa * k1 + RKb * k2, Case);
	j3 = h * F[Surface] / MASS;

	k4 = h * (p + RKc * j2 + RKd * j3) / MASS;
	F = CalcForce(x + RKc * k2 + RKd * k3, Case);
	j4 = h * F[Surface] / Mass;

	x = x + (k1 + 2 * RKb * k2 + 2 * RKd * k3 + k4) / 6;
	p = p + (j1 + 2 * RKb * j2 + 2 * RKd * j3 + j4) / 6;
}

void MDStep(double &x, double &p, double TimeStep, int Surface, std::string Case)
{

}