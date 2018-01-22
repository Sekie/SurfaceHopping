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
#include "Input.h"

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
	if (Case == "HO")
	{
		V = Eigen::MatrixXd::Zero(2, 2);
		V(0, 0) = 0.5 * x * x;
		V(1, 1) = 0.5 * x * x;
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
   
   To solve this, we must solve the equations
		y1(x + h) = y1(x) + 1/6 (k1 + 2b * k2 + 2d * k3 + k4)
		y2(x + h) = y2(x) + 1/6 (j1 + 2b * j2 + 2d * j3 + j4)
		k1 = h * f1(x, y1, y2)
		j1 = h * f2(...)
		k2 = h * f1(x + h / 2, y1 + k1 / 2, y2 + j1 / 2)
		j2 = h * f2(...)
		k3 = h * f1(x + h / 2, y1 + a * k1 + b * k2, y2 + a * j1 + b * j2)
		j3 = h * f2(...)
		k4 = h * f1(x + h, y1 + c * k2, y2 + d * k3)
		j4 = h * f2(...)
	which in terms of F = ma, replace the variables
		x  ---> t
		y1 ---> x
		y2 ---> v
	and 
		f1(t, x, v) = v
		f2(t, x, v) = F[x] / m
	so there is no explicit t dependence and we don't consider changes in t.
*/
void RKUpdate(double &x, double &v, int Surface, InputObj Input)
{
	// First, we define the constants in 4th order RK.
	double RKa = (sqrt(2.0) - 1.0) / 2.0;
	double RKb = (2.0 - sqrt(2.0)) / 2.0;
	double RKc = -1 * sqrt(2.0) / 2.0;
	double RKd = (2 + sqrt(2.0)) / 2.0;
	double k1, k2, k3, k4, j1, j2, j3, j4;

	k1 = Input.TimeStep * v;
	std::vector<double> F = CalcForce(x, Input.Case);
	j1 = Input.TimeStep * F[Surface] / Input.Mass;

	k2 = Input.TimeStep * (v + j1 / 2);
	F = CalcForce(x + k1 / 2, Input.Case);
	j2 = Input.TimeStep * F[Surface] / Input.Mass;

	k3 = Input.TimeStep * (v + RKa * j1 + RKb * j2);
	F = CalcForce(x + RKa * k1 + RKb * k2, Input.Case);
	j3 = Input.TimeStep * F[Surface] / Input.Mass;

	k4 = Input.TimeStep * (v + RKc * j2 + RKd * j3);
	F = CalcForce(x + RKc * k2 + RKd * k3, Input.Case);
	j4 = Input.TimeStep * F[Surface] / Input.Mass;

	x = x + (k1 + 2 * RKb * k2 + 2 * RKd * k3 + k4) / 6;
	v = v + (j1 + 2 * RKb * j2 + 2 * RKd * j3 + j4) / 6;
}