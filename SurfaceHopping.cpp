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
#include <complex>
#include "Input.h"

#define HBAR 1.0

void RKUpdate(double &x, double &v, int Surface, InputObj Input);
Eigen::MatrixXd Potential(double x, std::string Case);

/* Updates the density matrix using Eq (11). It is assumed that we are in the diabatic basis and the derivative coupling
   terms are zero. The update of interest is reduced to: 
		i * hbar * da / dt = Va - aV 
   and b is updated using
		b = 2 / hbar * Im[a^* cwiseProduct V] */
void UpdateA(Eigen::MatrixXcd &a, Eigen::MatrixXd &b, Eigen::MatrixXd V, double TimeStep)
{
	std::complex<double> IM_UNIT(0.0, 1.0);
	 Eigen::MatrixXcd da = (TimeStep / (HBAR * IM_UNIT)) * (V * a - a * V);
	 a = a + da;
	 b = (2.0 / HBAR) * ((a.conjugate()).cwiseProduct(V)).imag();
}

int SurfaceHopping(std::string Case, int NumWalkers)
{
	// Definitions for future generalizations.
	int NumSurface = 2;

	// Parameters related to surface hopping. The variables are named the same way as in Tully's publication.
	std::vector<double> aPrev;
	std::vector<double> a; // Diagonal Elements of a
	std::vector< std::vector< double > > b; // Matrix for b
	double TimeStep = 0.1;
	// Initialize the vectors.
	for (int i = 0; i < NumSurface; i++)
	{
		aPrev.push_back(0);
		a.push_back(0);

		std::vector<double> tmpVec;
		for (int j = 0; j < NumSurface; j++)
		{
			tmpVec.push_back(0.1); // Inialize to something not zero?
		}
		b.push_back(tmpVec);
	}

	return 0;
}

int main()
{
	InputObj InitParam;
	// InitParam.Init();
	InitParam.Default();

	std::ofstream Output(InitParam.OutputName.c_str());

	srand(0); // Seed for RNG.

	/* Classical molecular dynamics is ran for each walker, with hopping considered after each step. */
	for (int Walker = 0; Walker < InitParam.NumWalkers; Walker++) // Loops over each walker
	{
		/***** STEP 1 *****/
		/* The initial values are set including the momentum and density matrix */

		// MD Definitions
		int CurrentSurface = 0; // Which PES the walker is on.
		double Position = InitParam.Position; // Current position of the walker.
		double Velocity = InitParam.Momentum / InitParam.Mass; // Current velcity of the walker;
																	   
		// Surface Hopping Definitions
		Eigen::MatrixXcd a = Eigen::MatrixXcd::Zero(InitParam.NumSurfaces, InitParam.NumSurfaces); // This is the density matrix.
		a(0, 0) = 1; // Initialize the density matrix so that initially, only the first state is populated.
		Eigen::MatrixXd b = Eigen::MatrixXd::Zero(InitParam.NumSurfaces, InitParam.NumSurfaces); // b Matrix that controls transition probabilities.

		for (int Step = 0; Step < InitParam.NumSteps; Step++)
		{
			/***** STEP 2 *****/
			/* Position and velocity is propagated along the current PES classically. */
			// Run one step of MD.
			RKUpdate(Position, Velocity, CurrentSurface, InitParam);

			/***** STEP 3 *****/
			/* Calculate probability of surface hopping and attempt a surface hop. */

			// Update V using new position and use that to update a and b.
			Eigen::MatrixXd V; // = Eigen::MatrixXd::Zero(InitParam.NumSurfaces, InitParam.NumSurfaces);
			V = Potential(Position, InitParam.Case);
			UpdateA(a, b, V, InitParam.TimeStep);
			//std::cout << "a\n" << a << std::endl;
			//std::cout << "b\n" << b << std::endl;
			//std::cout << "x = " << Position << std::endl;

			// This is Equation (19) without the denominator.
			Eigen::MatrixXd g = InitParam.TimeStep * b.transpose();
			// We isolate one row of g, the one that matters for hopping considerations.
			Eigen::VectorXd gRow = g.row(CurrentSurface) / std::real(a(CurrentSurface, CurrentSurface));
			
			// Now we find a random number between 0 and 1.
			int z_int = rand();
			double z = (double)z_int / RAND_MAX;

			// And then we figure out which surface this corresponds to.
			int NewSurface = CurrentSurface;
			for (int i = 0; i < gRow.size(); i++)
			{
				if (gRow[i] < 0)
				{
					// We set negative elements to zero, which means z cannot be less than it and we can move on.
					continue;
				}

				if (z > gRow[i])
				{
					continue;
				}
				else
				{
					NewSurface = i;
				}
			}

			/***** STEP 4 *****/
			/* We rescale the energy if there is a surface hopping. This is done by simply varying the magnitude of velocity without
			   changing anything about the direction. I am not sure how to treat direction in the case of zero nonadibatic coupling */
			if (NewSurface != CurrentSurface)
			{
				std::cout << "Hopping attempt on step " << Step << " for walker " << Walker << ". Surface " << CurrentSurface << " to " << NewSurface << std::endl;
				double ChangeInEnergy = V(CurrentSurface, CurrentSurface) - V(NewSurface, NewSurface);
				if (ChangeInEnergy < 0) // If this is a decrease, we need to check that we have enough KE to make the change.
				{
					if (fabs(ChangeInEnergy) < 0.5 * InitParam.Mass * Velocity * Velocity) // Means we have enough KE to hop
					{
						std::cout << "Success" << std::endl;
						CurrentSurface = NewSurface;
						// Subtract from current velocity. The last term is to account for sign.
						Velocity = Velocity - (sqrt(2 * InitParam.Mass * fabs(ChangeInEnergy)) * Velocity / fabs(Velocity));
					} // Otherwise, nothing should happen.
				}
				else // Otherwise, make the change and rescale velocity.
				{
					std::cout << "Success" << std::endl;
					CurrentSurface = NewSurface;
					// Add to current velocity. The last term is to account for sign. Delta E should be positive, but just to be sure...
					Velocity = Velocity + (sqrt(2 * InitParam.Mass * fabs(ChangeInEnergy)) * Velocity / fabs(Velocity));
				}
			}
		}
	}
	system("pause");
	return 0;
}