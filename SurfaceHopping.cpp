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

#define HBAR 1.0

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

//int main()
//{
//	double TimeStep = 0.1;
//	Eigen::MatrixXd V = Eigen::MatrixXd::Zero(2, 2);
//	V(0, 0) = 1;
//	V(1, 1) = -1;
//	V(0, 1) = 0.5;
//	V(1, 0) = 0.5;
//	Eigen::MatrixXcd a = Eigen::MatrixXcd::Zero(2, 2);
//	a(0, 0) = 1.0;
//	Eigen::MatrixXcd b = Eigen::MatrixXcd::Zero(2, 2);
//	
//	std::cout << "Before Update:\na\n" << a.real() << "\n" << a.imag() << "\nb\n" << b.real() << "\nV\n" << V << std::endl;
//	UpdateA(a, b, V, TimeStep);
//	std::cout << "After Update:\na\n" << a.real() << "\n" << a.imag() << "\nb\n" << b.real() << "\nV\n" << V << std::endl;
//	system("pause");
//
//	return 0;
//}

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