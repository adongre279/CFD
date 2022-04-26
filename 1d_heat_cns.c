////////////////////////////////////////////////////////////////////////////////
//
// 1d_heat_cns.c
//
// Description :
// C program to solve 1D heat conduction equation.
// Using the Crank-Nicolson method
// Tridag subroutine used to solve matrix system referred from the book 'Numerical Recipes in C'
// TDMA subroutine referred from the book 'Numerical Heat Transfer and Fluid Flow' by Dr. Ghoshdatidar
//
// Metrics:
// Date/Time           : 12-09-2020 14:28:36
// Hostname            : colossus.it.mtu.edu
// OS                  : RHEL Workstation release 7.7 (Maipo)
// Kernel              : 3.10.0-1062.el7.x86_64
// RAM                 : 502 GB
// CPU model           : Intel(R) Xeon(R) Gold 6130 CPU @ 2.10GHz
// CPU/Core count      : 64
// Author information  : Akshay Dongre (adongre@mtu.edu)
// Source code license : GPL v3 (https://www.gnu.org/licenses/gpl-3.0.txt)
// Software/Language   : Dev-C++/ Linux C Compiler
// Version             : 5.11
// Pre-req/Dependency  :
// Compilation command : gcc -Wall -g -lm 1d_heat_cns.c -o 1d_heat_cns.x
// Compilation time    : A few seconds
// Execution command   : ./1d_heat_cns.x
//
// Note                :
// Execution time      : Apprx. 0.22 seconds
//
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
//#include <curses.h>     // For compiling on Linux machines
#include <conio.h> // Comment this line while compiling on Linux machine
#include <string.h>
#include <math.h>

float tridag(float a[], float b[], float c[], float r[], float u[], unsigned long n);
float TDMA(int Ii, int Nn, float Aa[], float Bb[], float Cc[], float Rr[], float Uu[]);

#define pi 3.141592653589793

#define alpha 0.02 // Thermal diffusivity (m^2/hr)
#define c0 100.0   // Initial temperature (Celcius)
#define L 1.0	   // Length of the rod (m)

int main()
{
	float dx, dt, t;
	float D;
	int nx, i, j, k;
	int tl;

	// User inputs

	printf("Enter the spatial step increment dx (eg. 0.1): ");
	scanf("%f", &dx);

	printf("Enter time step dt (eg. 0.1): ");
	scanf("%f", &dt);

	printf("Enter final time in hrs. (eg. 10.0): ");
	scanf("%f", &t);

	// Problem formulation

	D = alpha * (dt / (dx * dx));

	nx = round(L / dx); // Total number of nodes

	tl = round(t / dt); // Number of loops to run the code => (final time / time step)

	printf("\n");
	printf("D = %0.2f, nx = %d, tl = %d\n\n", D, nx, tl);

	float T0[nx + 1], T[nx + 1], rhs[nx + 1];
	float x_ax[nx + 1];
	float T_ana[nx + 1];
	float subdia[nx + 1], dia[nx + 1], supdiag[nx + 1], rgn[nx + 1], Temp[nx + 1];

	for (i = 1; i <= (nx + 1); i++)
	{
		x_ax[i] = (i - 1) * dx;

		// Solution initialization
		T0[i] = c0 * sin(pi * x_ax[i] / L);

		// Analytical solution
		T_ana[i] = c0 * exp((-alpha * pi * pi * t) / (L * L)) * sin(pi * x_ax[i] / L);
	}

	// Crank-Nicolson  method

	// Forming the sub-diagonal matrix

	for (i = 2; i <= nx; i++)
	{
		subdia[i] = -D / 2;
	}

	subdia[nx + 1] = 0.0;

	// Forming the super-diagonal matrix

	supdiag[1] = 0.0;

	for (i = 2; i <= nx; i++)
	{
		supdiag[i] = -D / 2;
	}

	// Forming the diagonal matrix

	dia[1] = 1.0;

	for (i = 2; i <= nx; i++)
	{
		dia[i] = 1 + D;
	}

	dia[nx + 1] = 1.0;

	// Solving the tri-diagonal matrix system using the TDMA algorithm

	for (j = 1; j <= tl; j++) // Time stepping
	{
		// Computing the RHS for Crank-Nicolson Method
		for (i = 2; i <= nx; i++)
		{
			rhs[i] = ((D / 2) * T0[i - 1]) + ((1 - D) * T0[i]) + ((D / 2) * T0[i + 1]);
		}

		// Solving the LHS for Crank-Nicolson Method using the TDMA algorithm
		tridag(subdia, dia, supdiag, rhs, T, nx);
		// TDMA(1, nx + 1, subdia, dia, supdiag, rhs, T);

		// Updating the RHS for time stepping
		for (i = 2; i <= nx; i++)
		{
			T0[i] = T[i];
		}
	}

	for (i = 1; i <= nx + 1; i++)
	{
		printf("i = %d, j = %d, D = %0.2f, tl = %d, dx = %0.2f, x_ax = %0.2f, T_ana = %0.2f, T = %0.2f\n", i, j, D, tl, dx, x_ax[i], T_ana[i], T[i]);
	}

	// Writing output data file

	FILE *fp;
	char fname[80];

	sprintf(fname, "phw1_cns_nodes_%d.dat", (nx + 1));

	fp = fopen(fname, "w");

	for (i = 1; i <= nx + 1; i++)
	{
		fprintf(fp, "%f,%f\n", x_ax[i], T[i]);
	}

	fclose(fp);

	return 0;
}

// tridag(a, b, c, r, u, n) => Tridiagonal matrix system solver algorithm that takes following inputs;
// a -> Sub-diagonal element
// b -> Diagonal element
// c -> Super-diagonal element
// r -> Right-hand side element array
// u -> Solution vector
// n -> Total number of equations
// Reference => Numerical Recipes in C

float tridag(float a[], float b[], float c[], float r[], float u[], unsigned long n)
{
	unsigned long j;
	float bet, gam[n];

	if (b[1] == 0.0)
		printf("Error 1 in tridag\n");

	bet = b[1];

	u[1] = r[1] / b[1];

	// Decomposition and forward substitution

	for (j = 2; j <= n; j++)
	{
		gam[j] = c[j - 1] / bet;
		bet = b[j] - a[j] * gam[j];

		if (bet == 0.0)
			printf("Error 2 in tridag\n");

		u[j] = (r[j] - a[j] * u[j - 1]) / bet;
	}

	// Back substitution

	for (j = (n - 1); j >= 1; j--)
		u[j] -= gam[j + 1] * u[j + 1];
}

// TDMA(Ii, Nn, Aa, Bb, Cc, Rr, Uu) => Tridiagonal matrix system solver algorithm that takes following inputs;
// Ii -> Initial number of equation
// Nn -> Total number of nodes
// Aa -> Sub-diagonal element
// Bb -> Diagonal element
// Cc -> Super-diagonal element
// Rr -> Right-hand side element array
// Uu -> Solution vector
// Reference => Numerical Heat Transfer and Fluid Flow by Dr. Ghoshdatidar

float TDMA(int Ii, int Nn, float Aa[], float Bb[], float Cc[], float Rr[], float Uu[])
{
	unsigned long j, k;
	float beta[Nn], gama[Nn];

	// Compute intermediate arrays beta and gama

	beta[1] = Bb[1];
	gama[1] = Rr[1];

	int i1 = Ii + 1;

	for (j = i1; j <= Nn; j++)
	{
		beta[j] = Bb[j] - ((Aa[j] * Cc[j - 1]) / beta[j - 1]);
		gama[j] = (Rr[j] - (Aa[j] * gama[j - 1])) / beta[j];
	}

	// Compute the solution vector Uu

	Uu[Nn] = gama[Nn];

	int n1 = Nn - 1;

	for (k = 1; k <= n1; k++)
	{
		j = Nn - k;

		Uu[j] = gama[j] - Cc[j] * Uu[j + 1] / beta[j];
	}
}