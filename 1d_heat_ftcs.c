////////////////////////////////////////////////////////////////////////////////
//
// 1d_heat_ftcs.c
//
// Description :
// C program to solve 1D heat conduction equation.
// Using the forward time-central space (Euler Explicit) method
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
// Compilation command : gcc -Wall -g -lm 1d_heat_ftcs.c -o 1d_heat_ftcs.x
// Compilation time    : A few seconds
// Execution command   : ./1d_heat_ftcs.x
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

#define pi 3.141592653589793

#define alpha 0.02 // Thermal diffusivity (m^2/hr)
#define c0 100.0   // Initial temperature (Celcius)
#define L 1.0	   // Length of the rod (m)

int main()
{
	float dx, dt, t, temp;
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

	temp = t;

	D = alpha * (dt / (dx * dx));

	nx = round(L / dx); // Total number of nodes

	tl = round(t / dt); // Number of loops to run the code => (final time / time step)

	printf("\n");
	printf("D = %0.2f, nx = %d, t = %f, tl = %d\n\n", D, nx, t, tl);

	float T0[nx + 1], T[nx + 1];
	float x_ax[nx + 1], x_an[1002];
	float T_ana[nx + 1], Tan[1002];

	// Computing the analytical solution

	for (i = 1; i <= 1002; i++)
	{
		x_an[i] = (i - 1) * (dx / 100);

		// Analytical solution
		Tan[i] = c0 * exp((-alpha * pi * pi * temp) / (L * L)) * sin(pi * x_an[i] / L);
	}

	// Computing numerical solution

	for (i = 1; i <= (nx + 1); i++)
	{
		x_ax[i] = (i - 1) * dx;

		// Solution initialization
		T0[i] = c0 * sin(pi * x_ax[i] / L);

		// Analytical solution
		T_ana[i] = c0 * exp((-alpha * pi * pi * temp) / (L * L)) * sin(pi * x_ax[i] / L);
	}

	// Boundary Conditions

	T[1] = 0.0;
	T[nx + 1] = 0.0;

	// Euler Explicit Method

	for (j = 1; j <= tl; j++) // Time stepping
	{
		for (i = 2; i <= nx; i++)
		{
			T[i] = T0[i] + D * (T0[i - 1] - (2 * T0[i]) + T0[i + 1]);
		}

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

	FILE *fp, *fp1;
	char fname[80], fname1[80];

	sprintf(fname, "phw1_ftcs_nodes_%d.dat", (nx + 1));
	sprintf(fname1, "phw1_analytical.dat");

	fp = fopen(fname, "w");
	fp1 = fopen(fname1, "w");

	for (i = 1; i <= nx + 1; i++)
	{
		fprintf(fp, "%f,%f\n", x_ax[i], T[i]);
	}

	for (i = 1; i <= 1001; i++)
	{
		fprintf(fp1, "%f,%f\n", x_an[i], Tan[i]); // Writes the analytical solution for 1D heat equation
	}

	fclose(fp);
	fclose(fp1);

	return 0;
}