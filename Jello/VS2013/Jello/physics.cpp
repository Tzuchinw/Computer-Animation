/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1
  Tzu-Chin Wang (Rin)

*/

#include "jello.h"
#include "physics.h"

#define cubeL (double)1.00/(double)7.00

#define cubeL3 cubeL * sqrt((double)3)	//length for 3D diagonal line
#define cubeL2 cubeL * sqrt((double)2)	//length for 2D diagonal line

//	Calculate structure force
inline void structForce(struct world * jello, const double &ke, const double &de, const int &i, const int &j, const int &k, point *result)
{
	point tmp;
	point dx, dv;

	double stringL, stringLs, tempF, tempD;

	//i-axis
	if(i != 0)
	{
		dx = jello->p[i-1][j][k] - jello->p[i][j][k];
		dv = jello->v[i-1][j][k] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;

		if(i != 7)
		{
			dx = jello->p[i+1][j][k] - jello->p[i][j][k];
			dv = jello->v[i+1][j][k] - jello->v[i][j][k];
		
			stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
			stringL = sqrt(stringLs);

			tempF = ke * (stringL - cubeL);
			tmp += dx * tempF / stringL;

			tempD = de / stringLs * dv.dot(dx);
			tmp += dx * tempD;
		}
	}
	else  //i==0
	{
		dx = jello->p[i+1][j][k] - jello->p[i][j][k];
		dv = jello->v[i+1][j][k] - jello->v[i][j][k];
		
		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}

	//j-axis
	if(j != 0)
	{
		dx = jello->p[i][j-1][k] - jello->p[i][j][k];
		dv = jello->v[i][j-1][k] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;

		if(j != 7)
		{
			dx = jello->p[i][j+1][k] - jello->p[i][j][k];
			dv = jello->v[i][j+1][k] - jello->v[i][j][k];

			stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
			stringL = sqrt(stringLs);

			tempF = ke * (stringL - cubeL);
			tmp += dx * tempF / stringL;

			tempD = de / stringLs * dv.dot(dx);
			tmp += dx * tempD;
		}
	}
	else  //j==0
	{
		dx = jello->p[i][j+1][k] - jello->p[i][j][k];
		dv = jello->v[i][j+1][k] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}

	//k-axis
	if(k != 0)
	{
		dx = jello->p[i][j][k-1] - jello->p[i][j][k];
		dv = jello->v[i][j][k-1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;

		if(k != 7)
		{
			dx = jello->p[i][j][k+1] - jello->p[i][j][k];
			dv = jello->v[i][j][k+1] - jello->v[i][j][k];
		
			stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
			stringL = sqrt(stringLs);

			tempF = ke * (stringL - cubeL);
			tmp += dx * tempF / stringL;

			tempD = de / stringLs * dv.dot(dx);
			tmp += dx * tempD;
		}
	}
	else  //k==0
	{
		dx = jello->p[i][j][k+1] - jello->p[i][j][k];
		dv = jello->v[i][j][k+1] - jello->v[i][j][k];
		
		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}

	*result += tmp;
}

//	Calculate bend force
inline void bendForce(struct world * jello, const double &ke, const double &de, const int &i, const int &j, const int &k, point *result)
{
	point tmp;
	point dx, dv;

	double stringL, stringLs, tempF, tempD;

	//i-axis
	if(i > 1)
	{
		dx = jello->p[i-2][j][k] - jello->p[i][j][k];
		dv = jello->v[i-2][j][k] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL * 2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;

		if(i < 6)
		{
			dx = jello->p[i+2][j][k] - jello->p[i][j][k];
			dv = jello->v[i+2][j][k] - jello->v[i][j][k];
		
			stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
			stringL = sqrt(stringLs);

			tempF = ke * (stringL - cubeL * 2);
			tmp += dx * tempF / stringL;

			tempD = de / stringLs * dv.dot(dx);
			tmp += dx * tempD;
		}
	}
	else  //i==0
	{
		dx = jello->p[i+2][j][k] - jello->p[i][j][k];
		dv = jello->v[i+2][j][k] - jello->v[i][j][k];
		
		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL * 2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	
	//j-axis
	if(j > 1)
	{
		dx = jello->p[i][j-2][k] - jello->p[i][j][k];
		dv = jello->v[i][j-2][k] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL * 2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;

		if(j < 6)
		{
			dx = jello->p[i][j+2][k] - jello->p[i][j][k];
			dv = jello->v[i][j+2][k] - jello->v[i][j][k];

			stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
			stringL = sqrt(stringLs);

			tempF = ke * (stringL - cubeL * 2);
			tmp += dx * tempF / stringL;

			tempD = de / stringLs * dv.dot(dx);
			tmp += dx * tempD;
		}
	}
	else  //j==0
	{
		dx = jello->p[i][j+2][k] - jello->p[i][j][k];
		dv = jello->v[i][j+2][k] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL * 2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}

	//k-axis
	if(k > 1)
	{
		dx = jello->p[i][j][k-2] - jello->p[i][j][k];
		dv = jello->v[i][j][k-2] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL * 2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;

		if(k < 6)
		{
			dx = jello->p[i][j][k+2] - jello->p[i][j][k];
			dv = jello->v[i][j][k+2] - jello->v[i][j][k];
		
			stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
			stringL = sqrt(stringLs);

			tempF = ke * (stringL - cubeL * 2);
			tmp += dx * tempF / stringL;

			tempD = de / stringLs * dv.dot(dx);
			tmp += dx * tempD;
		}
	}
	else  //k==0
	{
		dx = jello->p[i][j][k+2] - jello->p[i][j][k];
		dv = jello->v[i][j][k+2] - jello->v[i][j][k];
		
		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL * 2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	
	*result += tmp;
}

//	Calculate shear force
inline void shearForce(struct world * jello, const double &ke, const double &de, const int &i, const int &j, const int &k, point *result)
{
	point tmp;
	point dx, dv;
	double stringL, stringLs, tempF, tempD;
	
	//8
	if(i!=7&&j!=7&&k!=7)
	{
		dx = jello->p[i+1][j+1][k+1] - jello->p[i][j][k];
		dv = jello->v[i+1][j+1][k+1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL3);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	if(i!=7&&j!=7&&k!=0)
	{
		dx = jello->p[i+1][j+1][k-1] - jello->p[i][j][k];
		dv = jello->v[i+1][j+1][k-1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL3);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	if(i!=7&&j!=0&&k!=7)
	{
		dx = jello->p[i+1][j-1][k+1] - jello->p[i][j][k];
		dv = jello->v[i+1][j-1][k+1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL3);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	if(i!=7&&j!=0&&k!=0)
	{
		dx = jello->p[i+1][j-1][k-1] - jello->p[i][j][k];
		dv = jello->v[i+1][j-1][k-1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL3);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	if(i!=0&&j!=7&&k!=7)
	{
		dx = jello->p[i-1][j+1][k+1] - jello->p[i][j][k];
		dv = jello->v[i-1][j+1][k+1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL3);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	if(i!=0&&j!=7&&k!=0)
	{
		dx = jello->p[i-1][j+1][k-1] - jello->p[i][j][k];
		dv = jello->v[i-1][j+1][k-1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL3);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	if(i!=0&&j!=0&&k!=7)
	{
		dx = jello->p[i-1][j-1][k+1] - jello->p[i][j][k];
		dv = jello->v[i-1][j-1][k+1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL3);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	if(i!=0&&j!=0&&k!=0)
	{
		dx = jello->p[i-1][j-1][k-1] - jello->p[i][j][k];
		dv = jello->v[i-1][j-1][k-1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL3);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}

	//12
	//fix
	if(j!=0&&k!=0)
	{
		dx = jello->p[i][j-1][k-1] - jello->p[i][j][k];
		dv = jello->v[i][j-1][k-1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	if(j!=0&&k!=7)
	{
		dx = jello->p[i][j-1][k+1] - jello->p[i][j][k];
		dv = jello->v[i][j-1][k+1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	if(j!=7&&k!=0)
	{
		dx = jello->p[i][j+1][k-1] - jello->p[i][j][k];
		dv = jello->v[i][j+1][k-1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	if(j!=7&&k!=7)
	{
		dx = jello->p[i][j+1][k+1] - jello->p[i][j][k];
		dv = jello->v[i][j+1][k+1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	//fix
	if(i!=0&&k!=0)
	{
		dx = jello->p[i-1][j][k-1] - jello->p[i][j][k];
		dv = jello->v[i-1][j][k-1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	if(i!=0&&k!=7)
	{
		dx = jello->p[i-1][j][k+1] - jello->p[i][j][k];
		dv = jello->v[i-1][j][k+1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	if(i!=7&&k!=0)
	{
		dx = jello->p[i+1][j][k-1] - jello->p[i][j][k];
		dv = jello->v[i+1][j][k-1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	if(i!=7&&k!=7)
	{
		dx = jello->p[i+1][j][k+1] - jello->p[i][j][k];
		dv = jello->v[i+1][j][k+1] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	//fix
	if(i!=0&&j!=0)
	{
		dx = jello->p[i-1][j-1][k] - jello->p[i][j][k];
		dv = jello->v[i-1][j-1][k] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	if(i!=0&&j!=7)
	{
		dx = jello->p[i-1][j+1][k] - jello->p[i][j][k];
		dv = jello->v[i-1][j+1][k] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	if(i!=7&&j!=0)
	{
		dx = jello->p[i+1][j-1][k] - jello->p[i][j][k];
		dv = jello->v[i+1][j-1][k] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	if(i!=7&&j!=7)
	{
		dx = jello->p[i+1][j+1][k] - jello->p[i][j][k];
		dv = jello->v[i+1][j+1][k] - jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = ke * (stringL - cubeL2);
		tmp += dx * tempF / stringL;

		tempD = de / stringLs * dv.dot(dx);
		tmp += dx * tempD;
	}
	
	*result += tmp;
}

//	Calculate force of force field
inline void fieldForce(struct world * jello, int &i, int &j, int &k, point *result)
{
	int resol = jello->resolution;
	double ffBlockL = 4.0 / (resol-1);

	//Obtain position(world space) of the point
	point pos = jello->p[i][j][k];
	pos.x += 2;
	pos.y += 2;
	pos.z += 2;

	//Convert position to corresponding force field sample point(least position value)
	int fi = (pos.x) / ffBlockL;
	int fj = (pos.y) / ffBlockL;
	int fk = (pos.z) / ffBlockL;
	//Next sample point(max position value)
	int finxt;
	int fjnxt;
	int fknxt;

	//boundary check
	if(fi < (resol-1))
	{
		if(fi < 0)	fi = 0;
		finxt = fi+1;
	}
	else
		finxt = fi = (resol-1);
		
	if(fj < (resol-1))
	{
		if(fj < 0)	fj = 0;
		fjnxt = fj+1;
	}
	else
		fjnxt = fj = (resol-1);

	if(fk < (resol-1))
	{
		if(fk < 0)	fk = 0;
		fknxt = fk+1;
	}
	else
		fknxt = fk = (resol-1);

	//Trilinear interpolation along j k i
	//interpolate along j 4
	point intj1, intj2, intj3, intj4;

	intj1 = jello->forceField[fi * resol * resol + fj * resol + fk] * (pos.y - fj * ffBlockL) +
			jello->forceField[fi * resol * resol + fjnxt * resol + fk] * (fjnxt * ffBlockL - pos.y);
	intj1 /= ffBlockL;

	intj2 = jello->forceField[fi * resol * resol + fj * resol + fknxt] * (pos.y - fj * ffBlockL) +
			jello->forceField[fi * resol * resol + fjnxt * resol + fknxt] * (fjnxt * ffBlockL - pos.y);
	intj2 /= ffBlockL;

	intj3 = jello->forceField[finxt * resol * resol + fj * resol + fk] * (pos.y - fj * ffBlockL) +
			jello->forceField[finxt * resol * resol + fjnxt * resol + fk] * (fjnxt * ffBlockL - pos.y);
	intj3 /= ffBlockL;

	intj4 = jello->forceField[finxt * resol * resol + fj * resol + fknxt] * (pos.y - fj * ffBlockL) +
			jello->forceField[finxt * resol * resol + fjnxt * resol + fknxt] * (fjnxt * ffBlockL - pos.y);
	intj4 /= ffBlockL;

	//interpolate along k 2
	intj1 = intj1 * (pos.z - fk * ffBlockL) +
			intj2 * (fknxt * ffBlockL - pos.z);
	intj1 /= ffBlockL;

	intj3 = intj3 * (pos.z - fk * ffBlockL) +
			intj4 * (fknxt * ffBlockL - pos.z);
	intj3 /= ffBlockL;

	//interpolate along i
	intj1 = intj1 * (pos.x - fi * ffBlockL) +
			intj3 * (finxt * ffBlockL - pos.x);
	intj1 /= ffBlockL;

	*result += intj1;
}

//	collision string force (against bounding cube)
inline void collision0(struct world * jello, const double &kc, const double &dc, int &i, int &j, int &k, point *result)
{
	point tmp;
	point dx, dv;
	double stringL, stringLs, tempF, tempD;

	if(jello->p[i][j][k].x > 2)
	{	
		point wall(jello->p[i][j][k]);
		wall.x = 2;
			
		dx = wall - jello->p[i][j][k];
		dv = jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = kc * (stringL);
		tmp += dx * tempF / stringL;

		tempD = dc / stringLs * dv.dot(dx);
		tmp -= dx * tempD;
	}
	else if(jello->p[i][j][k].x < -2)	
	{
		point wall(jello->p[i][j][k]);
		wall.x = -2;
			
		dx = wall - jello->p[i][j][k];
		dv = jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = kc * (stringL);
		tmp += dx * tempF / stringL;

		tempD = dc / stringLs * dv.dot(dx);
		tmp -= dx * tempD;
	}

	if(jello->p[i][j][k].y > 2)
	{	
		point wall(jello->p[i][j][k]);
		wall.y = 2;
			
		dx = wall - jello->p[i][j][k];
		dv = jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = kc * (stringL);
		tmp += dx * tempF / stringL;

		tempD = dc / stringLs * dv.dot(dx);
		tmp -= dx * tempD;
	}
	else if(jello->p[i][j][k].y < -2)	
	{
		point wall(jello->p[i][j][k]);
		wall.y = -2;
			
		dx = wall - jello->p[i][j][k];
		dv = jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = kc * (stringL);
		tmp += dx * tempF / stringL;

		tempD = dc / stringLs * dv.dot(dx);
		tmp -= dx * tempD;
	}
	if(jello->p[i][j][k].z > 2)
	{	
		point wall(jello->p[i][j][k]);
		wall.z = 2;
			
		dx = wall - jello->p[i][j][k];
		dv = jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = kc * (stringL);
		tmp += dx * tempF / stringL;

		tempD = dc / stringLs * dv.dot(dx);
		tmp -= dx * tempD;
	}
	else if(jello->p[i][j][k].z < -2)	
	{
		point wall(jello->p[i][j][k]);
		wall.z = -2;
			
		dx = wall - jello->p[i][j][k];
		dv = jello->v[i][j][k];

		stringLs = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
		stringL = sqrt(stringLs);

		tempF = kc * (stringL);
		tmp += dx * tempF / stringL;

		tempD = dc / stringLs * dv.dot(dx);
		tmp -= dx * tempD;
	}

	*result += tmp;
}

//	position based collision (against bounding cube)
void collision(struct world * jello, int &i, int &j, int &k)
{
	if(jello->p[i][j][k].x > 2)
	{	
		jello->v[i][j][k].x = -(jello->p[i][j][k].x - 2)/jello->dt;
		jello->p[i][j][k].x = 2;
	}
	else if(jello->p[i][j][k].x < -2)	
	{
		jello->v[i][j][k].x = -(jello->p[i][j][k].x + 2)/jello->dt;
		jello->p[i][j][k].x = -2;
	}
	if(jello->p[i][j][k].y > 2)
	{	
		jello->v[i][j][k].y = -(jello->p[i][j][k].y - 2)/jello->dt;
		jello->p[i][j][k].y = 2;
	}
	else if(jello->p[i][j][k].y < -2)	
	{
		jello->v[i][j][k].y = -(jello->p[i][j][k].y + 2)/jello->dt;
		jello->p[i][j][k].y = -2;
	}
	if(jello->p[i][j][k].z > 2)
	{	
		jello->v[i][j][k].z = -(jello->p[i][j][k].z - 2)/jello->dt;
		jello->p[i][j][k].z = 2;
	}
	else if(jello->p[i][j][k].z < -2)	
	{
		jello->v[i][j][k].z = -(jello->p[i][j][k].z + 2)/jello->dt;
		jello->p[i][j][k].z = -2;
	}
}
//	position based collision (against Ball)
void collisionBall(struct world * jello, int &i, int &j, int &k, struct ball * sphere)
{
	point dx = jello->p[i][j][k] - sphere->p;
	double lengths = (dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
	double length = sqrt(lengths);

	if(length < sphere->r)
	{
		jello->p[i][j][k] = sphere->p + dx/length * sphere->r;
		jello->v[i][j][k] = dx/length * sphere->r/jello->dt;
	}

}
//	position based collision (against inclined plane)
void collisionPlane(struct world * jello, int &i, int &j, int &k)
{
	point n(jello->a,jello->b,jello->c);
	
	point po(0,0,-jello->d/jello->c);
	point vector = jello->p[i][j][k] - po;
	n.normalize();

	double pointv = n.dot(jello->p[i][j][k]) + jello->d;
	double planev = n.dot(n) + jello->d;

	double distance =  n.dot(vector);
	if(distance < 0)	//check if it's in another side of plane
	{
		jello->p[i][j][k] = jello->p[i][j][k] - n * distance;
		jello->v[i][j][k] = n * distance / jello->dt * -1;
	}
}


/*
//	failed
inline void dragforce(struct world * jello, const double &kc, const double &dc, int &i, int &j, int &k, point *result)
{
	//F = 1/2 * (density of the fluid) * v^2 * (cross sectional area) * drag coefficient
	point tmp;
	point vel = jello->v[i][j][k];
	double den = 1.0;
	double coef = 0.001;

	point velsqr(vel.x *vel.x,vel.x *vel.x,vel.x *vel.x);
	tmp = point(0,0,0)-velsqr*den*coef;
	
	*result += tmp;
}
*/
/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
	// I did this!
	double mass = jello->mass;
	double ke = jello->kElastic;
	double de = jello->dElastic;
	double kc = jello->kCollision;
	double dc = jello->dCollision;
	
	for(int i = 0; i < 8 ; i++)
	{
		for(int j = 0; j < 8 ; j++)
		{
			for(int k = 0; k < 8 ; k++)
			{
				point ftmp;

				structForce(jello, ke, de, i, j, k, &ftmp);
				bendForce(jello, ke, de, i, j, k, &ftmp);
				shearForce(jello, ke, de, i, j, k, &ftmp);
				fieldForce(jello, i, j, k, &ftmp);

				if(!positionBasedCollision)	//check positionBasedCollision flag to decide which method to apply
					collision0(jello, kc, dc, i, j, k, &ftmp);

				//bonus 1: pushing
				if(g_iLeftMouseButton)
				{
					ftmp -= point(cos(Phi) * cos (Theta), sin(Phi) * cos (Theta), sin (Theta));
					/*
						glReadBuffer(GL_FRONT);
						glReadPixels(g_vMousePos[0], g_vMousePos[1], 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth_z);
					*/
				}

				a[i][j][k] = ftmp/mass;
			}
		}
	}
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}
