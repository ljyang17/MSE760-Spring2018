#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// create lattice points
void createLattice(int xdim, int ydim, int zdim,
                   double initialx, double initialy, double initialz,
                   double* arrayX, double* arrayY, double* arrayZ, double a_red) {

  // define coordinates of lattice points
  // i is x axis, j is y axis, k is z axis
  int i = 0;
  int j = 0;
  int k = 0;

  for (i = 0; i < xdim; i++) {
    for (j = 0; j < ydim; j++) {
      for (k = 0; k < zdim; k++) {
        int offset = k + zdim * j + ydim * zdim * i;
        arrayX[offset] = initialx + offset * a_red;
        arrayY[offset] = initialy + offset * a_red;
        arrayZ[offset] = initialz + offset * a_red;
      }
    }
  } // end of for loop
} // end of function

// output lattice points to text file
void outToFile(FILE* fp, int xdim, int ydim, int zdim,
               double* arrayX, double* arrayY, double* arrayZ,
               double sig, int atomID) {
   int i = 0;
   int j = 0;
   int k = 0;
  for (i = 0; i < xdim; i++) {
    for (j = 0; j < ydim; j++) {
      for (k = 0; k < zdim; k++) {
        int offset = k + zdim * j + ydim * zdim * i;

        fprintf(fp, "cell ID: %d, atom %d, x=%.6e, ", offset, atomID, sig * arrayX[offset]);
        fprintf(fp, "y=%.6e, ", sig * arrayY[offset]);
        fprintf(fp, "z=%.6e\n", sig * arrayZ[offset]);
      }
    }
  } // end of for loop
} // end of function

// merge two arrays into one
void mergeArray(int xdim, int ydim, int zdim,
                double* arrayX1, double* arrayY1, double* arrayZ1,
                double* arrayX2, double* arrayY2, double* arrayZ2,
                double* arrayX, double* arrayY, double* arrayZ) {

  int halfNumOfAtoms = xdim * ydim * zdim;
  int i = 0;

  for (i = 0; i < halfNumOfAtoms; i++) {
    arrayX[i] = arrayX1[i];
    arrayX[i + halfNumOfAtoms] = arrayX2[i];
    arrayY[i] = arrayY1[i];
    arrayY[i + halfNumOfAtoms] = arrayY2[i];
    arrayZ[i] = arrayZ1[i];
    arrayZ[i + halfNumOfAtoms] = arrayZ2[i];
  } // end of for loop
} // end of function

// calculate distance between 2 points
double distanceSquare(double x1, double y1, double z1,
                      double x2, double y2, double z2){
  return pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2);
}

// calculate LJ energy
double LennardJones(int xdim, int ydim, int zdim,
                    double* arrayX, double* arrayY, double* arrayZ){

  int numOfAtoms = 2 * xdim * ydim * zdim;
  double distSquare = 0.0;
  double LJenergy = 0.0;
  int i = 0;
  int j = 0;

  for (i = 0; i < numOfAtoms; i++) {
    for (j = i + 1; j < numOfAtoms; j++) {
      distSquare = distanceSquare(arrayX[i], arrayY[i],arrayZ[i], arrayX[j], arrayY[j],arrayZ[j]);
      LJenergy += 4 * (pow((1 / distSquare), 6) - pow((1 / distSquare), 3));
    }
  } // end of for loop
  return LJenergy;
} // end of function

// start main
int main(int argc, char const *argv[]) {

  // define base constants
  double sig = 3.4e-10;
  double eps = 1.04e-2;
  double mass = 6.6e-26;

  // define actual constants
  double a = 5.26e-10;
  double m = 6.6e-26;

  // define reduced-unit constants
  double a_red = a / sig;
  double m_red = m / mass;

  // define initial offset of Argon
  double b1x = 0.0;
  double b1y = 0.0;
  double b1z = 0.0;

  double b2x = a_red / 2.0;
  double b2y = a_red / 2.0;
  double b2z = a_red / 2.0;
  int xdim = 0;

  for (xdim = 20; xdim <= 50; xdim += 5) {
    int ydim = xdim;
    int zdim = xdim;
    // allocate memory for storage of 3D lattice points
    double *arrayX = malloc(2 * xdim * ydim * zdim * sizeof(double));
    double *arrayY = malloc(2 * xdim * ydim * zdim * sizeof(double));
    double *arrayZ = malloc(2 * xdim * ydim * zdim * sizeof(double));

    double *arrayX1 = malloc(xdim * ydim * zdim * sizeof(double));
    double *arrayY1 = malloc(xdim * ydim * zdim * sizeof(double));
    double *arrayZ1 = malloc(xdim * ydim * zdim * sizeof(double));

    double *arrayX2 = malloc(xdim * ydim * zdim * sizeof(double));
    double *arrayY2 = malloc(xdim * ydim * zdim * sizeof(double));
    double *arrayZ2 = malloc(xdim * ydim * zdim * sizeof(double));

    // call a function to create lattice points for Argon
    createLattice(xdim, ydim, zdim, b1x, b1x, b1x, arrayX1, arrayY1, arrayZ1, a_red);
    createLattice(xdim, ydim, zdim, b2x, b2x, b2x, arrayX2, arrayY2, arrayZ2, a_red);

    // // output to a text file
    // FILE *fp = fopen("output.txt", "w");
    //
    // outToFile(fp, xdim, ydim, zdim, arrayX1, arrayY1, arrayZ1, sig, 1);
    // outToFile(fp, xdim, ydim, zdim, arrayX2, arrayY2, arrayZ2, sig, 2);
    //
    // fclose(fp);

    // mergy two arrays into one
    mergeArray(xdim, ydim, zdim, arrayX1, arrayY1, arrayZ1, arrayX2, arrayY2, arrayZ2, arrayX, arrayY, arrayZ);

    // calculate LJ energy
    double LJenergy = LennardJones(xdim, ydim, zdim, arrayX, arrayY, arrayZ);

    // printf("LJ energy in Reduced unit is %.10e\n", LJenergy);
    // printf("LJ energy in SI unit is %.10e eV\n", eps * LJenergy);

    // cohesive energy energy/atom in SI units
    double cohesiveEnergy = eps * LJenergy / (2 * xdim * ydim * zdim);
    printf("%.4f\n", cohesiveEnergy);

    // free memory
    free(arrayX);
    free(arrayY);
    free(arrayZ);

    free(arrayX1);
    free(arrayY1);
    free(arrayZ1);

    free(arrayX2);
    free(arrayY2);
    free(arrayZ2);

  } // end of for loop

  return 0;
} // end of main
