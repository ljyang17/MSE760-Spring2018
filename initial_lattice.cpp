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
  for (int i = 0; i < xdim; i++) {
    for (int j = 0; j < ydim; j++) {
      for (int k = 0; k < zdim; k++) {
        int offset = k + zdim * j + ydim * zdim * i;
        arrayX[offset] = (initialx + offset) * a_red;
        arrayY[offset] = (initialy + offset) * a_red;
        arrayZ[offset] = (initialz + offset) * a_red;
      }
    }
  } // end of for loop
} // end of function

// output lattice points to text file
void outToFile(FILE* fp, int xdim, int ydim, int zdim,
               double* arrayX, double* arrayY, double* arrayZ,
               double sig, int atomID) {

  for (int i = 0; i < xdim; i++) {
    for (int j = 0; j < ydim; j++) {
      for (int k = 0; k < zdim; k++) {
        int offset = k + zdim * j + ydim * zdim * i;

        fprintf(fp, "cell ID: %d, atom %d, x=%.6e, ", offset, atomID, sig * arrayX[offset]);
        fprintf(fp, "y=%.6e, ", sig * arrayY[offset]);
        fprintf(fp, "z=%.6e\n", sig * arrayZ[offset]);
      }
    }
  } // end of for loop
} // end of function

// merge two arrays into one
// void mergeArray(int xdim, int ydim, int zdim,
//                 double* arrayX1, double* arrayY1, double* arrayZ1,
//                 double* arrayX2, double* arrayY2, double* arrayZ2,
//                 double* arrayX3, double* arrayY3, double* arrayZ3,
//                 double* arrayX4, double* arrayY4, double* arrayZ4,
//                 double* arrayX, double* arrayY, double* arrayZ) {
//
//   int quarterNumOfAtoms = xdim * ydim * zdim;
//
//   for (int i = 0; i < quarterNumOfAtoms; i++) {
//     arrayX[i] = arrayX1[i];
//     arrayX[i + quarterNumOfAtoms] = arrayX2[i];
//     arrayX[i + 2 * quarterNumOfAtoms] = arrayX3[i];
//     arrayX[i + 3 * quarterNumOfAtoms] = arrayX4[i];
//
//     arrayY[i] = arrayY1[i];
//     arrayY[i + quarterNumOfAtoms] = arrayY2[i];
//     arrayY[i + 2 * quarterNumOfAtoms] = arrayY3[i];
//     arrayY[i + 3 * quarterNumOfAtoms] = arrayY4[i];
//
//     arrayZ[i] = arrayZ1[i];
//     arrayZ[i + quarterNumOfAtoms] = arrayZ2[i];
//     arrayZ[i + 2 * quarterNumOfAtoms] = arrayZ3[i];
//     arrayZ[i + 3 * quarterNumOfAtoms] = arrayZ4[i];
//   } // end of for loop
// } // end of function

// calculate distance between 2 points
double distanceSquare(double x1, double y1, double z1,
                      double x2, double y2, double z2){
  return pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2);
}

// calculate LJ energy
double LennardJones(int xdim, int ydim, int zdim,
                    double* arrayX1, double* arrayY1, double* arrayZ1,
                    double* arrayX2, double* arrayY2, double* arrayZ2,
                    double* arrayX3, double* arrayY3, double* arrayZ3,
                    double* arrayX4, double* arrayY4, double* arrayZ4){

  int numOfAtoms = 4 * xdim * ydim * zdim;
  double distSquare[4];
  double LJenergy = 0.0;

  for (int i = 0; i < numOfAtoms - 1; i++) {
    for (int j = i + 1; j < numOfAtoms; j++) {
      distSquare[0] = distanceSquare(arrayX1[i], arrayY1[i],arrayZ1[i], arrayX1[j], arrayY1[j],arrayZ1[j]);
      distSquare[1] = distanceSquare(arrayX2[i], arrayY2[i],arrayZ2[i], arrayX2[j], arrayY2[j],arrayZ2[j]);
      distSquare[2] = distanceSquare(arrayX3[i], arrayY3[i],arrayZ3[i], arrayX3[j], arrayY3[j],arrayZ3[j]);
      distSquare[3] = distanceSquare(arrayX4[i], arrayY4[i],arrayZ4[i], arrayX4[j], arrayY4[j],arrayZ4[j]);
      LJenergy += 4 * (pow((1 / distSquare[0]), 6) - pow((1 / distSquare[0]), 3) +
                       pow((1 / distSquare[1]), 6) - pow((1 / distSquare[1]), 3) +
                       pow((1 / distSquare[2]), 6) - pow((1 / distSquare[2]), 3) +
                       pow((1 / distSquare[3]), 6) - pow((1 / distSquare[3]), 3));
      printf("%f\n", pow((1 / distSquare[0]), 6) - pow((1 / distSquare[0]), 3));       
    }
  } // end of for loop
  return (double) LJenergy / numOfAtoms;
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

  double b2x = 1.0 / 2.0;
  double b2y = 1.0 / 2.0;
  double b2z = 0.0;

  double b3x = 1.0 / 2.0;
  double b3y = 0.0;
  double b3z = 1.0 / 2.0;

  double b4x = 0.0;
  double b4y = 1.0 / 2.0;
  double b4z = 1.0 / 2.0;

  // input argument determine size of cell
  int xdim = atoi(argv[1]);
  int ydim = atoi(argv[2]);
  int zdim = atoi(argv[3]);

  // allocate memory for storage of 3D lattice points
  double *arrayX1 = (double *)malloc(xdim * ydim * zdim * sizeof(double));
  double *arrayY1 = (double *)malloc(xdim * ydim * zdim * sizeof(double));
  double *arrayZ1 = (double *)malloc(xdim * ydim * zdim * sizeof(double));

  double *arrayX2 = (double *)malloc(xdim * ydim * zdim * sizeof(double));
  double *arrayY2 = (double *)malloc(xdim * ydim * zdim * sizeof(double));
  double *arrayZ2 = (double *)malloc(xdim * ydim * zdim * sizeof(double));

  double *arrayX3 = (double *)malloc(xdim * ydim * zdim * sizeof(double));
  double *arrayY3 = (double *)malloc(xdim * ydim * zdim * sizeof(double));
  double *arrayZ3 = (double *)malloc(xdim * ydim * zdim * sizeof(double));

  double *arrayX4 = (double *)malloc(xdim * ydim * zdim * sizeof(double));
  double *arrayY4 = (double *)malloc(xdim * ydim * zdim * sizeof(double));
  double *arrayZ4 = (double *)malloc(xdim * ydim * zdim * sizeof(double));

  // call a function to create lattice points for Argon
  createLattice(xdim, ydim, zdim, b1x, b1y, b1z, arrayX1, arrayY1, arrayZ1, a_red);
  createLattice(xdim, ydim, zdim, b2x, b2y, b2z, arrayX2, arrayY2, arrayZ2, a_red);
  createLattice(xdim, ydim, zdim, b3x, b3y, b3z, arrayX3, arrayY3, arrayZ3, a_red);
  createLattice(xdim, ydim, zdim, b4x, b4y, b4z, arrayX4, arrayY4, arrayZ4, a_red);


  // output to a text file
  remove("output.txt");
	FILE *fp = fopen("output.txt", "w");

  outToFile(fp, xdim, ydim, zdim, arrayX1, arrayY1, arrayZ1, sig, 1);
  outToFile(fp, xdim, ydim, zdim, arrayX2, arrayY2, arrayZ2, sig, 2);
  outToFile(fp, xdim, ydim, zdim, arrayX3, arrayY3, arrayZ3, sig, 3);
  outToFile(fp, xdim, ydim, zdim, arrayX4, arrayY4, arrayZ4, sig, 4);

  fclose(fp);

  // calculate LJ energy
  double LJenergy = LennardJones(xdim, ydim, zdim,
                                 arrayX1, arrayY1, arrayZ1, arrayX2, arrayY2, arrayZ2,
                                 arrayX3, arrayY3, arrayZ3, arrayX4, arrayY4, arrayZ4);
  printf("LJ energy in Reduced unit is %.10f\n", LJenergy);
  printf("LJ energy in SI unit is %.10f eV\n", eps * LJenergy);

  // free memory
  free(arrayX1);
  free(arrayY1);
  free(arrayZ1);

  free(arrayX2);
  free(arrayY2);
  free(arrayZ2);

  free(arrayX3);
  free(arrayY3);
  free(arrayZ3);

  free(arrayX4);
  free(arrayY4);
  free(arrayZ4);

  return 0;
}
