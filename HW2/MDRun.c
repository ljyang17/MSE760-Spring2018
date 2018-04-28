#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#define OUTPOS false
#define OUTFORCE true
#define OUTENERGY true
#define OUTGR true

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#define SIGMA 3.4e-10
#define EPSILON 1.04e-2
#define MASS 6.63e-26
#define KB 8.617e-5
#define EPSILONS 1.6663e-21

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
        arrayX[offset] = (initialx + i) * a_red;
        arrayY[offset] = (initialy + j) * a_red;
        arrayZ[offset] = (initialz + k) * a_red;
      }
    }
  } // end of for loop
} // end of function

// output lattice points to text file
void outToFilePos(FILE* fp, int xdim, int ydim, int zdim,
               double* arrayX, double* arrayY, double* arrayZ, int atomID) {

  for (int i = 0; i < xdim; i++) {
    for (int j = 0; j < ydim; j++) {
      for (int k = 0; k < zdim; k++) {
        int offset = k + zdim * j + ydim * zdim * i;

        fprintf(fp, "%d %d %.6e ", offset, atomID, SIGMA * arrayX[offset]);
        fprintf(fp, "%.6e ", SIGMA * arrayY[offset]);
        fprintf(fp, "%.6e\n", SIGMA * arrayZ[offset]);
        // fprintf(fp, "%d %d %.6e ", offset, atomID, SIGMA * arrayX[offset]);
        // fprintf(fp, "%.6e ", SIGMA * arrayY[offset]);
        // fprintf(fp, "%.6e\n", SIGMA * arrayZ[offset]);
      }
    }
  } // end of for loop
} // end of function

// output force to text file
void outToFileForce(FILE* fp, int xdim, int ydim, int zdim,
                    double* arrayX, double* arrayY, double* arrayZ) {
    int numOfAtoms = 4 * xdim * ydim * zdim;
    int cutOff = 50;
    for (int i = 0; i < numOfAtoms; i++) {
        // fprintf(fp, "The force between atom %d and atom %d is: Fx=%.6e, ", i, j, EPSILON / SIGMA * arrayX[offset]);
        // fprintf(fp, "Fy=%.6e, ", EPSILON / SIGMA * arrayY[offset]);
        // fprintf(fp, "Fz=%.6e\n", EPSILON / SIGMA * arrayZ[offset]);
        fprintf(fp, "%d %.6e ", i, EPSILONS / SIGMA * arrayX[i]);
        fprintf(fp, "%.6e ", EPSILONS / SIGMA * arrayY[i]);
        fprintf(fp, "%.6e\n", EPSILONS / SIGMA * arrayZ[i]);
        // fprintf(fp, "%d %d %.6e ", offset, atomID, SIGMA * arrayX[offset]);
        // fprintf(fp, "%.6e ", SIGMA * arrayY[offset]);
        // fprintf(fp, "%.6e\n", SIGMA * arrayZ[offset]);
  } // end of for loop
} // end of function
// output force to text file
void outToFileEnergy(FILE* fp, int stepTotal, int xdim, int ydim, int zdim,
                    double* arrayX, double* arrayY, double* arrayZ) {
    double numOfAtoms;
    numOfAtoms = (double) 4 * xdim * ydim * zdim;
    for (int step = 0; step < stepTotal; step++) {

        // fprintf(fp, "The force between atom %d and atom %d is: Fx=%.6e, ", i, j, EPSILON / SIGMA * arrayX[offset]);
        // fprintf(fp, "Fy=%.6e, ", EPSILON / SIGMA * arrayY[offset]);
        // fprintf(fp, "Fz=%.6e\n", EPSILON / SIGMA * arrayZ[offset]);
        fprintf(fp, "%d %.6e ", step, EPSILON * arrayX[step] / numOfAtoms);
        fprintf(fp, "%.6e ", EPSILON * arrayY[step] / numOfAtoms);
        fprintf(fp, "%.6e\n", EPSILON * arrayZ[step] / numOfAtoms);
        // fprintf(fp, "%d %d %.6e ", offset, atomID, SIGMA * arrayX[offset]);
        // fprintf(fp, "%.6e ", SIGMA * arrayY[offset]);
        // fprintf(fp, "%.6e\n", SIGMA * arrayZ[offset]);
  } // end of for loop
} // end of function

// output force to text file
void outToFileGr(FILE* fp, double numOfBins, double* gr) {
    for (int i = 0; i < numOfBins; i++) {
        fprintf(fp, "%d %.6e\n", i, gr[i]);
  } // end of for loop
} // end of function

// merge two arrays into one
void mergeArray(int xdim, int ydim, int zdim,
                double* arrayX1, double* arrayY1, double* arrayZ1,
                double* arrayX2, double* arrayY2, double* arrayZ2,
                double* arrayX3, double* arrayY3, double* arrayZ3,
                double* arrayX4, double* arrayY4, double* arrayZ4,
                double* arrayX, double* arrayY, double* arrayZ) {

  int quarterNumOfAtoms = xdim * ydim * zdim;

  for (int i = 0; i < quarterNumOfAtoms; i++) {
    arrayX[i] = arrayX1[i];
    arrayX[i + quarterNumOfAtoms] = arrayX2[i];
    arrayX[i + 2 * quarterNumOfAtoms] = arrayX3[i];
    arrayX[i + 3 * quarterNumOfAtoms] = arrayX4[i];

    arrayY[i] = arrayY1[i];
    arrayY[i + quarterNumOfAtoms] = arrayY2[i];
    arrayY[i + 2 * quarterNumOfAtoms] = arrayY3[i];
    arrayY[i + 3 * quarterNumOfAtoms] = arrayY4[i];

    arrayZ[i] = arrayZ1[i];
    arrayZ[i + quarterNumOfAtoms] = arrayZ2[i];
    arrayZ[i + 2 * quarterNumOfAtoms] = arrayZ3[i];
    arrayZ[i + 3 * quarterNumOfAtoms] = arrayZ4[i];
  } // end of for loop
} // end of function

// calculate distance between 2 points
double distanceSquare(double x1, double y1, double z1){
  return pow(x1, 2) + pow(y1, 2) + pow(z1, 2);
}

// not applying PBC
void noPBC(int xdim, int ydim, int zdim, int currentDim, double a_red,
            double* posArray, double* posArrayij){
    int numOfAtoms = 4 * xdim * ydim * zdim;
    for (int i = 0; i < numOfAtoms; i++) {
      for (int j = 0; j < numOfAtoms; j++) {
              posArrayij[j + i * numOfAtoms] = posArray[i] - posArray[j];
      }
    } // end of for loop
}

// applying PBC
void addPBC(int xdim, int ydim, int zdim, int currentDim, double a_red,
            double* posArray, double* posArrayij){
    int numOfAtoms = 4 * xdim * ydim * zdim;
    // upper triangular matrix
    for (int i = 0; i < numOfAtoms ; i++) {
      for (int j = 0; j < numOfAtoms; j++) {
        posArrayij[j + i * numOfAtoms] = posArray[i] - posArray[j];
        if (posArrayij[j + i * numOfAtoms] < - (currentDim * a_red) / 2.0)
          posArrayij[j + i * numOfAtoms] = posArrayij[j + i * numOfAtoms] + (currentDim * a_red);

        else if (posArrayij[j + i * numOfAtoms] > (currentDim * a_red) / 2.0)
          posArrayij[j + i * numOfAtoms] = posArrayij[j + i * numOfAtoms] - (currentDim * a_red);
        // else
        //   posArrayij[j + i * numOfAtoms] = posArrayij[j + i * numOfAtoms];
        }
    } // end of for loop
    // // diagonal and lower triangular matrix
    // for (int i = 0; i < numOfAtoms; i++) {
    //   for (int j = 0; j <= i; j++) {
    //         posArrayij[j + i * numOfAtoms] = posArray[j] - posArray[i];
    //         if (posArrayij[j + i * numOfAtoms] < - (currentDim * a_red) / 2.0)
    //           posArrayij[j + i * numOfAtoms] = posArrayij[j + i * numOfAtoms] + (currentDim * a_red);
    //
    //         else if (posArrayij[j + i * numOfAtoms] > (currentDim * a_red) / 2.0)
    //           posArrayij[j + i * numOfAtoms] = posArrayij[j + i * numOfAtoms] - (currentDim * a_red);
    //         else
    //           posArrayij[j + i * numOfAtoms] = posArrayij[j + i * numOfAtoms];
    //         }
    // } // end of for loop
}

// generate 3 random numbers
void randomVelGenerate(int xdim, int ydim, int zdim,
                       double* arrayVx, double* arrayVy, double* arrayVz){
  srand(time(NULL));   // should only be called once
  int numOfAtoms = 4 * xdim * ydim * zdim;

  for (int i = 0; i < numOfAtoms; i++) {
    double zeta1, zeta2;
    double zetaSquare = 1.0;
    do {
      zeta1 = (double) rand() / RAND_MAX;
      zeta2 = (double) rand() / RAND_MAX;
      // printf("%.10f %.10f\n", zeta1, zeta2);
      zeta1 = 1 - 2 * zeta1;
      zeta2 = 1 - 2 * zeta2;
      zetaSquare = zeta1 * zeta1 + zeta2 * zeta2;
    } while(zetaSquare >= 1.0);

      arrayVx[i] = 2 * zeta1 * sqrt(1 - zetaSquare);
      arrayVy[i] = 2 * zeta2 * sqrt(1 - zetaSquare);
      arrayVz[i] = 1 - 2 * zetaSquare;
    }
    // for (int i = 0; i < numOfAtoms; i++) {
    //     printf("%f %f %f\n", arrayVx[i], arrayVy[i], arrayVz[i]);
    // } // end of for loop

}

// calculate average temperature
void calcAveTemp(int xdim, int ydim, int zdim,
                  double givenTemp_red,
                  double* arrayVx, double* arrayVy, double* arrayVz){
  int numOfAtoms = 4 * xdim * ydim * zdim;
  double v0 = sqrt(3 * givenTemp_red);
  double sumTemp = 0.0;

  for (int i = 0; i < numOfAtoms; i++) {
    double Vx = v0 * arrayVx[i];
    double Vy = v0 * arrayVy[i];
    double Vz = v0 * arrayVz[i];
    double VSquare = Vx * Vx + Vy * Vy + Vz * Vz;
    sumTemp += VSquare;
    arrayVx[i] = Vx;
    arrayVy[i] = Vy;
    arrayVz[i] = Vz;
    }
    sumTemp /= 3 * numOfAtoms;
    if (fabs(sumTemp - givenTemp_red) > 1e-3){
      printf("Error! Caculated average temperature is not the same as given temperature\n");
      printf("Caculated average temperature is %f, given temperature is %f\n", sumTemp, givenTemp_red);
    }
    printf("Caculated average temperature in reduced unit is %f\n", sumTemp);
    printf("Caculated average temperature in SI unit is %f (K)\n", sumTemp * EPSILON / KB);
}
// calculate average temperature
void calcAveTempStep(FILE* fp, int xdim, int ydim, int zdim,
                  double givenTemp_red, int step, double* tempArray,
                  double* arrayVx, double* arrayVy, double* arrayVz){
  int numOfAtoms = 4 * xdim * ydim * zdim;
  double sumTemp = 0.0;
    for (int i = 0; i < numOfAtoms; i++) {
      double Vx = arrayVx[i];
      double Vy = arrayVy[i];
      double Vz = arrayVz[i];
      double VSquare = Vx * Vx + Vy * Vy + Vz * Vz;
      sumTemp += VSquare;
      }
      sumTemp /= 3 * numOfAtoms;
      tempArray[step] = sumTemp * EPSILON / KB;
      fprintf(fp, "%d %.6e\n", step, tempArray[step]);
      // printf("Caculated average temperature in reduced unit is %f\n", sumTemp);
      // printf("Caculated average temperature in SI unit is %f (K)\n", sumTemp * EPSILON / KB);
}

// calculate LJ energy
void LJEnergyForce(int xdim, int ydim, int zdim, double a_red, int PBC, int step,
                   double* arrayX, double* arrayY, double* arrayZ,
                   double* arrayXij, double* arrayYij, double* arrayZij,
                   double* forceXij, double* forceYij, double* forceZij, double* LJenergy){

  int numOfAtoms = 4 * xdim * ydim * zdim;
  double distSquare = 0.0;
  double LJ = 0.0;

  if (PBC)
  {
  addPBC(xdim, ydim, zdim, xdim, a_red, arrayX, arrayXij);
  addPBC(xdim, ydim, zdim, ydim, a_red, arrayY, arrayYij);
  addPBC(xdim, ydim, zdim, zdim, a_red, arrayZ, arrayZij);
  }

  else
  {
  noPBC(xdim, ydim, zdim, xdim, a_red, arrayX, arrayXij);
  noPBC(xdim, ydim, zdim, ydim, a_red, arrayY, arrayYij);
  noPBC(xdim, ydim, zdim, zdim, a_red, arrayZ, arrayZij);
  }
  // for (int i = 0; i < numOfAtoms; i++) {
  //   for (int j = 0; j < numOfAtoms; j++) {
  //     int offset = j + i * numOfAtoms;
  //     printf("%f ", arrayXij[offset]);
  //   }
  //   printf("\n");
  // } // end of for loop

  for (int i = 0; i < numOfAtoms; i++) {
      double forceX = 0.0;
      double forceY = 0.0;
      double forceZ = 0.0;
      for (int j = 0; j < numOfAtoms; j++) {
          double Xij = arrayXij[j + i * numOfAtoms];
          double Yij = arrayYij[j + i * numOfAtoms];
          double Zij = arrayZij[j + i * numOfAtoms];
          distSquare = distanceSquare(Xij, Yij, Zij);

          if (distSquare != 0.0){
            forceX += 48.0 * Xij * (pow((1.0 / distSquare), 7) - 1.0 / 2.0 * pow((1.0 / distSquare), 4));
            forceY += 48.0 * Yij * (pow((1.0 / distSquare), 7) - 1.0 / 2.0 * pow((1.0 / distSquare), 4));
            forceZ += 48.0 * Zij * (pow((1.0 / distSquare), 7) - 1.0 / 2.0 * pow((1.0 / distSquare), 4));
          }

      // if (distSquare<0.5){
      //   printf("%d %d %.6f\n", i, j, distSquare);
      //     printf("%.10f\n", 4.0 * (pow((1.0 / distSquare), 6) - pow((1.0 / distSquare), 3)));}
      // printf("%.10f\n", LJenergy);

      // printf("%.10f\n", 4.0 * (pow((1.0 / distSquare), 6) - pow((1.0 / distSquare), 3)));
      // printf("%.10f\n", pow((1.0 / distSquare), 6) - pow((1.0 / distSquare), 3));
      }
      forceXij[i] = forceX;
      forceYij[i] = forceY;
      forceZij[i] = forceZ;
  } // end of for loop

  for (int i = 0; i < numOfAtoms - 1; i++) {
    for (int j = i + 1; j < numOfAtoms; j++) {
      double Xij = arrayXij[j + i * numOfAtoms];
      double Yij = arrayYij[j + i * numOfAtoms];
      double Zij = arrayZij[j + i * numOfAtoms];
      distSquare = distanceSquare(Xij, Yij, Zij);

      LJ += 4.0 * (pow((1.0 / distSquare), 6) - pow((1.0 / distSquare), 3));
      // if (distSquare<0.5){
      //   printf("%d %d %.6f\n", i, j, distSquare);
      //     printf("%.10f\n", 4.0 * (pow((1.0 / distSquare), 6) - pow((1.0 / distSquare), 3)));}

      // printf("%.10f\n", 4.0 * (pow((1.0 / distSquare), 6) - pow((1.0 / distSquare), 3)));
      // printf("%.10f\n", pow((1.0 / distSquare), 6) - pow((1.0 / distSquare), 3));
    }
  } // end of for loop
  LJenergy[step] = LJ;
} // end of function

// velocity Verlet algorithm
void velVerletIntegration(int xdim, int ydim, int zdim, double dT, double a_red, int PBC, int step,
                          double* LJenergy, double* kinEnergy, double* totEnergy,
                          double* arrayXij, double* arrayYij, double* arrayZij,
                          double* forceXij, double* forceYij, double* forceZij,
                          double* posArrayMX, double* posArrayMY, double* posArrayMZ,
                          double* velArrayMX, double* velArrayMY, double* velArrayMZ){
  int numOfAtoms = 4 * xdim * ydim * zdim;
  double kinetic = 0.0;

  for (int i = 0; i < numOfAtoms; i++) {
     // need to fix index issue--fixed
     // update half step velocity in place
    //  printf("vx before %f \n", velArrayMX[i]);
    //  printf("second before %f \n", forceXij[i] * dT);
    //  printf("second before %f \n", dT);
    //  printf("second before %f \n", forceXij[i]);
     velArrayMX[i] = velArrayMX[i] + 1.0 / 2.0 * forceXij[i] * dT;
    //  printf("vx after %f \n", velArrayMX[i]);
     velArrayMY[i] = velArrayMY[i] + 1.0 / 2.0 * forceYij[i] * dT;
     velArrayMZ[i] = velArrayMZ[i] + 1.0 / 2.0 * forceZij[i] * dT;
     // update position in place
     posArrayMX[i] = posArrayMX[i] + velArrayMX[i] * dT;
     posArrayMY[i] = posArrayMY[i] + velArrayMY[i] * dT;
     posArrayMZ[i] = posArrayMZ[i] + velArrayMZ[i] * dT;
   }
   // check BC for updated position, update arrayXij, arrayYij, arrayZij
   if (PBC){
     addPBC(xdim, ydim, zdim, xdim, a_red, posArrayMX, arrayXij);
     addPBC(xdim, ydim, zdim, ydim, a_red, posArrayMY, arrayYij);
     addPBC(xdim, ydim, zdim, zdim, a_red, posArrayMZ, arrayZij);
   }
   // calculate force based on position
   LJEnergyForce(xdim, ydim, zdim, a_red, PBC, step,
                 posArrayMX, posArrayMY, posArrayMZ,
                 arrayXij, arrayYij, arrayZij,
                 forceXij, forceYij, forceZij, LJenergy);
                //  for (int i = 0; i < numOfAtoms; i++) {
                //      printf("half %f ", forceXij[i]);
                //     printf("\n");
                //  } // end of for loop
   // update full step velocity in place
   for (int i = 0; i < numOfAtoms; i++) {
      velArrayMX[i] = velArrayMX[i] + 1.0 / 2.0 * forceXij[i] * dT;
      velArrayMY[i] = velArrayMY[i] + 1.0 / 2.0 * forceYij[i] * dT;
      velArrayMZ[i] = velArrayMZ[i] + 1.0 / 2.0 * forceZij[i] * dT;
      // printf("vx after %f \n", velArrayMX[i]);
      kinetic += 1 / 2.0 * (velArrayMX[i] * velArrayMX[i] + velArrayMY[i] * velArrayMY[i] + velArrayMZ[i] * velArrayMZ[i]);
      // printf("kinetic %f ", kinetic);
    }
    kinEnergy[step] =  kinetic;
    // printf("kinetic %f ", kinEnergy[step]);
    totEnergy[step] = kinEnergy[step] + LJenergy[step];
 }
 // calculate radial distribution function
 void calcRadDis(int xdim, int ydim, int zdim, double a_red, int PBC, double cutoff, double dr,
                double* posArrayMX, double* posArrayMY, double* posArrayMZ,
                 double* arrayXij, double* arrayYij, double* arrayZij, int* H){
   // check BC for updated position, update arrayXij, arrayYij, arrayZij
   if (PBC){
     addPBC(xdim, ydim, zdim, xdim, a_red, posArrayMX, arrayXij);
     addPBC(xdim, ydim, zdim, ydim, a_red, posArrayMY, arrayYij);
     addPBC(xdim, ydim, zdim, zdim, a_red, posArrayMZ, arrayZij);
   }
   int bin;
   int numOfAtoms = 4 * xdim * ydim * zdim;

   for (int i = 0; i < numOfAtoms - 1; i++) {
     for (int j = i + 1; j < numOfAtoms; j++) {
       double Xij = arrayXij[j + i * numOfAtoms];
       double Yij = arrayYij[j + i * numOfAtoms];
       double Zij = arrayZij[j + i * numOfAtoms];
       double distSquare = distanceSquare(Xij, Yij, Zij);
       if (distSquare < cutoff * cutoff) {
         bin = (int)(sqrt(distSquare) / dr);
         H[bin] += 2;
       }
   }
 }
 }

// start main
int main(int argc, char const *argv[]) {

  //=========================parameter definition================================
  // define system size & time step from input arguments
  int xdim = atoi(argv[1]);
  int ydim = atoi(argv[2]);
  int zdim = atoi(argv[3]);
  double dT = atof(argv[4]);
  int PBC = 1;

  // define actual constants
  double a = 5.7e-10;
  double m = 6.63e-26;

  // system termperature
  double givenTemp = 500.0;

  // time step
  // double dT = 0.001;
  double giventimeTotal = 2.2e-13;
  double dr = 0.1;

  // define reduced-unit constants
  double a_red = a / SIGMA;
  double m_red = m / MASS;
  double givenTemp_red = KB * givenTemp / EPSILON;
  // double real_time = sqrt((1.0 * SIGMA * SIGMA) / EPSILON) * dT;
  double giventimeTotal_red = giventimeTotal * sqrt(EPSILONS / (m * SIGMA * SIGMA));
  // printf("%f\n", giventimeTotal_red);
  int stepTotal = (int) ceil(giventimeTotal_red / dT);
  // printf("Total simulation time steps: %d\n", stepTotal);
  // printf("%f\n", giventimeTotal_red * 1000 / (dT*1000));
  // int stepTotal = 50;

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

  //=========================system initialization===============================
  // input argument determine size of cell
  int numberOfCells = xdim * ydim * zdim;
  int numOfAtoms = 4 * numberOfCells;
  int sumOfAtoms = (4 * numberOfCells) * (4 * numberOfCells);
  double cutoff = zdim * a_red / 2.0;
  double totalVolume = numberOfCells * pow(a_red, 3);
  double totalMass = numOfAtoms * m_red;
  double rho = totalMass / totalVolume;
  int numOfBins = (int)(cutoff / dr) + 1;
  // int sumOfAtomsForce = (4 * numberOfCells) * (4 * numberOfCells);

  // allocate memory for storage of position
  double *posArrayMX = (double *)malloc(4 * numberOfCells* sizeof(double));
  double *posArrayMY = (double *)malloc(4 * numberOfCells* sizeof(double));
  double *posArrayMZ = (double *)malloc(4 * numberOfCells* sizeof(double));

  // allocate memory for distance of atoms
  double *arrayXij = (double *)calloc(sumOfAtoms , sizeof(double));
  double *arrayYij = (double *)calloc(sumOfAtoms , sizeof(double));
  double *arrayZij = (double *)calloc(sumOfAtoms , sizeof(double));

  // allocate memory for storage of 3D lattice points
  double *arrayX1 = (double *)malloc(numberOfCells* sizeof(double));
  double *arrayY1 = (double *)malloc(numberOfCells* sizeof(double));
  double *arrayZ1 = (double *)malloc(numberOfCells* sizeof(double));

  double *arrayX2 = (double *)malloc(numberOfCells* sizeof(double));
  double *arrayY2 = (double *)malloc(numberOfCells* sizeof(double));
  double *arrayZ2 = (double *)malloc(numberOfCells* sizeof(double));

  double *arrayX3 = (double *)malloc(numberOfCells* sizeof(double));
  double *arrayY3 = (double *)malloc(numberOfCells* sizeof(double));
  double *arrayZ3 = (double *)malloc(numberOfCells* sizeof(double));

  double *arrayX4 = (double *)malloc(numberOfCells* sizeof(double));
  double *arrayY4 = (double *)malloc(numberOfCells* sizeof(double));
  double *arrayZ4 = (double *)malloc(numberOfCells* sizeof(double));

  // allocate memory for force arrays
  double *forceXij = (double *)calloc(numOfAtoms , sizeof(double));
  double *forceYij = (double *)calloc(numOfAtoms , sizeof(double));
  double *forceZij = (double *)calloc(numOfAtoms , sizeof(double));

  // allocate memory for velocity arrays
  double *velArrayMX = (double *)calloc(numOfAtoms , sizeof(double));
  double *velArrayMY = (double *)calloc(numOfAtoms , sizeof(double));
  double *velArrayMZ = (double *)calloc(numOfAtoms , sizeof(double));

  // allocate memory for energy arrays
  double *LJEnergy = (double *)calloc(stepTotal , sizeof(double));
  double *kinEnergy = (double *)calloc(stepTotal , sizeof(double));
  double *totEnergy = (double *)calloc(stepTotal , sizeof(double));
  double *tempArray = (double *)calloc(stepTotal , sizeof(double));

  // allocate memory for bins
  int* H = (int*)calloc(numOfBins,sizeof(int));
  double* gr = (double*)calloc(numOfBins,sizeof(double));

  // call a function to create lattice points for Argon
  createLattice(xdim, ydim, zdim, b1x, b1y, b1z, arrayX1, arrayY1, arrayZ1, a_red);
  createLattice(xdim, ydim, zdim, b2x, b2y, b2z, arrayX2, arrayY2, arrayZ2, a_red);
  createLattice(xdim, ydim, zdim, b3x, b3y, b3z, arrayX3, arrayY3, arrayZ3, a_red);
  createLattice(xdim, ydim, zdim, b4x, b4y, b4z, arrayX4, arrayY4, arrayZ4, a_red);

  // output to a text file
  if (OUTPOS){
    remove("outputpos.txt");
  	FILE *fp = fopen("outputpos.txt", "w");

    outToFilePos(fp, xdim, ydim, zdim, arrayX1, arrayY1, arrayZ1, 1);
    outToFilePos(fp, xdim, ydim, zdim, arrayX2, arrayY2, arrayZ2, 2);
    outToFilePos(fp, xdim, ydim, zdim, arrayX3, arrayY3, arrayZ3, 3);
    outToFilePos(fp, xdim, ydim, zdim, arrayX4, arrayY4, arrayZ4, 4);
    fclose(fp);
  }
  // initialize atom position
  mergeArray(xdim, ydim, zdim, arrayX1, arrayY1, arrayZ1, arrayX2, arrayY2, arrayZ2,
             arrayX3, arrayY3, arrayZ3, arrayX4, arrayY4, arrayZ4, posArrayMX, posArrayMY, posArrayMZ);

  //=========================start working=======================================
  // static problems
  // calculate LJ energy
  LJEnergyForce(xdim, ydim, zdim, a_red, PBC, 0,
                posArrayMX, posArrayMY, posArrayMZ,
                arrayXij, arrayYij, arrayZij,
                forceXij, forceYij, forceZij, LJEnergy);

  if (OUTFORCE){
   remove("outputforce.txt");
   FILE *fp = fopen("outputforce.txt", "w");

   outToFileForce(fp, xdim, ydim, zdim, forceXij, forceYij, forceZij);
   fclose(fp);
  }
  calcRadDis(xdim, ydim, zdim, a_red, PBC, cutoff, dr,
             posArrayMX, posArrayMY, posArrayMZ,
             arrayXij, arrayYij, arrayZij, H);

  for (int i = 0; i < numOfBins; i++) {
    // double r = dr * (i + 0.5);
    double binVolume = ((i + 1) * (i + 1) * (i + 1) - i * i * i) * dr * dr * dr;
    double norm = (4.0 / 3.0) * M_PI * binVolume * rho;
    gr[i] = H[i] / norm;
    // printf("%f\n", gr[i]);
    // printf("norm %f\n", norm);
  }
  if (OUTGR){
   remove("outputgr.txt");
   FILE *fp = fopen("outputgr.txt", "w");

   outToFileGr(fp, numOfBins, gr);
   fclose(fp);
  }

  // printf("LJ energy in Reduced unit is %.10f\n", LJenergy[0]);
  // printf("LJ energy in SI unit is %.10f eV\n", EPSILON * LJenergy[0]);

  randomVelGenerate(xdim, ydim, zdim, velArrayMX, velArrayMY, velArrayMZ);
  calcAveTemp(xdim, ydim, zdim, givenTemp_red, velArrayMX, velArrayMY, velArrayMZ);
  for (int i = 0; i < numOfAtoms; i++) {
     kinEnergy[0] += 1.0 / 2.0 * (velArrayMX[i] * velArrayMX[i] + velArrayMY[i] * velArrayMY[i] + velArrayMZ[i] * velArrayMZ[i]);
   }
  totEnergy[0] = LJEnergy[0] + kinEnergy[0];
  // dynamic problems
  FILE *fp = fopen("temp", "w");
  for (int step = 1; step < stepTotal; step++) {
    velVerletIntegration(xdim, ydim, zdim, dT, a_red, PBC, step,
                         LJEnergy, kinEnergy, totEnergy,
                         arrayXij, arrayYij, arrayZij,
                         forceXij, forceYij, forceZij,
                         posArrayMX, posArrayMY, posArrayMZ,
                         velArrayMX, velArrayMY, velArrayMZ);
    printf("step: %d total energy %.10e (eV)\n", step, totEnergy[step] * EPSILON);
    calcAveTempStep(fp, xdim, ydim, zdim, givenTemp_red, step, tempArray,
                    velArrayMX, velArrayMY, velArrayMZ);
  }
  fclose(fp);
    if (OUTENERGY){
     remove("outputenergy.txt");
     FILE *fp = fopen("outputenergy.txt", "w");

     outToFileEnergy(fp, stepTotal, xdim, ydim, zdim, kinEnergy, LJEnergy, totEnergy);
     fclose(fp);
    }
    // calcRadDis(xdim, ydim, zdim, a_red, PBC, cutoff, dr,
    //            posArrayMX, posArrayMY, posArrayMZ,
    //            arrayXij, arrayYij, arrayZij, H);
    //
    // for (int i = 0; i < numOfBins; i++) {
    //   // double r = dr * (i + 0.5);
    //   double binVolume = ((i + 1) * (i + 1) * (i + 1) - i * i * i) * dr * dr * dr;
    //   double norm = (4.0 / 3.0) * M_PI * binVolume * rho;
    //   gr[i] = H[i] / norm;
    //   // printf("%f\n", gr[i]);
    //   // printf("norm %f\n", norm);
    // }
    // if (OUTGR){
    //  remove("outputgr.txt");
    //  FILE *fp = fopen("outputgr.txt", "w");
    //
    //  outToFileGr(fp, numOfBins, gr);
    //  fclose(fp);
    // }
    // // printf("%f %f\n", totalMass, totalVolume);
    // // givenTemp_red /= (250.0 / 250.0);
    // double sumTemp = 0.0;
    // for (int i = 0; i < numOfAtoms; i++) {
    //   double Vx = velArrayMX[i];
    //   double Vy = velArrayMY[i];
    //   double Vz = velArrayMZ[i];
    //   double VSquare = Vx * Vx + Vy * Vy + Vz * Vz;
    //   sumTemp += VSquare;
    //   }
    //   sumTemp /= 3 * numOfAtoms;
    //   if (fabs(sumTemp - givenTemp_red) > 1e-3){
    //     printf("Error! Caculated average temperature is not the same as given temperature\n");
    //     printf("Caculated average temperature is %f, given temperature is %f\n", sumTemp, givenTemp_red);
    //   }
    //   printf("Caculated average temperature in reduced unit is %f\n", sumTemp);
    //   printf("Caculated average temperature in SI unit is %f (K)\n", sumTemp * EPSILON / KB);


  //===========================wrap up===========================================
  // free memory
  free(posArrayMX);
  free(posArrayMY);
  free(posArrayMZ);

  free(arrayXij);
  free(arrayYij);
  free(arrayZij);

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

  free(forceXij);
  free(forceYij);
  free(forceZij);

  free(velArrayMX);
  free(velArrayMY);
  free(velArrayMZ);

  free(LJEnergy);
  free(kinEnergy);
  free(totEnergy);

  free(H);
  free(gr);

  return 0;
}
