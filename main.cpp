#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <algorithm>
#include <math.h>

#define SIGMA0 3.4e-10
#define EPSILON0 0.0104



void create_x(double**, double, const int, const int,int ,int );
void compute_f(double **x,double **f, const int natom, double rcutoff,double epsilon,double sigma, double L,int flag_pbc);
double compute_potential(double **x, const int natom, double rcutoff,double epsilon,double sigma,double L,int flag_pbc);
void update_x(double **x, double **v, const int natom, double dt);
void update_v(double **v, double **f, const int natom, double mass, double dt);

int main(void){
  double  potential;
  double a, rcutoff, epsilon, sigma, mass, dt, L;
  int ntimestep;
  int flag_pbc,flag_file;
  int natom;
  FILE * fp;

  flag_pbc = 1;
  flag_file = 1;
  a = 5.26e-10/SIGMA0;
  mass = 1.0;
  rcutoff=4.0;
  dt = 0.1;
  ntimestep = 1000;
  epsilon = 1.0;
  sigma = 1.0;
  if (flag_pbc==0)
   fp= fopen("Result_wo_pbc.txt","w");
  else
   fp= fopen("Result_w_pbc.txt","w");
  fprintf(fp, "#  ncell potential-energy-per-atom\n");
  for(int s = 0; s< 1; s++){
    double **x;
    const int ncell_perdim= 5+s*1;

    if(flag_pbc==0)
      natom = pow(ncell_perdim,3)*4 + pow(ncell_perdim,2)*2*3 + (ncell_perdim)*3 + 1;
    else
      natom = pow(ncell_perdim,3)*4;

    printf("natom=%d\n",natom );
    L = a*ncell_perdim;
    x = new double * [natom];

    for (int i = 0;i<natom;i++){
      x[i] = new double[3]();

    }
    create_x(x, a, natom, ncell_perdim, flag_file, flag_pbc);
    potential = compute_potential(x, natom, rcutoff, epsilon, sigma, L, flag_pbc);
    potential *= EPSILON0;
    printf(" %d %f\n", ncell_perdim*ncell_perdim*ncell_perdim, potential);
    fprintf(fp, " %d %f\n", ncell_perdim*ncell_perdim*ncell_perdim, potential);

    for(int i=0; i<natom; i++){
      delete[] x[i];
    }
    delete[] x;
  }

  fclose(fp);

  return 0;
}

//---------------------------------------------------------
void create_x(double** x, double a, const int natom, const int ncell_perdim, int flag_file, int flag_pbc){
  int index;

  for(int i=0; i<ncell_perdim; i++){
    for(int j=0; j<ncell_perdim; j++){
      for(int k=0; k<ncell_perdim; k++){
        index = (i*ncell_perdim*ncell_perdim + j*ncell_perdim + k)*4;

        x[index  ][0] = i*a;
        x[index  ][1] = j*a;
        x[index  ][2] = k*a;

        x[index+1][0] = i*a;
        x[index+1][1] = j*a + 0.5*a;
        x[index+1][2] = k*a + 0.5*a;

        x[index+2][0] = i*a + 0.5*a;
        x[index+2][1] = j*a;
        x[index+2][2] = k*a + 0.5*a;

        x[index+3][0] = i*a + 0.5*a;
        x[index+3][1] = j*a + 0.5*a;
        x[index+3][2] = k*a;
      }
    }
  }
  if(flag_pbc==0){
    // i= ncell_perdim
    for(int j=0; j< ncell_perdim; j++){
      for(int k=0; k<ncell_perdim; k++){
        index = pow(ncell_perdim,3)*4 + (j*ncell_perdim + k)*2;
        x[index  ][0] = ncell_perdim*a;
        x[index  ][1] = j*a;
        x[index  ][2] = k*a;

        x[index+1][0] = ncell_perdim*a;
        x[index+1][1] = j*a + 0.5*a;
        x[index+1][2] = k*a + 0.5*a;
      }
    }
    // j= ncell_perdim
    for(int i=0; i< ncell_perdim; i++){
      for(int k=0; k<ncell_perdim; k++){
        index = pow(ncell_perdim,3)*4 + pow(ncell_perdim,2)*2*2 + (i*ncell_perdim+k)*2;
        x[index][0] = i*a;
        x[index][1] = ncell_perdim*a;
        x[index][2] = k*a;

        x[index+1][0] = i*a + 0.5*a;
        x[index+1][1] = ncell_perdim*a;
        x[index+1][2] = k*a + 0.5*a;
      }
    }
    // k= ncell_perdim
    for(int i=0; i< ncell_perdim; i++){
      for(int j=0; j<ncell_perdim; j++){
        index = pow(ncell_perdim,3)*4 + pow(ncell_perdim,2)*2 + (i*ncell_perdim+j)*2;
        x[index][0] = i*a;
        x[index][1] = j*a;
        x[index][2] = ncell_perdim*a;

        x[index+1][0] = i*a + 0.5*a;
        x[index+1][1] = j*a + 0.5*a;
        x[index+1][2] = ncell_perdim*a;
      }
    }
    //j = K = ncell_perdim
    for(int i=0; i < ncell_perdim; i++){
      index = pow(ncell_perdim,3)*4 + pow(ncell_perdim,2)*2*3 + i;
      x[index][0] = i*a;
      x[index][1] = ncell_perdim*a;
      x[index][2] = ncell_perdim*a;
    }
    //i = k = ncell_perdim
    for(int j=0; j < ncell_perdim; j++){
      index = pow(ncell_perdim,3)*4 + pow(ncell_perdim,2)*2*3 + ncell_perdim + j;
      x[index][0] = ncell_perdim*a;
      x[index][1] = j*a;
      x[index][2] = ncell_perdim*a;
    }
    //i = j = ncell_perdim
    for(int k=0; k < ncell_perdim+1; k++){
      index = pow(ncell_perdim,3)*4 + pow(ncell_perdim,2)*2*3 + ncell_perdim*2 + k;
      x[index][0] = ncell_perdim*a;
      x[index][1] = ncell_perdim*a;
      x[index][2] = k*a;
    }
  }

  if(flag_file==1){
    FILE *f = fopen("coordinate.txt","w");
    for(int i=0; i<natom; i++){
        fprintf(f, "%d %d %e %e %e\n", i+1, 1, x[i][0]*SIGMA0, x[i][1]*SIGMA0, x[i][2]*SIGMA0 );
    }
    fclose(f);
  }

}

//---------------------------------------------------------
void compute_f(double **x,double **f, const int natom, double rcutoff,double epsilon,double sigma,double L,int flag_pbc){
  double xtmp,ytmp,ztmp;
  double delx,dely,delz;
  double rsq,fij,r2inv,r6inv;
  double rcutsq,lj1,lj2,fpair,forcelj;
  rcutsq = rcutoff*rcutoff;
  lj1 = 48.0 * epsilon * pow(sigma,12.0);
  lj2 = 24.0 * epsilon * pow(sigma,6.0);

  //f clear
  for(int s=0; s<natom; s++ ){
    f[s][0] = 0.0; f[s][1] = 0.0; f[s][2] = 0.0;
  }

  for(int i = 0; i < natom-1; i++){
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    for(int j = i+1; j < natom; j++){
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      // //pbc
      if(flag_pbc==1){
        if (delx > L/2.0) delx = delx - L;
        else if (delx < -L/2.0) delx = delx + L;
        if (dely > L/2.0) dely = dely - L;
        else if (dely < -L/2.0) dely = dely + L;
        if (delz > L/2.0) delz = delz - L;
        else if (delz < -L/2.0) delz = delz + L;
      }

      rsq = delx*delx + dely*dely + delz*delz;
      if(rsq <= rcutsq){
        //'''Fij = 48/rs^13 - 24/r^7, return vector of LJForce'''
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        forcelj = r6inv * (lj1*r6inv - lj2);
        fpair = forcelj*r2inv;
        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        f[j][0] -= delx*fpair;
        f[j][1] -= dely*fpair;
        f[j][2] -= delz*fpair;

      }

    }
  }
}
//---------------------------------------------------------
double compute_potential(double **x, const int natom, double rcutoff,double epsilon,double sigma,double L,int flag_pbc){
  double xtmp,ytmp,ztmp;
  double delx,dely,delz;
  double rsq,r2inv,r6inv;
  double rcutsq,lj3,lj4,potential,offset,ratio;
  int index = 0;
  rcutsq = rcutoff*rcutoff;
  lj3 = 4.0 * epsilon * pow(sigma,12.0);
  lj4 = 4.0 * epsilon * pow(sigma,6.0);
  potential = 0.0;
  ratio = sigma / rcutoff;
  offset = 4.0 * epsilon * (pow(ratio,12.0) - pow(ratio,6.0));
  for(int i = 0; i < natom-1; i++){
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    for(int j = i+1; j < natom; j++){
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      // //pbc
      if(flag_pbc==1){
        if (delx > L/2.0) delx = delx - L;
        else if (delx < -L/2.0) delx = delx + L;
        if (dely > L/2.0) dely = dely - L;
        else if (dely < -L/2.0) dely = dely + L;
        if (delz > L/2.0) delz = delz - L;
        else if (delz < -L/2.0) delz = delz + L;
      }

      rsq = delx*delx + dely*dely + delz*delz;
      if(rsq <= rcutsq && rsq>0.0 ){
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        potential += r6inv*(lj3*r6inv-lj4) - offset;
        // printf("i=%d,j=%d,rsq=%f,delx=%f,dely=%f,delz=%f\n",i,j,rsq,delx,dely,delz);
        // printf("i=%d,j=%d,rsq=%f,lj3=%f,lj4=%f,potential_temp=%f,offset=%f\n",i,j,rsq,lj3,lj4,r6inv*(lj3*r6inv-lj4),offset );
      }
      if(rsq < 0.00001){
        index++;
        printf("rsq=0,i=%d,j=%d,x=%f,y=%f,z=%f\n",i,j,x[i][0],x[i][1],x[i][2] );

      }

    }
  }
  return potential/(natom-index);
}
//-----------------------------------------------------------------------------
void update_v(double **v, double **f, const int natom, double mass, double dt){
  double halfdtpm = 0.5*dt/mass;
  for(int i =0; i < natom; i++){
    v[i][0] += halfdtpm*f[i][0];
    v[i][1] += halfdtpm*f[i][1];
    v[i][2] += halfdtpm*f[i][2];
  }
}
//-----------------------------------------------------------------------------
void update_x(double **x, double **v, const int natom, double dt){

  for(int i =0; i < natom; i++){
    x[i][0] += dt*v[i][0];
    x[i][1] += dt*v[i][1];
    x[i][2] += dt*v[i][2];
  }
}
