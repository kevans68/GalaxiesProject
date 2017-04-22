#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <time.h>
#include "nrutil.h"
#include "nr.h"
#define SCALE 50
#define SCALETWO .3
// #define G1 .4
// #define G2 .8
// #define OS1 1 
// #define OS2 -1

double randR(double min, double max)
{
  return min + (double) (rand() / (double) (RAND_MAX) * (max - min));
}


int main() 
{
  char xy[20];
  double u, v, w, x, y, z, xshift, yshift, zshift, vx, vy, vz;

  sprintf(xy,"xy.dat");

  FILE* xy_f = fopen(xy, "w");
  double i;

  // cacluate 500 points within the circle of radius r
  double c=0;
  int j;

  // radius of circle
  double r = 1;

  //radius of sphere 
  double rsphere = .125;
  // radius of Andromeda 
  //double r = 27; 


  // define density
  double A = M_PI*pow(r,2);
  double p_r = exp(-r);
  double * G1;
  double * shift; 
  double * sign;
  sign = dvector(1,2);
  G1 = dvector(1,2); 
  //G1[1] = .4;
  G1[1] =0;
  G1[2] = 0;
  shift = dvector(1,2); 
  shift[1] = -20; 
  shift[2] = 20;
  sign[1]= -1; 
  sign[2] = 1; 



  for(j=1;j<=2;j++)
  {
    while(c<=400)
    {
      // create exponential distribution of random points x and y 
      u = randR(-1*SCALE,1*SCALE);
      v = randR(-1*SCALE,1*SCALE);
      w = randR(-1*SCALE,1*SCALE);
      x = sqrt(-2*log(u))-1;
      y = sqrt(-2*log(v))-1;
      z = sqrt(-2*log(w))-1;

      //calculate the area
      //double bound = (M_PI*pow(r,2))*SCALE;
      double bound = sqrt( pow(x,2) + pow(y,2));

      // if the point is within the area of the circle the get other initial conditions 
      // and add to file initial conditions file. 
      if(bound>rsphere && bound<=r)// && c<=150) 
      {
        //double rcalc = pow(x,2) +pow(y,2);
        //double area = rcalc*M_PI;

        x = x*27;
        y=y*27;
        double r = sqrt(pow(x,2) +pow(y,2));

        xshift = (x+shift[j]);
        yshift = (y+shift[j]);
        z = (G1[j]*x);

        v = sign[j]*sqrtf(p_r*A/r);
        vy =(v*yshift)/r;
        vx = (v*yshift)/r;

        if(G1[j] != 0)
        {
          vz = 5/2 *vx ;
        }
        else 
        {
          double vz = 0;
        }

        c = c+1; 
        printf("%lf\n",c);


        fprintf(xy_f,"%lf %lf %lf %lf %lf %lf \n",xshift,yshift,z,vx,vy,vz);
      }

      while(c>400 && c<= 499)
      {

        x = randR(-1*SCALETWO,1*SCALETWO);
        y = randR(-1*SCALETWO,1*SCALETWO);
        z = randR(-1*SCALETWO,1*SCALETWO);

        double boundsphere = sqrt(pow(x,2) + pow(y,2) + pow(z,2));

        if(boundsphere<=rsphere)
        {

          x = x*27;
          y=y*27;
          z=z*27; 
          double r = sqrt(pow(x,2) + pow(y,2) + pow(z,2));

          xshift = (x+shift[j]);
          yshift = (y+shift[j]);

          v = sign[j]*sqrtf(p_r*A/r);
          vy =(v*yshift)/r;
          vx = (v*yshift)/r;

          if(G1[j] != 0)
          {
            vz = 5/2 *vx ;
          }
          else 
          {
            double vz = 0;
          }
          c =c+1;

          printf("%lf\n",c);
          fprintf(xy_f,"%lf %lf %lf %lf %lf %lf \n",xshift,yshift,z,vx,vy,vz);


        }


      }

      //double vy = 3.7*sign[j];
      //printf("%lf\n", rcalc); 
      //fprintf(xy_f,"%lf %lf %lf %lf %lf %lf \n",x,vx,y,vy,z,vz);
      //fprintf(xy_f,"%lf %lf %lf %lf %lf %lf \n",xshift,vx,yshift,vy,z,vz);
    }
    c = 0;

  }
  return 0;

}
