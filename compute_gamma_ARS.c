
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"
#include "tags.h"

void compute_gamma(double abh2, double *pgamma, double entropy_init, double gamma_init, double dens)
{
double mu, muh2, muin, muh2in, h2frac;
int j=0;
double e, x, T, UnitEnergy_in_cgs, UnitTime_in_s, UnitMass_in_g, UnitVelocity_in_cm_per_s, UnitLength_in_cm;
double u, a3inv;

 if(All.ComovingIntegrationOn)
    {
      /* Factors for comoving integration of hydro */
      a3inv = 1 / (All.Time * All.Time * All.Time);
   }
  else
    a3inv = 1;

mu=1.2195;
h2frac=2.0*abh2;
 
muh2in=(0.24/4.0) + ((1.0-h2frac)*0.76) + (h2frac*.76/2.0);
muh2=pow(muh2in, -1.0);

if(muh2 >= 1.22){
   mu=muh2;
}

UnitLength_in_cm= 3.085678e21;
UnitVelocity_in_cm_per_s= 1.0e5;
UnitMass_in_g= 1.989e43;  
UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);

u=(entropy_init/(gamma_init-1.0))*pow(dens*a3inv,(gamma_init-1.0));

T=mu*PROTONMASS/BOLTZMANN * (gamma_init-1.0)*u * UnitEnergy_in_cgs/ UnitMass_in_g;

e=2.7182818;
x=6100.0/T;

*pgamma=1.0 + (1.0/mu)*pow((1.5*0.24/4.0)+(1.5*0.76*(1.0-h2frac)) + ((2.5+pow(x,2.0)*pow(e,x)/pow((pow(e,x)-1.0),2.0))*(.76/2.0)*h2frac), -1.e0);

}

