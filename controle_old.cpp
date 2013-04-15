#include <iostream>
#include <stdint.h>
using namespace std;

extern "C" {
void controle_ (int32_t &nx,int32_t &ny,float *x,float *y,float *u,float *v,
		int32_t &nxc,float *xc,float *uc,float *vc,float &xa,float &xb);}

void controle_ (int32_t &nx,int32_t &ny,float *x,float *y,float *u,float *v,
		int32_t &nxc,float *xc,float *uc,float *vc,float &xa,float &xb)
{
  //cout << "Hello World!";
  cout << "nx" << nx << "ny="<< ny <<endl; 
  //for (int i = 0; i < nx; i++)
  // cout << x[i]<< "\t" <<y[i] << "\t" <<u[i+10]<< "\t" <<v[i+10] <<endl; 
  xa=23.;
  xb=26.;
  nxc=5;
  for (int i = 0; i < nxc; i++)
   {
     xc[i]=xa+(float)i*(xb-xa)/(nxc-1);
     //cout << xc[i]  <<endl;
     uc[i]=xc[i]/2.;
     vc[i]=xc[i]/4.;
   }


}
