#include "ncepgrids.h"

template <class T>
class lakes : public lambert<T> {
  public :
    lakes(void);
};
template <class T>
lakes<T>::lakes(void) {
  this->a = parameters::a; // value from IPLIB, spherical earth
  
  this->nx = 2345 ;
  this->ny = 1597 + 2000;
  this->lat1   = 19.229;
  this->lon1   = 233.723448 - 360. ;
  this->orient =  -95; 
  this->reflat = 25.0; 
  this->dx = 2539.703;
  this->dy = 2539.703;
  this->h  = 1;

  this->grid = new T[this->nx*this->ny];

//Derived parameters:
  double dpr = parameters::degrees_per_radian;
  this->an = sin(this->h*this->reflat/dpr);
  this->de = this->a * cos(this->reflat/dpr)*
            pow(tan( (this->h*this->reflat + 90.)/2. /dpr), this->an) / this->an;
  this->dr = this->de / pow(tan( (this->h*this->lat1 + 90.)/2. /dpr), this->an);
  this->dlon1 = (this->lon1 - this->orient);
  this->x0 =  - this->h*sin(this->an*this->dlon1 /dpr) * this->dr / this->dx;
  this->y0 =  + cos(this->an*this->dlon1 /dpr ) * this->dr / this->dy;
  this->antr = 1./2./this->an;
  this->de2 = this->de*this->de;
}

int main(void) {
  FILE *fin;
  lakes<unsigned char> maskout;
  global_12th<unsigned char> maskin;
  ijpt loc, loc2;
  latpt ll;
  float lmax = -90;
  
  fin = fopen("seaice_gland5min.map","r");
  maskin.binin(fin);
  fclose(fin);
  printf("input  %f %3d %3d\n",lmax, maskin.gridmax(), maskin.gridmin());

  for (loc.j = 0; loc.j < maskout.ypoints(); loc.j++) {
  for (loc.i = 0; loc.i < maskout.xpoints(); loc.i++) {
    ll = maskout.locate(loc);
    loc2 = maskin.locate(ll);
    maskout[loc] = maskin[loc2];
    if (ll.lat > lmax) lmax = ll.lat;
    if (maskout[loc] == 157) maskout[loc] = 1;
    if (maskout[loc] == 195) maskout[loc] = 2;
  }
  }
  printf("max lat = %f %3d %3d %3d %3d\n",lmax, maskout.gridmax(), maskout.gridmin(), maskout.average(), maskout.rms() );

  palette<unsigned char> alpha(3);
  alpha.set_color(0, 255,255,255);
  alpha.set_color(1, 0,0,0);
  alpha.set_color(2, 0,0,255);

  maskout.xpm("lmask.xpm",1,alpha);


  return 0;
}
