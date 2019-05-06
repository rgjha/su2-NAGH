#include "read_in.h"

void read_in(Gauge_Field &u){

ifstream f_read;
int LDUM, TDUM;
 double BETADUM, DTDUM;
double  NCOLORDUM;
Umatrix dummy,dummy3;
Lattice_Vector x;
int site;



f_read.open("config");
if(f_read.bad()){
cout << "error opening config file to read\n" << flush;}

f_read >> LDUM >> TDUM >> BETADUM >>NCOLORDUM >> DTDUM;
cout << "parameters of input config\n";

cout << "electroweak coupling-1 is " << BETADUM << "\n";
cout << "number of electroweak colors-1 " << NCOLORDUM << endl;
cout << "time step used for input config " << DTDUM << "\n" << flush;

site=0;
while(loop_over_lattice(x,site)){
for(int j=0;j<D;j++){
f_read >> dummy;
u.set(x,j,dummy);}
}



f_read.close();
return;
}

