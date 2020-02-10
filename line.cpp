#include "line.h"

// Computes Polyakov lines

void line(const Gauge_Field  &U){
int site=0,t,i,j;
Lattice_Vector x,y,e_mu,e_nu;
Umatrix prod;
static int first_time=1,count=0;
static ofstream f_line1;
Complex p;

count++;

if(first_time){
f_line1.open("lines");

if(f_line1.bad()){
cerr << "failed to open lines_s file" << endl;exit(1);}
first_time=0;
}

p=Complex();

site=0;
while(loop_over_lattice(x,site)){

e_mu=Lattice_Vector(D-1);

prod=Umatrix(1);

y=x;
for(t=1;t<=T;t++){
prod=prod*U.get(y,D-1);
y=y+e_mu;
}

p=p+(1.0/(NCOLOR*SITES))*Tr(prod);

}

f_line1  << p.real() << endl;

return;
}
