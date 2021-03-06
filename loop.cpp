#include "loop.h"
// Computes wilson loops

void loop(const Gauge_Field  &U){

int site=0,r,t,R,M,mu,nu;
Lattice_Vector x,y,e_mu,e_nu;
Umatrix prod;
static int first_time=1,count=0;
static ofstream f_loop;
double wilson[L][L];

Gauge_Field Udag,Wdag;
Udag=Adj(U);

count++;

if(first_time){
f_loop.open("loops");
if(f_loop.bad()){
cerr << "Failed to open loops file" << "\n";exit(1);}

first_time=0;
}

for(R=1;R<=L/2;R++)
for(M=1;M<=L/2;M++){

wilson[R][M]=0.0;

site=0;
while(loop_over_lattice(x,site)){

for(mu=0;mu<D;mu++)
for(nu=0;nu<mu;nu++){

e_mu=Lattice_Vector(mu);
e_nu=Lattice_Vector(nu);

prod=Umatrix(1);

y=x;

for(r=1;r<=R;r++){
prod=prod*U.get(y,mu);
y=y+e_mu;
}
for(t=1;t<=M;t++){
prod=prod*U.get(y,nu);
y=y+e_nu;
}
y=y-e_mu;
for(r=1;r<=R;r++){
prod=prod*Udag.get(y,mu);
y=y-e_mu;
}
y=y+e_mu-e_nu;
for(t=1;t<=M;t++){
prod=prod*Udag.get(y,nu);
y=y-e_nu;
}

wilson[R][M]=wilson[R][M]+(1.0/NCOLOR)*Tr(prod).real();

}
}
}

for(R=1;R<=(L/2);R++){
for(M=1;M<=(L/2);M++){
f_loop << wilson[R][M]/SITES << "\t" << flush;
}
} 
f_loop  << "\n" << flush;



return;
}
