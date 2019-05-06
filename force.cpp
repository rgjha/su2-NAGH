#include "force.h"


void force(const Gauge_Field &U, Gauge_Field &f_U){
  Lattice_Vector x,e_mu,e_nu,xx=Lattice_Vector();
int sites,mu,nu;
 Umatrix staple,staple3;
 Gauge_Field Udag,utmp,Wdag,wtmp;

//cout << "in force" << endl;

Udag=Adj(U);
f_U=Gauge_Field();

// gauge force - Wilson contribution

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<D;mu++){
e_mu=Lattice_Vector(mu);

staple=Umatrix();

for(nu=0;nu<D;nu++){
e_nu=Lattice_Vector(nu);
if(nu==mu) continue;

staple=staple+U.get(x+e_mu,nu)*Udag.get(x+e_nu,mu)*Udag.get(x,nu)+
       Udag.get(x+e_mu-e_nu,nu)*Udag.get(x-e_nu,mu)*U.get(x-e_nu,nu);
}

staple=U.get(x,mu)*staple-Adj(U.get(x,mu)*staple);
staple=staple-(1.0/NCOLOR)*Tr(staple)*Umatrix(1);

//cout << "trace is " << Tr(staple) << "\n" << flush;

f_U.set(x,mu,-(BETA/(2*NCOLOR))*staple);
}
}

sites=0;
while(loop_over_lattice(x,sites)){
for(mu=0;mu<D;mu++){
staple3=U.get(x,mu)-Udag.get(x,mu);
staple3=staple3-(1.0/NCOLOR)*Tr(staple3)*Umatrix(1);
f_U.set(x,mu,f_U.get(x,mu)-(KAPPA/4)*staple3);
}}


return;
}
