#include "action.h"


double action(const Gauge_Field &U){

Lattice_Vector x,e_mu,e_nu,xx=Lattice_Vector();
double act_g=0.0,act_s=0.0;
int mu,nu,site;
Gauge_Field Udag;
Umatrix s1,w1;



Udag=Adj(U);

// Wilson gauge action
site=0;
while(loop_over_lattice(x,site)){

for(mu=0;mu<D;mu++){
e_mu=Lattice_Vector(mu);

for(nu=mu+1;nu<D;nu++){
e_nu=Lattice_Vector(nu);

	s1 = U.get(x,mu)*U.get(x+e_mu,nu)*Udag.get(x+e_nu,mu)*Udag.get(x,nu);
	act_g = act_g - BETA/NCOLOR*Tr(s1).real();
}
}
}

site=0;
w1=0.0;

while(loop_over_lattice(x,site)){
for(mu=0;mu<D;mu++){
  act_s=act_s-(KAPPA/2)*Tr(U.get(x,mu)).real();}}


return(act_g+act_s);
}
