#include "action.h"

void obs(const Gauge_Field &U, double & act_g, double & act_s){

Lattice_Vector x,e_mu,e_nu;

int mu,nu,site;
act_g=0.0;
// Wilson gauge action
site=0;
while(loop_over_lattice(x,site)){

for(mu=0;mu<D;mu++){
    e_mu=Lattice_Vector(mu);
    for(nu=mu+1;nu<D;nu++){
        e_nu=Lattice_Vector(nu);

Umatrix s1 = U.get(x,mu)*U.get(x+e_mu,nu)*Adj(U.get(x+e_nu,mu))*Adj(U.get(x,nu));
act_g = act_g -(1.0/2)*Tr(s1).real();
}
}
}
  
site=0;
act_s=0.0;
while(loop_over_lattice(x,site)){
for(mu=0;mu<D;mu++){
act_s=act_s-(1.0/2)*Tr(U.get(x,mu)).real();
}}

return;
}
