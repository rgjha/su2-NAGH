#include "kinetic_energy.h"

double kinetic_energy(const Gauge_Field &p_U){
Complex dum=Complex();
int sites,mu,b;
Lattice_Vector x;

// Note minus sign to take into account
// anti-hermitian (AH) generators (Lambdas)
 
 
 sites=0;
 while(loop_over_lattice(x,sites)){
 for(mu=0;mu<D;mu++){
 dum=dum-0.5*Tr(p_U.get(x,mu)*p_U.get(x,mu));
 }}
 

if(fabs(dum.imag())>0.0000000001){cout << "k.e not real!" << dum.imag()<<"\n" << flush;}
 
   
return(dum.real());
	
} 
