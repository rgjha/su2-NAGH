#include "update.h"

void update(Gauge_Field &U){

Gauge_Field p_U,old_U,old_f_U;
double H_old,K_old,K_new,S_new,H_new,hmc_test;
Lattice_Vector xx=Lattice_Vector();

static Gauge_Field f_U;

static double S_old=0.0;
static int first_time=1,accept=0,no_calls=0;
static ofstream f_hmc;


if(first_time){
f_hmc.open("hmc_test");
if(f_hmc.bad()){
cout << "failed to open hmc_test file\n" << flush;}


S_old=action(U);
force(U,f_U);
cout << "Trajectory length " << TRAJECTORY_LENGTH << "\n" ;
}	



if((no_calls%50==0)&&(!first_time)){
cout << "acceptance rate " << (double)accept/(double)no_calls << "\n" <<
flush;
no_calls=0;
accept=0;
}

first_time=0;
no_calls++;


// Refresh momenta
p_U=Gauge_Field(1);


K_old=kinetic_energy(p_U);
 

//Save copies of fields
old_U=U;
old_f_U=f_U;

H_old=S_old+K_old;

cout << "old H " << H_old << "\n" << flush;

//Classical evolution

for(int i=0;i<TRAJECTORY_LENGTH;i++)
evolve_fields(U,p_U,f_U);

 S_new=action(U);

 K_new=kinetic_energy(p_U);

H_new=S_new+K_new;
cout << "new H " << H_new << "\n" << flush;

hmc_test=exp(-H_new+H_old);

f_hmc << hmc_test << "\n" << flush;
      	    
//Metropolis test
	if(myrandom()<exp(H_old-H_new)){
		S_old=S_new;
		accept++;
        cout << "ACCEPT " << endl;
        
		return;
	}
	else{
    //cout << "hmc_test " << hmc_test << " failed\n" << flush;
    cout << "REJECT " << endl;
	// if fails copy back fields
	U=old_U;
	f_U=old_f_U;
	
	return;
	}	   



return;	    
}	    
	    
	    
	    
	    
	    
