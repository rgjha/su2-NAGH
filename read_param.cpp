#include "read_param.h"

#include <stdio.h>
#include <sys/time.h>

void read_param(void){
ifstream f_in("parameters");
double LAMBDA;
if(!f_in.good()){
	cout << "\ncan't open file parameters to read data!\n";
	exit(1);}
 f_in>>SWEEPS>>THERM>>GAP>>C>>KAPPA>>DT>>READIN;

TRAJECTORY_LENGTH=(int)(1.0/DT);

BETA=SITES*C;

cout << "Dimension " << D << "\n";
cout << "Number of colors " << NCOLOR <<  "\n";
cout << "Spatial lattice size " << L << "\n";
cout << "Temporal extent " << T << "\n"<<"\n";
cout << "Gauge coupling " << BETA << "\n";
cout << "Matter coupling " << KAPPA << "\n";
cout << "Thermalization sweeps " << THERM << "\n";
cout << "Number of sweeps " << SWEEPS << "\n";
cout << "Sweeps between measurements " << GAP << "\n";
cout << "Time step in leapfrog eqs " << DT << "\n";
    


cout << "Reading initial config: (1 for yes, 0 for no) " << READIN << "\n";

if (PBC==1.0) {cout << "periodic temporal bc for fermions" << "\n";}
else{cout << "antiperiodic temporal bc for fermions" << "\n";}

srand(random_seed());
//srand(0);

setup();

return;
}
