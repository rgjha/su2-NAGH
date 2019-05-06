#include "ym.h"

int SWEEPS,GAP,THERM,SEED,READIN,TRAJECTORY_LENGTH;
double BETA,KAPPA,C,DT,TIME;

Umatrix Lambda[NUMGEN];

int main(int argc, char *argv[]){
int sweep,number=0;
Gauge_Field U;

read_param();

if(READIN){
read_in(U);
//check(U);
}
else{
U=Gauge_Field(2);
}


cout << "Warming up" << "\n" << flush;
for(sweep=1;sweep<=THERM;sweep++){
clock_t time=clock();
update(U);
cout << "time for sweep " << (double)(clock()-time)/CLOCKS_PER_SEC << endl;
write_out(U,number);}

cout << "Commencing measurement sweeps" << "\n" << flush;

for(sweep=1;sweep<=SWEEPS;sweep++){
clock_t time=clock();
update(U);
cout << "time for sweep " << (double)(clock()-time)/CLOCKS_PER_SEC << endl;

//        check(U);
 write_out(U,number);

	//  measure config
        cout << "sweep no. " << sweep << "\n" << flush;
	if(sweep%GAP==0){
        number++;
        measure(U);
//        write_out(U,V,sigma,F,number);
	}
}

return(0);
}
