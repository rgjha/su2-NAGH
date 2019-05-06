#include "setup.h"

void setup(void){
int i,j,k;
double f[NUMGEN][NUMGEN][NUMGEN];

// note: Tr(TaTb)=-delta_ab

cout << "Computing generators for SU(NCOLOR)\n" << flush;
(void)my_gen();

Complex trace;
// test tracelessness
for(i=0;i<NUMGEN;i++){
trace=Tr(Lambda[i]);
if(trace.norm()>0.000001){
cout << "TrT_"<<i<<"is " << trace << "\n";}
}

// test orthogonality

for(i=0;i<NUMGEN;i++){
for(j=0;j<NUMGEN;j++){
trace=Tr(Lambda[i]*Lambda[j]);
if(trace.norm()>0.000001){
cout << "TrT_"<< i << "T_" << j << "= " <<  trace << "\n" <<
flush;}
}}

return;
}

