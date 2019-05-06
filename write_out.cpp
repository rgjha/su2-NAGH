#include "write_out.h"

void write_out(const Gauge_Field &u, const int num){

ofstream f_write;
Lattice_Vector x;
int site;
char s[100];

sprintf(s,"config%d",num);

//f_write.open(s);
f_write.open("config");

if(f_write.bad()){
cout << "error opening config file\n" << flush;}

f_write << L << "\t" << T << "\t" << BETA;
f_write << "\t" << NCOLOR << "\t" << DT << endl;


site=0;
while(loop_over_lattice(x,site)){
for(int mu=0;mu<D;mu++){
f_write << u.get(x,mu) << "\n";
}}


f_write.close();
return;
}

