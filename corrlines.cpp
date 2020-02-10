#include "corrlines.h"

void corrlines(const Gauge_Field &U){

static int first_time=1;
int site=0,site2=0,t,i;
Complex corr[L];
static ofstream f_line2;
Lattice_Vector x,y,e_mu;
Gauge_Field prod;

 
  	if(first_time==1){
        f_line2.open("corrlines");
	if(f_line2.bad()){
	cerr << "failed to open corrlines file" << "\n";exit(1);}
	first_time=0;}



for(i=0;i<L;i++){
corr[i]=Complex();
}

prod=Gauge_Field(0);
site=0;
while(loop_over_lattice(x,site)){

e_mu=Lattice_Vector(D-1);
y=x;
for(t=1;t<=T;t++){
prod.set(x,D-1,prod.get(x,D-1)*U.get(y,D-1));
y=y+e_mu;
}

}


for(int mu=0;mu<(D-1);mu++){
	
site=0;
while(loop_over_lattice(x,site)){
site2=0;
while(loop_over_lattice(y,site2)){
if((x.get(D-1)-y.get(D-1))==0){

t=(x.get(mu)-y.get(mu)+L)%L;
corr[t]=corr[t]+(1.0/((D-1)*NCOLOR*NCOLOR))*Tr(prod.get(x,D-1))*Tr(Adj(prod.get(y,D-1)));
}

}}}


for(i=0;i<L;i++){
f_line2 << i << "\t" <<  corr[i].real()*(1.0/(SITES))<< "\n";}
f_line2 << flush;

return;
}
