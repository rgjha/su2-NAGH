#include "measure.h"  

void measure(Gauge_Field &U){

    static int first_time=1,meas=0;
    static ofstream f_gauge,f_matter;
    double act_g,act_s;

    meas++;

    if(first_time){

    f_gauge.open("gauge");
    if(f_gauge.bad()){ 
    cout << "Failed to open gauge file\n" << flush ;}
    f_matter.open("matter");
    if(f_matter.bad()){
    cout << "Failed to open matter file\n"<< flush;}
    
    first_time=0;
	}	

	cout << "In measure \n" << flush;
	
	obs(U, act_g, act_s);
	loop(U);
	line(U);
    corrlines(U);
	     
	f_gauge << act_g << "\n" << flush;
    f_matter << act_s << "\n" << flush;
	
	return;
}
