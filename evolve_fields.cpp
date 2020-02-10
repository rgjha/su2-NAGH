#include "evolve_fields.h"


void evolve_fields(Gauge_Field &U,
		   Gauge_Field &p_U,
		   Gauge_Field &f_U
		   ){

Lattice_Vector x;
int mu,site;

update_gauge_momenta(p_U,f_U,0.5*DT);
update_gauge_field(U,p_U,DT);
force(U,f_U);
update_gauge_momenta(p_U,f_U,0.5*DT);


return;
}
