#include "utilities.h"

// Complex type
Complex::Complex(void){re=0.0;im=0.0;}
Complex::Complex(double x, double y){re=x;im=y;}
double Complex::real(void) const{return(re);}
double Complex::imag(void) const{return(im);}
double Complex::norm(void){return(sqrt(re*re+im*im));}
void Complex::print(void) const {cout << "("<< re << ", " << im << ")";}



ostream& operator<<(ostream& out,Complex c){
	out<<c.real()<<"\t"<<c.imag();
	return out;}
istream& operator>>(istream& in,Complex & c){
	double x,y;
	in>>x>>y;
	c=Complex(x,y);
	return in;}
		 
Complex operator /(const Complex &o1, const Complex &o2){
	Complex dum;
	double norm;
	norm=o2.real()*o2.real()+o2.imag()*o2.imag();
	dum=Complex((o1.real()*o2.real()+o1.imag()*o2.imag())/norm,
		(o1.imag()*o2.real()-o1.real()*o2.imag())/norm);
	return(dum);}
Complex pow(const Complex &o1, const int o2){
	Complex c(1,0);
	for (int i=0;i<o2;i++) c=c*o1;
	return c;}
  

// unitary matrix type Umatrix
Umatrix::Umatrix(void){
	for(int i=0;i<NCOLOR;i++){
		for(int j=0;j<NCOLOR;j++) {mat[i][j]=Complex();}}
	}
	
Umatrix::Umatrix(int k){

         if(k==1){
         for(int i=0;i<NCOLOR;i++){
	        for(int j=0;j<NCOLOR;j++) {mat[i][j]=Complex();}
		}
	 for(int n=0;n<NCOLOR;n++){
	 mat[n][n]=Complex(1.0,0.0);}
	 }	 
}

	 


Umatrix::Umatrix(Complex m[NCOLOR][NCOLOR]){
	for(int i=0;i<NCOLOR;i++)
		for(int j=0;j<NCOLOR;j++) mat[i][j]=m[i][j];}
		

Complex Umatrix::get(const int i, const int j) const {return(mat[i][j]);}
void Umatrix::set(const int i, const int j, const Complex &o){mat[i][j]=o;}
void Umatrix::print(void){
	for(int i=0;i<NCOLOR;i++){
		for(int j=0;j<NCOLOR;j++){mat[i][j].print();cout << "\t";}
		cout << "\n";}
		}
			
		
ostream& operator<<(ostream& out,Umatrix s){
for(int i=0;i<NCOLOR;i++)
for(int j=0;j<NCOLOR;j++){
	out<<s.get(i,j)<<'\t';}
	return out;}		
istream& operator>>(istream& in, Umatrix & s){
	Complex v[NCOLOR][NCOLOR];
	for(int i=0;i<NCOLOR;i++)
	for(int j=0;j<NCOLOR;j++){
	in>>v[i][j];}
	s=Umatrix(v);
	return in;
	}	
		
Umatrix Adj(const Umatrix &u){
        Umatrix res;
        for(int i=0;i<NCOLOR;i++)
	for(int j=0;j<NCOLOR;j++){
	res.set(i,j,conjug(u.get(j,i)));}
	return(res);}	
	
Umatrix Cjg(const Umatrix &u){
        Umatrix res;
        for(int i=0;i<NCOLOR;i++)
	for(int j=0;j<NCOLOR;j++){
	res.set(i,j,conjug(u.get(i,j)));}
	return(res);}	
    
Umatrix gaussU(void){
        Umatrix tmp=Umatrix();
        for(int i=0;i<NUMGEN;i++){
	tmp=tmp+gasdev()*Lambda[i];}
	return(tmp);
	}

Umatrix cgaussU(void){
        Umatrix tmp=Umatrix();
        for(int i=0;i<NUMGEN;i++){
	tmp=tmp+1.0/sqrt(2.0)*(Complex(1.0,0.0)*gasdev()+Complex(0.0,1.0)*gasdev())*Lambda[i];}
	return(tmp);
	}
	
Umatrix egaussU(void){
        return(exp(0.1*gaussU()));}

Umatrix operator *(const Umatrix &o1, const Umatrix &o2){
         Umatrix r;
	 Complex dum;
	 for(int i=0;i<NCOLOR;i++)
	 for(int j=0;j<NCOLOR;j++){
	 dum=Complex();
	 for(int k=0;k<NCOLOR;k++){
	 dum=dum+o1.get(i,k)*o2.get(k,j);}
	 r.set(i,j,dum);}
	 return(r);}
	 
	 
        Umatrix operator *(const Umatrix &o1, const Complex &o2){
        Umatrix dum;
	for(int i=0;i<NCOLOR;i++)
	for(int j=0;j<NCOLOR;j++){
	dum.set(i,j,o1.get(i,j)*o2);}
	return(dum);}
        Umatrix operator *(const Complex &o2, const Umatrix &o1){
	Umatrix dum;
	for(int i=0;i<NCOLOR;i++)
	for(int j=0;j<NCOLOR;j++){
	dum.set(i,j,o1.get(i,j)*o2);}
	return(dum);}
	Umatrix operator *(const Umatrix &o1, const double o2){
        Umatrix dum;
	for(int i=0;i<NCOLOR;i++)
	for(int j=0;j<NCOLOR;j++){
	dum.set(i,j,o1.get(i,j)*o2);}
	return(dum);}
	Umatrix operator *(const double o2, const Umatrix &o1){
        Umatrix dum;
	for(int i=0;i<NCOLOR;i++)
	for(int j=0;j<NCOLOR;j++){
	dum.set(i,j,o1.get(i,j)*o2);}
	return(dum);}
	Umatrix operator +(const Umatrix &x, const Umatrix &y){
	Umatrix dum;
	for(int i=0;i<NCOLOR;i++)
	for(int j=0;j<NCOLOR;j++)
	dum.set(i,j,x.get(i,j)+y.get(i,j));
	return(dum);
	}
        Umatrix operator -(const Umatrix &x, const Umatrix &y){
	Umatrix dum;
	for(int i=0;i<NCOLOR;i++)
	for(int j=0;j<NCOLOR;j++)
	dum.set(i,j,x.get(i,j)-y.get(i,j));
	return(dum);
	}
	


    Umatrix exp(const Umatrix &u){
	Umatrix c,del,prod;
	double fac=1.0;
        int i=1;
        prod=Umatrix(1);
        c=Umatrix(1);
        static int sum=0,counter=0;
	
	do{
        fac=fac*(double)i;
        prod=prod*u;
        del=prod*(1.0/fac);
        c=c+del;
        i++;}
        while(sqrt(Tr(del*Adj(del)).real())>GAUGETOL);

        sum+=i;
        counter++;
        if(counter==100000){
        cout << "mean no. of terms in U exp() " 
        << (double)sum/counter << "\n" << flush;
        counter=0;sum=0;}
        
	return(c);}
	 
					
	
Complex Tr(const Umatrix &o){
         Complex dum=Complex();
	 for(int i=0;i<NCOLOR;i++)
	 dum=dum+o.get(i,i);
	 return(dum);}	 
 	

			
// Lattice vector on torus
	
Lattice_Vector::Lattice_Vector(void){for(int i=0;i<D;i++)coords[i]=0;}
Lattice_Vector::Lattice_Vector(int mu){for(int i=0;i<D;i++){
coords[i]=0;}coords[mu]=1;}
void Lattice_Vector::set(int i, int a){
coords[i]=a;
return;}
int Lattice_Vector::get(int i) const{
return(coords[i]);}
void Lattice_Vector::print(void) const {
for(int i=0;i<D;i++) {cout << coords[i] << "\t";} cout << "\n";}

Lattice_Vector operator +(const Lattice_Vector &x, const Lattice_Vector &y){
Lattice_Vector dum;
for(int i=0;i<(D-1);i++)
dum.set(i,(x.get(i)+y.get(i))%L);
dum.set(D-1,(x.get(D-1)+y.get(D-1))%T);
return(dum);}

Lattice_Vector operator -(const Lattice_Vector &x, const Lattice_Vector &y){
Lattice_Vector dum;
for(int i=0;i<(D-1);i++)
dum.set(i,(x.get(i)-y.get(i)+L)%L);
dum.set(D-1,(x.get(D-1)-y.get(D-1)+T)%T);
return(dum);}

Lattice_Vector operator -(const Lattice_Vector &x){
Lattice_Vector dum;
for(int i=0;i<D;i++)
dum.set(i,-1*x.get(i));
return(dum);
}

double BC(const Lattice_Vector &x, const Lattice_Vector &y){

if(x.get(D-1)+y.get(D-1)<0) return(PBC);
if(x.get(D-1)+y.get(D-1)>(T-1))return(PBC);

return(1.0);
}

Lattice_Vector operator *(const int a, const Lattice_Vector &x){
Lattice_Vector dum;
for(int i=0;i<D;i++){
dum.set(i,a*x.get(i));}
return(dum);
}

int loop_over_lattice(Lattice_Vector &x, int &site){
int test,i,current;
static int Lattice_Map[D];
static int first_time=1;
static int ixav[];
static int ixbv[];
    

if(first_time){
for(i=0;i<D;i++)
Lattice_Map[i]=(int)pow((double)L,(double)i);
first_time=0;
}
    
if(first_time){
for(i=0;i<Nxa;i++)
        ixav[i] = i
for(j=0;j<Nxb;i++)
    ixbv[j] = Nxa + j
}


current=site;
for(i=D-1;i>=0;i--){
x.set(i,current/Lattice_Map[i]);
current=current-Lattice_Map[i]*x.get(i);
}

if(current!=0){cout << "Error in loop_over_lattice" << "\n";}

if(site==SITES)
test=1;
else
test=0;

site++;
return(!test);
}




// Gauge Field

Gauge_Field::Gauge_Field(void){
for(int i=0;i<SITES;i++)
for(int j=0;j<D;j++){
link[i][j]=Umatrix();}
}

Gauge_Field::Gauge_Field(int hot){

// eg for momenta
if(hot==1){
for(int i=0;i<SITES;i++){
for(int j=0;j<D;j++){
link[i][j]=gaussU();}}
}

// Warm start etc.
if(hot==2){
for(int i=0;i<SITES;i++){
for(int j=0;j<D;j++){
link[i][j]=egaussU();}}
}

// Unit matrix
if(hot==0){
for(int i=0;i<SITES;i++)
for(int j=0;j<D;j++){
link[i][j]=Umatrix(1);}
}


return;
}


Umatrix Gauge_Field::get(const Lattice_Vector &x, const int mu) const{
int site=0,i;
static int first_time=1;
static int Lattice_Map[D];

if(first_time){
for(i=0;i<D;i++)
Lattice_Map[i]=(int)pow((double)L,(double)i);
first_time=0;
}

for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

return(link[site][mu]);
}

void Gauge_Field::set(const Lattice_Vector &x, const int mu, const Umatrix &u){
int site=0,i;
static int first_time=1;
static int Lattice_Map[D];

if(first_time){
for(i=0;i<D;i++)
Lattice_Map[i]=(int)pow((double)L,(double)i);
first_time=0;
}

for(i=0;i<D;i++)
site=site+x.get(i)*Lattice_Map[i];

link[site][mu]=u;
return;
}


	 
Gauge_Field Adj(const Gauge_Field & U){
    int site;
    Gauge_Field W;
    Lattice_Vector x;
    for(int mu=0;mu<D;mu++){
        site=0;
        while(loop_over_lattice(x,site))
            W.set(x,mu,Adj(U.get(x,mu)));
    }
    return W;
}

Gauge_Field Cjg(const Gauge_Field & U){
    int site;
    Gauge_Field W;
    Lattice_Vector x;
    for(int mu=0;mu<D;mu++){
        site=0;
        while(loop_over_lattice(x,site))
            W.set(x,mu,Cjg(U.get(x,mu)));
    }
    return W;
}




