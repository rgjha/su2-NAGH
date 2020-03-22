#!/usr/bin/python3
# This script does Monte Carlo
# for SU(2) non-Abelian gauge Higgs model in 2d. 
# Work in progress!

import numpy as np
import cmath
import sys 
import datetime as dt
from itertools import product
from math import * 
from scipy import *
from scipy.linalg import expm
from numpy import linalg as LA
import os 

# .....................................
NCOL = 2   # Rank of gauge group 
LEN = 4  # Sites
BETA = 6.0
KAPPA = 1.0 
EPS = 0.01  # 
GENS = NCOL**2 - 1 
TRAJ_LENGTH = int(0.01/EPS)
Niters_sim = 100
DIM = 2 
HAM = []
expDS = [] 
ACT = [] 
MOM = []
Nxa = LEN-2
Nxb = LEN-2
Nxa2 = Nxa-1
Nxb2 = Nxb-1
Nrep = 2
alpha = 0.2
Ntr = Nrep*LEN 

 
# Average of column file
# awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' gauge 

link = mom  = U = L = U_bak = f_U = np.zeros((LEN, LEN, DIM, NCOL, NCOL), dtype=complex)
Lambda = np.zeros((GENS, NCOL, NCOL), dtype=complex)
ixav = ixbv = ixav2 = ixbv2 = np.zeros((Nxa), dtype=integer)


if len(sys.argv) < 3:
  print("Usage:", str(sys.argv[0]), "READIN " " SAVE_or_NOT ? ")
  sys.exit(1)

READIN = int(sys.argv[1])
SAVE = int(sys.argv[2])


def pretty_print_matrix(matrix):
    for row in matrix:
       for col in row:
           print("{:8.3f}".format(col), end=" ")
print("")


# .......Boundary conditions...........
def bc(x):
    return x%LEN 
# .....................................


# .......Map local to global...........
def do_map(ix):

    for i in range (Nxa):
        ixav[i] = ix
    for j in range (Nxb):
        ixbv[i] = Nxa + ix 
    for k in range (Nxa2):
        ixav2[i] = ix
    for l in range (Nxb2):
        ixbv2[i] = Nxa2 + ix

    return None  
# .....................................


# .....................................
# Set up generator matrices with normalization, Tr(TaTb) = -delta_ab
# Switching to Tr(TaTb) = delta_ab requires multiplying all Lambda's by "i"
# Nice reference arXiv:1310.5353

def setup_generators():
    
    inv_sqrt = complex(1.0 / math.sqrt(2.0), 0.0)
    i_inv_sqrt = complex(0.0, 1.0 / math.sqrt(2.0))
    count = 0
    
    for i in range(NCOL):
        for j in range(i + 1, NCOL):
            for k in range(NCOL):
                for l in range(0, NCOL):
                    
                    if k == i and l == j:
                        Lambda[count][k][l] += i_inv_sqrt
                        Lambda[count + 1][k][l] += inv_sqrt
                    
                    elif k == j and l == i:
                        Lambda[count][k][l] += i_inv_sqrt
                        Lambda[count + 1][k][l] -= inv_sqrt
            count += 2

    if count != NCOL * (NCOL - 1):
        print ("Count is", count)
        print('Wrong number of off-diagonal generators')
        sys.exit(1)
    
    for i in range(NCOL - 1):
        j = NCOL * (NCOL - 1) + i
        k = i + 1
        
        i_inv_sqrt = complex(0.0, 1.0 / math.sqrt(k * (k + 1.0)))
        
        for p in range(k + 1):
            Lambda[j][p][p] = i_inv_sqrt
        
        Lambda[j][k][k] *= -k


    for j in range (GENS):
        
        if np.trace(np.dot(Lambda[j],Lambda[j]))-1.0 > 1e-14:
            print(("Tr[TaTb] is incorrect ", np.trace(np.dot(Lambda[j],Lambda[j]))))

    return Lambda

# .....................................


# .....................................
def dagger(a):

    return np.transpose(a).conj()
# .....................................


# .....................................
def cold_start():

    for dir in range (DIM): 
        for i in range (LEN):
            for j in range (LEN):
                link[i][j][dir] = np.eye(NCOL)

    return link 
# .....................................


# ...Return product of many matrices.....
def mm(*mats):
    product = np.eye(NCOL) 
    for i in mats:
        product = np.dot(product, i)
    return product
# .....................................


# .....................................
def copy_fields(a):  
    
    for t in range (LEN):
        for x in range (LEN):
            for dir in range (DIM): 

                U_bak[x][t][dir]  = a[x][t][dir]  
        
    return U_bak

# .....................................


# .....................................
def rejected_go_back_old_fields(a):
    
    for t in range (LEN):
        for x in range (LEN):
            for dir in range (DIM): 
                
                U[x][t][dir] = a[x][t][dir]


    return U
# .....................................


def refresh_mom():
    
     
    for i in range (LEN):
        for j in range (LEN):
            for dir in range (DIM):

                mom[i][j][dir] = random_anti_hermitian2()

                # Check AH (if needed!)
                # print (mom[i][j][dir] + dagger(mom[i][j][dir]))
        

    return mom

# .....................................

# Construct random Anti-hermitian (AH) matrix

def random_anti_hermitian(): 

    tmp = np.zeros((NCOL, NCOL), dtype=complex)
 
    for i in range (NCOL):
        for j in range (NCOL):

            if i==j:
                tmp[i][j] = 0.0 
            else:
                tmp[i][j] = complex(np.random.normal(0,1), np.random.normal(0,1))/sqrt(2) 

    tmp = (tmp + dagger(tmp))*0.50

    for i in range (NCOL):
        for j in range (NCOL):

            tmp[i][j] = complex(tmp[i][j].imag,tmp[i][j].real)


    return tmp

# .....................................



def random_anti_hermitian2():

    tmp = np.zeros((NCOL, NCOL), dtype=complex)

    for i in range (GENS):

        tmp = tmp + Lambda[i]* np.random.normal(0, 1)


    return tmp

# .....................................


def kinetic_energy(mom):

    s1 = 0.0  

    for t in range (LEN):
        for x in range (LEN):
            for dir in range (DIM): 
                s1 -= 0.50 * np.trace(np.dot(mom[x][t][dir], mom[x][t][dir]))

    return (s1.real)  

# .....................................


def force(U):


    tmp = np.zeros((LEN, LEN, DIM, NCOL, NCOL), dtype=complex)
    staple = np.zeros((LEN, LEN, DIM, NCOL, NCOL), dtype=complex)
    staple3 = np.zeros((LEN, LEN, DIM, NCOL, NCOL), dtype=complex)

    for t in range (LEN):
        for x in range(LEN):
            for dir in range (DIM):  # \mu direction

                for dir2 in range (DIM):  # \nu direction 
                    if dir2 == dir:
                        continue
                    else:
                         
                        staple[x][t][dir] += mm(U[bc(x+dir2)][bc(t+dir)][dir2],dagger(U[bc(x+dir)][bc(t+dir2)][dir]),dagger(U[bc(x)][bc(t)][dir2])) \
                        + mm(dagger(U[bc(x+dir2-dir)][bc(t-dir2+dir)][dir2]),dagger(U[bc(x-dir)][bc(t-dir2)][dir]),U[bc(x-dir)][bc(t-dir2)][dir2])

                                    
                staple[x][t][dir] = np.dot(U[bc(x)][bc(t)][dir],staple[x][t][dir]) - dagger(np.dot(U[bc(x)][bc(t)][dir],staple[x][t][dir]))                
                staple[x][t][dir] = staple[x][t][dir] - (1/NCOL)*np.trace(staple[x][t][dir])*np.eye(NCOL) 
                f_U[x][t][dir] = (-BETA/(2.0*NCOL)) * staple[x][t][dir]
                #print ("Force @ ", x, t , dir, pretty_print_matrix(f_U[x][t][dir])) 
                # For debugging force


    for t in range (LEN):
        for x in range(LEN):
            for dir in range (DIM):


                staple3[x][t][dir] = U[bc(x)][bc(t)][dir] - dagger(U[bc(x)][bc(t)][dir])
                staple3[x][t][dir] = staple3[x][t][dir] - (1/NCOL)*np.trace(staple3[x][t][dir])*np.eye(NCOL) 
                f_U[x][t][dir] = f_U[x][t][dir] - (KAPPA/4.0)*staple3[x][t][dir]
 
    return f_U

# .....................................

def leapfrog(U,f_U,mom,eps):


    mom = update_gauge_momenta(mom, f_U, eps)
    U = update_gauge_field(U, mom, eps)
    f_U = force(U)
    mom = update_gauge_momenta(mom, f_U, eps)


    return U, mom, f_U

# .....................................

def update_gauge_momenta(mom, f_U, eps):


    for t in range (LEN):
        for x in range(LEN):
            for dir in range (DIM):  
                mom[x][t][dir] += eps*f_U[x][t][dir]

    return mom 

# .....................................



# .....................................

def update_gauge_field(U, mom, eps):


    for t in range (LEN):
        for x in range(LEN):
            for dir in range (DIM):  

                U[x][t][dir] = np.dot(expm(mom[x][t][dir]*eps), U[x][t][dir]) 

    return U 

# .....................................


# .....................................
# Calculates the action 
# .....................................
def action(L):


    tmp = 0.0 
    tmp1 = 0.0 

    for t in range (LEN):
        for x in range(LEN):
          
            tmp -= np.real(np.trace(mm(L[bc(x)][bc(t)][0],L[bc(x+1)][bc(t)][1],dagger(L[bc(x)][bc(t+1)][0]),dagger(L[bc(x)][bc(t)][1]))))

            for dir in range (DIM):
                tmp1 -= np.real(np.trace(L[bc(x)][bc(t)][dir]))
                #print ("Action  @ ", x ,  t , (BETA/NCOL)*tmp + ((KAPPA/2.0)*tmp1))  

                # Dictionary from this code to C++ code on GitHub 
                # L[bc(x)][bc(t)][0] <-> U.get(x,mu)
                # L[bc(x+1)][bc(t)][1] <-> U.get(x+e_mu,nu)
                # L[bc(x)][bc(t+1)][0]^dagger <-> Udag.get(x+e_nu,mu)
                # L[bc(x)][bc(t)][1]^dagger <-> Udag.get(x,nu)


    return ((BETA/NCOL)*tmp) + ((KAPPA/2.0)*tmp1)  
# .....................................


def update(U):

    mom = refresh_mom()
    KE = kinetic_energy(mom)
    ba = action(U)
    print ("Action is:", ba)
    
    f_U = force(U)
    sys.exit(1)
    start_act =  ba + KE
    #print("start action: kinetic+comm. = " , ba , "and gauge+scalar momenta =" , KE)
    U_bak = copy_fields(U) 

    for i in range (TRAJ_LENGTH):
        U, mom, f_U = leapfrog(U, f_U, mom, EPS) 
    

    KE = kinetic_energy(mom)
    ba = action(U)
    end_act = ba + KE
    #print("end action: kinetic+comm. = " , ba , "and gauge+scalar momenta =" , KE)
    change = end_act - start_act

    HAM.append(abs(change))
    expDS.append(np.exp(change))


    expDS.append(np.exp(-1.0*change))

    if np.exp(-change) < random.uniform(0,1):
        U  = rejected_go_back_old_fields(U_bak)
        print(("REJECT: deltaS = " "%8.7f " " startS = " "%8.7f" " endS = " "%8.7f" % (change, start_act, end_act))) 
    else:   
        print(("ACCEPT: deltaS = " "%8.7f " "startS = " "%8.7f" " endS = " "%8.7f" % (change, start_act, end_act)))

    return U, ba 



if __name__ == "__main__":

    print ("2d Wilson type SU(2) action on " "%3.0f " " x" " %3.0f "  "lattice " %(LEN,LEN))
    print ("Gauge coupling (BETA) = " "%3.3f " "& " " Matter coupling (KAPPA) " " %3.3f " %(BETA,KAPPA))
    print ("Unitary gauge fixed for all links")

    Lambda = setup_generators() 

    if READIN ==1:

        file_length = 0

        with open('CONFIG.txt') as infp:
            for line in infp:
                if line.strip():
                    file_length += 1

        if (file_length != LEN*LEN*DIM):
            print ("Error in reading configuration file!")
            sys.exit(1) 
        else: 
            print ("Reading old config.")
            f = open("CONFIG.txt", "r")
        
            for t in range (LEN):
                for x in range (LEN):
                    for dir in range (DIM):

                        dum = f.readline().split()

                        for i in range (NCOL):
                            for j in range (NCOL):

                                L[x][t][dir][i][j] = complex(float(dum[(4*i)+(2*j)]), float(dum[(4*i)+(2*j)+1])) 

            f.close()


    else:
        L = cold_start() 


    
    t0 = dt.datetime.now()
    f1 = open("gauge", "w") 
    
    for i in range (Niters_sim):
        L, ba = update(L)

        if i%2 == 1:
                f1.write("%s\n" % ba)  

        if READIN ==1 and i%10 == 0:

            if os.path.exists('CONFIG.txt'):
                os.remove('CONFIG.txt')

            f3 = open("CONFIG.txt", "w")

            for t in range (LEN):
                for x in range (LEN):
                    for dir in range (DIM):
                        np.savetxt(f3, L[x][t][dir].view(double),fmt='%4e', newline="   ") 
                        f3.write("\n")

            f3.close() 


 
    f1.close()


    t1 = dt.datetime.now()
    print('\nTime total: %10.3f seconds' %(t1 - t0).total_seconds())


