
import matplotlib.pyplot as plt
import math as m
import numpy as np
import random

a = 5.26e-10; # meter
sig = 3.4e-10; # meter
a_1 = a/sig; # reduced a
eps = 0.0104; # eV

def initialize_lattice(X):
    N = (X**3)*4; # Number of atoms
    # coordinates of the 4 base atoms in a unit cell
    xunit = [0, 0.5, 0.5, 0];
    yunit = [0, 0.5, 0, 0.5];
    zunit = [0, 0, 0.5, 0.5];

    # generating lattice
    r = np.zeros(shape=(N,3));
    x = 0;
    y = 0;
    z = 0;
    k = 0;
    n = 0;
    for x in range(X):
        for y in range(X):
            for z in range(X):
                for k in range(4):
                    if n < N:
                        r[n, 0] = (x + xunit[k])*a_1;
                        r[n, 1] = (y + yunit[k])*a_1;
                        r[n, 2] = (z + zunit[k])*a_1;
                        n = n + 1;
    # return the coordinated of the position of the atoms
    return r;
    
def forces(r, X):
    N = (X**3)*4; 
    i = 0;
    m = 0;
    s = 0;
    p = int((N*(N-1))/2);
    rij = np.zeros(shape=(p,3));
    rij_2 = np.zeros(shape=(p,1));
    alpha = np.zeros(shape=(p,3));
    a = np.zeros(shape=(N,3));
    for i in range(0, N-1):
        for j in range(i+1,N):
            for m in range(3):
                rij[s, m] = r[i, m] - r[j, m];
                # applying periodic boundary condition
                if rij[s, m] < (-(X*a_1)/2):
                    rij[s, m] = rij[s, m] + (X*a_1);
                if rij[s, m] >= ((X*a_1)/2):
                    rij[s, m] = rij[s, m] - ((X*a_1));
                else:
                    rij[s, m] = rij[s, m];
            rij_2[s,0] = np.dot(rij[s,:],rij[s,:]);
            alpha[s,:] = rij[s,:]*((48.0/rij_2[s,0])*(((1/rij_2[s,0])**6.) - (1/(2.0*(rij_2[s,0])**3.))));
            a[i,:] = a[i,:] + alpha[s,:];
            a[j,:] = a[j,:] - alpha[s,:];
            s = s + 1;

  #  rij_2 = np.zeros(shape=(p,1));
  #  uij = np.zeros(shape=(p,1));
  #  alpha = np.zeros(shape=(p,3));
  #  a = np.zeros(shape=(N,3));
 #   for i in range(0, N-1):
  #      for j in range(i+1,N):
   #         for c in range(p):
            #    uij[c,0] = 4.*((1/(rij_2[c,0]))**6. - (1/(rij_2[c,0]))**3.);
       # if c<50: #calculating forces for 50 atom pairs
           #     alpha[c,:] = rij[c,:]*((48/rij_2[c,0])*(((1/rij_2[c,0])**6.) - (1/(2*(rij_2[c,0])**3.))));
        #    a[i,:] = a[i,:] + alpha[]
 #   a_units = a*(eps/sig)*(1.60217733e-19);
 #   np.savetxt('forces.txt', a_units.reshape(50,3))
   # u = np.sum(uij, axis=0);
   # u_avg_pbc = (u*eps) / N;

    return a;

    
def initialize_velocity(X, Temp):
    N = (X**3)*4;
    kb = 1.3807e-23; # joule per kelvin 
    T_1 = (Temp*kb)/(1.65e-21);
    v = np.zeros(shape=(N,3));
    v_2 = np.zeros(shape=(N,1));
    com_v = np.zeros(shape=(1,3));
    xi = np.zeros(shape=(N,3));
    v_o = m.sqrt(3*T_1);
    n = 0;
    for i in range(0,N+1000):
        if n < N:
            xi_0 = 2.*random.uniform(0,1) - 1;
            xi_1 = 2.*random.uniform(0,1) - 1;
            s_2 = xi_0**2 + xi_1**2;
            if s_2 < 1.0:
                xi[n,0] = 2.0*m.sqrt(1.-s_2)*xi_0;
                xi[n,1] = 2.0*m.sqrt(1.-s_2)*xi_1;
                xi[n,2] = 1.-(2.0*s_2);
                v[n,:] = v_o*xi[n,:];
                #com_v = com_v + v[n,:];
                v_2[n,0] = np.dot(v[n,:],v[n,:]);
                n = n + 1;
    #for n in range(N):
        #v[n,:] = v[n,:] - (com_v/N);

    cur_ke = 0.0;
    for n in range(0,N):
        cur_ke = cur_ke + (np.dot(v[n,:],v[n,:]));
    cur_temp = (cur_ke)/(3.0*N);
    #scale_factor = m.sqrt(Temp/cur_temp);

    #new_ke = 0.0;
    #for n in range(N):
    #    v[n,:] = v[n,:]*scale_factor;
    #    new_ke = new_ke + (0.5*np.dot(v[n,:],v[n,:])); 
    #rescale_temp = (0.66*new_ke)/N;
    print(cur_temp*119.5)
    
    #T = (1/(3.0*N))*np.sum(v_2, axis=0)
    #print(T*119.5) # to check the temp is correct
    #np.savetxt('v.txt', v.reshape(N,3))
    return v;
    
def LJ_pbc(r, X):
    N = (X**3)*4; 
    i = 0;
    m = 0;
    s = 0;
    p = int((N*(N-1))/2);
    rij = np.zeros(shape=(p,3))
    for i in range(0, N-1):
        for j in range(i+1,N):
            for m in range(3):
                rij[s, m] = r[i, m] - r[j, m];
                # applying periodic boundary condition
                if rij[s, m] < (-(X*a_1)/2):
                    rij[s, m] = rij[s, m] + (X*a_1);
                if rij[s, m] >= ((X*a_1)/2):
                    rij[s, m] = rij[s, m] - ((X*a_1));
                else:
                    rij[s, m] = rij[s, m];
            s = s + 1;

    rij_2 = np.zeros(shape=(p,1));
    inv_rij_2 = np.zeros(shape=(p,1));
    uij = np.zeros(shape=(p,1));
    for c in range(p):
        rij_2[c,0] = np.dot(rij[c,:],rij[c,:]);
        inv_rij_2[c,0] = 1.0/(rij_2[c,0]);
        uij[c,0] = (inv_rij_2[c,0]**6. - inv_rij_2[c,0]**3.);
    
    u = np.sum(uij, axis=0);
    u_avg_pbc = (4.0*u*eps) / N;

    return u_avg_pbc;

def KE(v,X):
    N = (X**3)*4; 
    v_2 = np.zeros(shape=(N,1));
    for c in range(N):
        v_2[c,0] = np.dot(v[c,:],v[c,:]);
    
    k_tot = np.sum(v_2, axis=0);
    k_val = (0.5*k_tot*eps) / N;

    return k_val;
    
def MD(X, Temp, dt, t):
    n = int(t / dt);
    N = (X**3)*4;
    r = initialize_lattice(X);
    v = initialize_velocity(X, Temp);
    a = forces(r, X);
    u = np.zeros(shape=(n,1));
    k = np.zeros(shape=(n,1));
    e = np.zeros(shape=(n,1));
    for nstep in range(0,n):
        for i in range(N):
            v[i,:] = v[i,:] + (a[i,:]*(dt/2));
            r[i,:] = r[i,:] + (v[i,:]*dt);
            for b in range(3):
                if r[i, b] < (-(X*a_1)/2): # check
                    r[i, b] = r[i, b] + (X*a_1);
                if r[i, b] >= ((X*a_1)/2):
                    r[i, b] = r[i, b] - ((X*a_1));
        a = forces(r, X);
        for i in range(N):
            v[i,:] = v[i,:] + (a[i,:]*(dt/2.));
        u[nstep,:] = LJ_pbc(r,X);
        k[nstep,:] = KE(v,X);
        e[nstep,:] = k[nstep,:] + u[nstep,:];
    np.savetxt('e.txt', e.reshape(n,1))    
    time = np.zeros(shape=(n,1))
    for ncount in range(1,n):
        time[ncount,0] = ncount*dt;
    #plt.figure()
    pe, = plt.plot(time,u, label='Potential Energy') 
    ke, = plt.plot(time,k, label='Kinetic Energy')
    te, = plt.plot(time,e, label='Total Energy') 
    plt.xlabel('time (LJ units)')
    plt.ylabel('Energy (eps)')
    plt.legend(handles=[pe, ke, te], loc=1)
    plt.savefig('dt_001.jpeg')
    plt.close()  
    return e;      
u_001 = MD(3, 300, 0.001, 0.01);
