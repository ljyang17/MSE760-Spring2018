c---------------------------------------------------------
c     Including file for parameters and constants used
c     in the entire software packages.
c--------------------------------------------------------- 

       double precision a,sigma,eps,a_reduced,box_size
	   integer npmax,np
	   parameter (npmax = 20000)
	   double precision x(3,npmax)
	   integer lattice_size
	   
	   common /parameters/ a,sigma,eps,a_reduced,x,lattice_size,np
	   common /parameters/ box_size
	   