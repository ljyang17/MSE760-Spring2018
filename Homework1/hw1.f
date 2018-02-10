      program HW1
	  implicit none
	  include 'param.inc'
	  
	  character(len=64) :: cmd, filename
	  logical :: filefound
	  double precision total_potential,temp_potential,avg_potential
	  double precision temp
	  integer i,j,k
	  ! initialize parameters
	  
	  a = 5.26d-10
	  sigma = 3.4d-10
	  eps = 0.0104
	  a_reduced = a/sigma
	  lattice_size = 15
	  
	  ! initialize parameters
	  
	  ! delete old data file
	  
	  filename = "./position.vtk"
      inquire (file=filename, exist=filefound)
      if (filefound) then
        write (cmd, '("/bin/rm ", A)' ) trim (filename)
        call system (cmd)
      endif
	  
	  ! end delete old data file
	  call initialize_lattice()
	  call output_data()
	  
	  ! find box size
	  box_size = 0.0d0
	  do i = 2,np
	    temp = abs(x(1,1)-x(1,i))
		if(temp.gt.box_size) then
		  box_size = temp
		endif
	  enddo
	  
	  !write(*,*) box_size
	  !pause
	  
	  ! end find box size
	  
	  ! potential without PBC
	  
	  total_potential = 0.0d0
	  
	  do i = 1,np-1
	    do j = i+1,np
		  call LJ_potential(x(:,i),x(:,j),temp_potential)
		  total_potential = total_potential + temp_potential
		enddo
	  enddo
	  
	  total_potential = total_potential*eps
	  avg_potential = total_potential/np
	  
	  write(*,*) 'Total potential NO PBC', total_potential
	  write(*,*) 'Average potential NO PBC', avg_potential
	  
	  ! end potential without PBC
	  
	  ! potential with PBC
	  
	  total_potential = 0.0d0
	  
	  do i = 1,np-1
	    do j = i+1,np
		  call LJ_potential_PBC(x(:,i),x(:,j),temp_potential)
		  total_potential = total_potential + temp_potential
		enddo
	  enddo
	  
	  total_potential = total_potential*eps
	  avg_potential = total_potential/np
	  
	  write(*,*) 'Total potential PBC', total_potential
	  write(*,*) 'Average potential PBC', avg_potential
	  
	  ! end potential with PBC
	  
	  end
	  
	  subroutine initialize_lattice()
	  implicit none
	  include 'param.inc'
	  
	  ! define parameters
	  
	  double precision parents_points(3,4)
	  integer i,j,k,kk,counter
	  
	  ! end define parameters
	  
	  ! initialize parents points
	  
	  parents_points(1,1) = 0.0d0*a_reduced
	  parents_points(2,1) = 0.0d0*a_reduced
	  parents_points(3,1) = 0.0d0*a_reduced
	  
	  parents_points(1,2) = 0.5d0*a_reduced
	  parents_points(2,2) = 0.5d0*a_reduced
	  parents_points(3,2) = 0.0d0*a_reduced
	  
	  parents_points(1,3) = 0.0d0*a_reduced
	  parents_points(2,3) = 0.5d0*a_reduced
	  parents_points(3,3) = 0.5d0*a_reduced
	  
	  parents_points(1,4) = 0.5d0*a_reduced
	  parents_points(2,4) = 0.0d0*a_reduced
	  parents_points(3,4) = 0.5d0*a_reduced
	  
	  ! end initialize parents points
	  
	  counter = 0
	  do i = 1,lattice_size
	    do j = 1,lattice_size
		  do k = 1,lattice_size
		    do kk = 1,4
			  counter = counter + 1
			  x(1,counter) = parents_points(1,kk) + (i-1)*a_reduced
			  x(2,counter) = parents_points(2,kk) + (j-1)*a_reduced
			  x(3,counter) = parents_points(3,kk) + (k-1)*a_reduced
			enddo
		  enddo
		enddo
	  enddo
	  
	  np = counter
	  
	  if(counter.gt.npmax) then
	  
	    write(*,*) 'need bigger array for data storage'
		write(*,*) '# of points', counter
		write(*,*) 'npmax', npmax
		
	  endif
	  
	  
	  end
	  
	  subroutine output_data()
	  implicit none
	  include 'param.inc'
	  integer k
	  
	  open(79,file="./position.vtk",position="append")
	  write(*,*) "Writing data into file"
	  
	  write(79,*) '# vtk DataFile Version 2.0'
	  write(79,*) 'particlePositions'
	  write(79,*) 'ASCII'
	  write(79,*) 'DATASET POLYDATA'
	  write(79,*) 'POINTS ', np, ' float'
	  
	  do k = 1,np
	    write(79,*) x(1,k), x(2,k), x(3,k)
      enddo
	  
	  close(79)
	  end
	  
	  subroutine LJ_potential (p1, p2, potential)
	  implicit none
	  include 'param.inc'
	  
	  double precision p1(3),p2(3),potential
	  double precision rr
	  integer k
	  
	  rr = 0.0d0
	  do k = 1,3
	    rr = rr + (p1(k)-p2(k))**2.0d0
	  enddo
	  rr = sqrt(rr)
	  potential = 4.0d0*(rr**(-12.0d0)-rr**(-6.0d0))
	  
	  end
	  
	  subroutine LJ_potential_PBC(p1, p2, potential)
	  implicit none
	  include 'param.inc'
	  
	  double precision p1(3),p2(3),potential
	  double precision x_ij,y_ij,z_ij
	  double precision rr,rc,uc
	  double precision tot
	  integer k
	  
	  x_ij = (p1(1)-p2(1))
	  y_ij = (p1(2)-p2(2))
	  z_ij = (p1(3)-p2(3))
	  
	  rr = 0.0d0
	  rr = sqrt(x_ij**2.0d0+y_ij**2.0d0+z_ij**2.0d0)
	  rc = sqrt(3.0d0)/2.0d0*box_size
	  uc = 4.0d0*(rc**(-12.0d0)-rc**(-6.0d0))
	  
	  if(rr.gt.rc) then
	    potential = 0.0d0
	  else
	    ! apply PBC 
	  
	    tot = 1.0d0
	  
	    if(x_ij.gt.box_size/2.0d0*tot) then
	      x_ij = x_ij - box_size
	    endif
	  
	    if(x_ij.lt.-box_size/2.0d0*tot) then
	      x_ij = x_ij + box_size
        endif
	  
	    if(y_ij.gt.box_size/2.0d0*tot) then
	      y_ij = y_ij - box_size
	    endif
	  
	    if(y_ij.lt.-box_size/2.0d0*tot) then
	      y_ij = y_ij + box_size
        endif
	  
	    if(z_ij.gt.box_size/2.0d0*tot) then
	      z_ij = z_ij - box_size
	    endif
	  
	    if(z_ij.lt.-box_size/2.0d0*tot) then
	      z_ij = z_ij + box_size
        endif
	  
	    ! end apply PBC
	    rr = 0.0d0
	    rr = sqrt(x_ij**2.0d0+y_ij**2.0d0+z_ij**2.0d0)
	  
	    if(rr.lt.1.0d-6) then
	      write(*,*) 'rr wrong', rr
		  pause
        endif
	    potential = 4.0d0*(rr**(-12.0d0)-rr**(-6.0d0))-uc 
      endif
	  
	  end