
program main

implicit none

real*16 		gamma,PL,PR,rhoL,rhoR,uL,uR,grid_length,tf
integer			success,ncells,nsteps,alpha,BC,k,j,preset,t_out,do_ppm
parameter		(ncells = 1997)
parameter		(nsteps = 100)
real*16 		P(nsteps,ncells),rho(nsteps,ncells), &
				u(nsteps,ncells),x_c(nsteps,ncells), &
				vol(nsteps,ncells),xc_i(ncells),rho_i(ncells), &
				E(nsteps,ncells)

PL = 1.0d0
PR = 1.0d0
rhoL = 1.0d0
rhoR = 1.0d0
uL = 0.0d0
uR = 0.0d0
gamma = 1.33d0
alpha = 2
grid_length = 1.0d0
BC = 34
tf = 6000.0
preset = 1
do_ppm = 0


open(unit = 1, file = "C:\Users\Carly\Desktop\Fortran_out\Hydro_out.txt",action="write",status="replace")

call Hydro1d(do_ppm,preset,ncells,nsteps,alpha,gamma,grid_length,&
					BC,PL,PR,rhoL,rhoR,uL,uR,tf,t_out,P,rho,u,x_c,vol,E,xc_i,rho_i)

do k = 1,ncells
	write(1,*) x_c(t_out,k), P(t_out,k), rho(t_out,k), u(t_out,k), vol(t_out,k), E(t_out,k), xc_i(k), rho_i(k)
enddo

close(1)
end