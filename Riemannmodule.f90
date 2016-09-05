subroutine Riemann(gamma,PL,PR,rhoL,rhoR,uL,uR,Pstar,ustar,success)

	! Non-relativistic, exact Riemann solver
	! Takes input for right and left states (P,rho,u) respectively
	! Computes pressure flux iteratively using a Newton-Raphson
	! Computes velocity flux
	! Outputs P*,u*

	real 		gamma,PL,PR,rhoL,rhoR,AL,AR,BL,BR,csL,csR,G2
	real		initial_guess,Pstar,ustar
	integer		success

	call factors(gamma,PL,PR,rhoL,rhoR,AL,AR,BL,BR,csL,csR,G2) 
	initial_guess = guess(gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2) 
	call NR(initial_guess,gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2,Pstar,success)
	if (success .eq. 2) then
		call Bisection(gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2)
	endif
	ustar = u_star(Pstar,gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2)

	end


	subroutine factors(gamma,PL,PR,rhoL,rhoR,AL,AR,BL,BR,csL,csR,G2)
	! subroutine returns useful factors to be used in Riemann solver

	!input
	real 		gamma,PL,PR,rhoL,rhoR
	!output
	real		AL,AR,BL,BR,csL,csR,G2

	AL = 2.0d0/((gamma+1.0d0)*rhoL)
	AR = 2.0d0/((gamma+1.0d0)*rhoR)
	BL = PL*(gamma-1.0d0)/(gamma+1.0d0)
	BR = PR*(gamma-1.0d0)/(gamma+1.0d0)
	csL = (gamma*PL/rhoL)**0.5d0
	csR = (gamma*PR/rhoR)**0.5d0
	G2 = (gamma-1.0d0)/(2.0d0*gamma)

	return
	end

    real function getsign(y)
    ! returns the sign of a number
    real		y
    if (y .lt. abs(y)) then
    	getsign = -1.0d0
        return
    else
    	getsign = 1.0d0
        return 
    endif
    return
    end



	real function F(Pguess,gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2)
	! Find 0's of this function

	!input
	real 		gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2,Pguess
	!output
	real 		vR,vL

	if (Pguess.le.PR) then
	    ! Right rarefaction wave
	    vR = (2.0d0*csR/(gamma-1.0))*((Pguess/PR)**G2-1.0)
	else if (Pguess.gt.PR) then
	    ! Right shock wave
	    vR = (Pguess-PR)*(AR/(Pguess+BR))**0.5d0
	endif

	if (Pguess.le.PL) then
	    ! Left rarefaction wave
	    vL = (2.0d0*csL/(gamma-1.0))*((Pguess/PL)**G2-1.0)
	else if (Pguess.gt.PL) then
	    ! Left shock wave
	    vL = (Pguess-PL)*(AL/(Pguess+BL))**0.5d0
	endif

	F = vR+vL+uR-uL

	return
	end



	real function dFdP(Pguess,gamma,PL,PR,rhoL,rhoR,AL,AR,BL,BR,csL,csR,G2)
	! Find derivative of function F wrto Pguess

	!input
	real		gamma,PL,PR,rhoL,rhoR,AL,AR,BL,BR,csL,csR,G2,Pguess
	!output
	real		dvR,dvL

	if (Pguess.le.PR) then
	    !Right rarefaction wave derivative
	    dvR = (1.0d0/(rhoR*csR))*(Pguess/PR)**(-(gamma+1.0d0)/(2.0d0*gamma))
	else if (Pguess.gt.PR) then
	    !Right shock wave derivative
	    dvR = ((AR/(BR+Pguess))**0.5d0)*(1.0d0-(Pguess-PR)/(2.0d0*(BR+Pguess)))
	endif

	if (Pguess.le.PL) then
	    ! Left rarefaction wave derivative
	    dvL = (1.0d0/(rhoL*csL))*(Pguess/PL)**(-(gamma+1.0d0)/(2.0d0*gamma))
	else if (Pguess.gt.PL) then
	    ! Left shock wave derivative
	    dvL = ((AL/(BL+Pguess))**0.5d0)*(1.0d0-(Pguess-PL)/(2.0d0*(BL+Pguess)))
	endif

	dFdP = dvL+dvR

	return
	end



	real function guess(gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2)
	!Returns guess depending on conditions

	!input
	real 		gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2
	!intermediate
	real		Pmin,Pmax,Fmin,Fmax

	if (PL.gt.PR) then
		Pmin = PR
		Pmax = PL
	else if (PR.ge.PL) then
		Pmin = PL
		Pmax = PR
	endif

	Fmin = F(Pmin,gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2)
	Fmax = F(Pmax,gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2)

	if (Fmin.gt.0.0d0 .and. Fmax.gt.0.0d0) then
	    ! Two rarefaction waves
	    guess = ((csL+csR-0.5d0*(gamma-1.0d0)*(uR-uL))/(csL/PL**((gamma-1.0d0)/ &
	        (2.0d0*gamma))+csR/PR**((gamma-1.0d0)/(2.0d0*gamma))))**((2.0d0*gamma)/(gamma-1.0d0))
	else if (Fmin.le.0.0d0 .and. Fmax.gt.0.0d0) then
	    ! One rarefaction, one shock wave
	    if (Pmax.eq.PR) then
	        guess = Pmax
	    else if (Pmax.eq.PL) then
	        guess = Pmin
	    endif
	else if (Fmin.lt.0.0d0 .and. Fmax.lt.0.0d0) then
	    ! Two shock waves
	    guess = 1.5d0*Pmax
	endif

	return
	end


	subroutine NR(Pguess,gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2,Pstar,success)
	! Newton-Raphson algorithm for finding roots of F

	integer		i,n,success
	real		gamma,maxerr,Fi,dFi,Fip1,Pnew,Pguess,Pstar, &
				PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2

	maxerr = 1.0e-6

	n = 20

	do i = 1,n
		Fi = F(Pguess,gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2)
		dFi = dFdP(Pguess,gamma,PL,PR,rhoL,rhoR,AL,AR,BL,BR,csL,csR,G2)
	    Pnew = Pguess-Fi/dFi
	    Pguess = Pnew
	    Fip1 = F(Pguess,gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2)
	    if (abs(Fip1).lt.maxerr) then
	        Pstar = Pnew
	        success = 1
	        return
	    endif
	    if (Pguess .ne. Pguess) then
	    	!print*, "NaN in Pguess, going to bisection..."
	    	Pstar = 0.0d0
	    	success = 2
	    	return 
	    endif
	    if (i == n) then
	        Pstar = 0.0d0
	        success = 0
	        print*, "NR failed"
	        print*, "PL, PR, rhoL, rhoR, uL, uR"
	        print*,  PL, PR, rhoL, rhoR, uL, uR
	        print*, "Pguess = ", Pguess
	        print*, "F = ", Fip1
	        print*,	"STOPPING"
	        stop
	    endif
	enddo

	return
	end



	subroutine Bisection(gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2)
	! If NR fails, brute force a solution

	integer		i,m,success
	real		gamma,maxerr,Pstar, &
				PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2, &
				a,b,c,Fa,Fb,Fc

	maxerr = 1.0e-6
	m = 1e5

	a = 0.0
	b = 1.0e2

	Fa = F(a,gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2)
	Fb = F(b,gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2)

    if (getsign(Fa).ne.getsign(Fb)) then
        do i = 1,m
            c = (a+b)/2.0d0
            Fc = F(c,gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2)
            if (abs(Fc).le. maxerr) then
                Pstar = c
                return
            endif
            if (getsign(Fc).eq.getsign(Fa)) then
                a = c
            else if (getsign(Fc).eq.getsign(Fb)) then
                b = c
            endif
            if (i == m) then
                print*, "solution did not converge"
                print*, "PL, PR, rhoL, rhoR, uL, uR"
        		print*,  PL, PR, rhoL, rhoR, uL, uR
        		print*,	"STOPPING"
                stop
            endif
            Fa = F(a,gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2)
			Fb = F(b,gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2)
        enddo
    else
        print*, "root is not bracketed!"
        print*, "Fa,Fb"
        print*,	Fa,Fb
        print*, "PL, PR, rhoL, rhoR, uL, uR"
        print*,  PL, PR, rhoL, rhoR, uL, uR
        print*,	"STOPPING"
        stop
    endif

	return
	end



	real function u_star(Pguess,gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2)
	! get u_star

	!input
	real 		gamma,PL,PR,rhoL,rhoR,uL,uR,AL,AR,BL,BR,csL,csR,G2,Pguess
	!output
	real 		vR,vL

	if (Pguess.le.PR) then
	    ! Right rarefaction wave
	    vR = (2.0d0*csR/(gamma-1.0))*((Pguess/PR)**G2-1.0)
	else if (Pguess.gt.PR) then
	    ! Right shock wave
	    vR = (Pguess-PR)*(AR/(Pguess+BR))**0.5d0
	endif

	if (Pguess.le.PL) then
	    ! Left rarefaction wave
	    vL = (2.0d0*csL/(gamma-1.0))*((Pguess/PL)**G2-1.0)
	else if (Pguess.gt.PL) then
	    ! Left shock wave
	    vL = (Pguess-PL)*(AL/(Pguess+BL))**0.5d0
	endif

	u_star = 0.5d0*(uR+uL+vR-vL)

	return
	end