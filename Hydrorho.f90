subroutine Hydro1d(do_ppm,preset,ncells,nsteps,alpha,gamma,grid_length,&
					BC,PL,PR,rhoL,rhoR,uL,uR,tf,t_out,P,rho,u,x_c,vol,E, &
                    xc_i,rho_i)

    implicit none

    ! authors/papers to check out...
    ! omer bomberg
    ! couch & wheeler
    ! sari/nakar
    ! drout, chornock
    ! totheus SN 2009 ip
    ! 2011dh 

    ! TO DO: arbitrary density profile, shock speed vs. position
    ! Sedov blast
    ! shockspeed: try recording max_i(-div(v)), record with zone (get x_c(i,shock_loc))

    ! MOCK-UP:
        
    !               !==============!==============!==============!==============!
    !   Ghost L     ! (Pj,uj,rho)  !              !              !              !  Ghost R
    !               !      x       !      x       !     ...      !      x       !
    !  uj*,Pj*,Aj*  !    (j=0)     !    (j=1)     !              ! (j=ncells-1) ! uj*,Pj*,Aj*
    !    j = 0      !==============!==============!==============!==============!  j = ncells
    !                         uj*,Pj*,Aj*        ...         uj*,Pj*,Aj*
    !                            j = 1                       j = ncells-1       


    !==================================================================
    !  1D NON-RELATIVISTIC LAGRANGIAN HYDROCODE + EXACT RIEMANN SOLVER
    !==================================================================
    ! COMPILE WITH "RIEMANNMODULE.F90" USING COMMAND 
    ! "GFORTRAN CALLHYDRO.F90 HYDRO2.F90 RIEMANNMODULE.F90 PPM.F90 -O (OUTPUT PROGRAM NAME)"


    ! FEB 26: 0-INDEXING ADDED, PPM ADDED

    ! INPUTS
    integer		nsteps,ncells,alpha,BC,preset,success,do_ppm
    real*16		gamma,grid_length,PL,PR,rhoL,rhoR,uL,uR,tf,Pstar,ustar

    ! CONSTANTS
    real*16		CFL,gamm1,gamp1,gamfc,gamfac,m2

    ! SIMULATION VALUES
    integer		i,j,t_out,shocki
    real*16		cell_length,boundary,dt,tstamp

    real*16		dtdm(0:ncells-1),dm(0:ncells-1), &
    			rho(0:nsteps-1,0:ncells-1),P(0:nsteps-1,0:ncells-1), &
                u(0:nsteps-1,0:ncells-1),vol(0:nsteps-1,0:ncells-1), &
                E(0:nsteps-1,0:ncells-1),e_int(0:nsteps-1,0:ncells-1), &
    			cs(0:nsteps-1,0:ncells-1),S(0:nsteps-1,0:ncells-1), &
                g(0:nsteps-1,0:ncells-1),A(0:nsteps-1,0:ncells-1), &
                du(0:ncells-1),dxdcs(0:ncells-1),x_l(0:nsteps-1,0:ncells-1), &
                x_c(0:nsteps-1,0:ncells-1),x_r(0:nsteps-1,0:ncells-1),&
                shockloc(0:nsteps-1)

    real*16		P_star(0:ncells),u_star(0:ncells),vstar(0:ncells), &
                Astar(0:ncells),ppm_PL(0:ncells-1),ppm_PR(0:ncells-1), &
                ppm_rhoL(0:ncells-1),ppm_rhoR(0:ncells-1),ppm_uL(0:ncells-1), &
                ppm_uR(0:ncells-1),P_interp(0:ncells-1),rho_interp(0:ncells-1), &
                duL(0:ncells-1),duR(0:ncells-1),dx(0:ncells-1),cellsqueeze(0:ncells-1),&
                xc_i(0:ncells-1),rho_i(0:ncells-1),divlist(0:ncells-1),tlist(0:ncells-1), &
                xlist(0:ncells-1),r(0:ncells+2),rhostar(0:ncells+2),shockv(0:ncells-2)

    real*16		PLfix,PRfix,rhoLfix,rhoRfix,uLfix,uRfix,ELfix,ERfix, &
    			PLghost,PRghost,rhoLghost,rhoRghost,uLghost,uRghost, &
    			ELghost,ERghost,null,div_v,ap,am,a_int,vweight,getindex, &
                maxsqueeze,M,Gc

    !======================================================================

    open(unit = 3, file = "C:\Users\Carly\Desktop\Shock_out.txt",action="write",status="replace")
    open(unit = 4, file = "C:\Users\Carly\Desktop\BSG_rrho_reformat.txt",action="read")
    open(unit = 5, file = "C:\Users\Carly\Desktop\BSG_rrho_Lagrangian.txt",action="write",status="replace")
    open(unit = 2, file = "C:\Users\Carly\Desktop\Shock2.txt",action="write",status="replace")

    do j = 0,ncells+2
        read(4,*) r(j), rhostar(j)
    enddo
    close(4)


    ! look at removing these, and constants from Riemann solver
    ! (n, maxerr) and placing them all in a common file...

    CFL = 0.8d0
    vweight = 1.0d0 ! weights velocity in the CFL condition
    gamm1 = gamma-1.0d0
    gamp1 = gamma+1.0d0
    gamfc = (gamp1/2.0d0)
    gamfac = (gamp1/(2.0d0*gamma))
    m2 = (gamma-1.0d0)/(gamma+1.0d0)

    !====== SET-UP =====
    cell_length = grid_length/real(ncells)
    boundary = 0.5d0*grid_length

    do j = 0,ncells+2
        if (j .ge. 1 .and. j .le. (ncells)) then
            x_c(0,j-1) = r(j)
            rho(0,j-1) = rhostar(j)
        endif
    enddo

    do j = 0,ncells+2
        if (j.ge.1 .and. j.le.(ncells)) then
            x_l(0,j-1) = 0.5d0*(r(j)+r(j-1)) 
        endif 
    enddo

    do j = 0,ncells-1
        if (j == (ncells-1)) then
            x_r(0,j) = 0.5d0*(r(j+2)+r(j+1))
        else
            x_r(0,j) = x_l(0,j+1)
        endif
        ! print*, x_l(0,j),x_r(0,j)
    enddo

    do j = 0,ncells-1
        vol(0,j) = ((x_r(0,j)**(alpha+1.0d0))-(x_l(0,j)**(alpha+1.0d0)))/( &
                                    alpha+1.0)
        !print*, rho(0,j)
        dm(j) = rho(0,j)*vol(0,j) 
    enddo

    !M = sum(dm(:)) ! total mass
    !Gc = 6.674e-8     ! Newton's constant, cgs

    if (preset == 0) then
    	! SOD SHOCK TUBE
    	do j = 0,ncells-1
    	    if (j/real(ncells).le.boundary) then
    	        rho(0,j) = rhoL
    	        P(0,j) = PL
    	        u(0,j) = uL
    	    else
    	        rho(0,j) = rhoR
    	        P(0,j) = PR
    	        u(0,j) = uR
    	    endif
    	enddo
    else if (preset == 1) then
    	!	SEDOV BLAST
    	do j = 0,ncells-1
            if (j.le.2) then
                P(0,j) = 1.0d51
            else
    	       P(0,j) = 1.0d-50
            endif
    	    u(0,j) = 0.0d0
    	enddo
    endif

    xc_i(:) = x_c(0,:)
    rho_i(:) = rho(0,:)


    !======= GENERATE INITIAL CONDITIONS =====
    do j = 0,ncells-1
        A(0,j) = ((x_r(0,j)**(alpha+1.0d0))-(x_l(0,j)**(alpha+1.0d0)))/( &
                                (alpha+1.0d0)*(x_r(0,j)-x_l(0,j))) 
        cs(0,j) = (gamma*P(0,j)/rho(0,j))**0.5
        if (rho(0,j).ne.0.0d0) then
            E(0,j) = P(0,j)/((gamma-1.0d0)*rho(0,j))+0.5d0*u(0,j)**2
        else
            E(0,j) = 0.0d0
        endif
        S(0,j) = log(P(0,j)/(rho(0,j)**gamma)) 
        e_int(0,j) = E(0,j)-0.5d0*u(0,j)**2
        g(0,j) = 0.0d0
        divlist(j) = 0.0d0
        tlist(j) = 0.0d0
        xlist(j) = 0.0d0
    enddo


    !======= HANDLE FIXED BOUNDARY CONDITIONS HERE =====
    PLfix = P(0,0)
    PRfix = P(0,ncells-1)
    rhoLfix = rho(0,0)
    rhoRfix = rho(0,ncells-1)
    uLfix = u(0,0)
    uRfix = u(0,ncells-1)
    ELfix = E(0,0)
    ERfix = E(0,ncells-1)

    !================ MAIN LOOP STARTS HERE =====================

    tstamp = 0
    i = 0
    do while (tstamp/tf.le.1.0d0)
        print*, 100.0d0*tstamp/tf, "% done, nzones,i = ", ncells

        do j = 0,ncells-1
        	if (vol(i,j).lt.0.0d0) then
            	print*, "negative volume at ",j," -- STOPPING..."
            	stop
        	endif
        enddo

        if (tstamp.lt.0.0d0) then
            print*, "negative time step -- STOPPING..."
            stop
        endif

        ! Refresh if time array is full
        if (i.ne.0 .and. mod(i,(nsteps-1)).eq.0) then
            P(0,:) = P(i,:)
            rho(0,:) = rho(i,:)
            u(0,:) = u(i,:)
            x_l(0,:) = x_l(i,:)
            x_c(0,:) = x_c(i,:)
            x_r(0,:) = x_r(i,:)
            vol(0,:) = vol(i,:)
            A(0,:) = A(i,:)
            cs(0,:) = cs(i,:)
            E(0,:) = E(i,:) 
            S(0,:) = S(i,:) 
            e_int(0,:) = e_int(i,:)
            shockloc(0) = shockloc(i)
            i = 0
        endif

        do j = 0,ncells-1
            dx(j) = x_r(i,j)-x_l(i,j)
            if (j.gt.0 .and. j.lt.(ncells-1)) then
                duL(j) = abs(u(i,j)-u(i,j-1))
                duR(j) = abs(u(i,j)-u(i,j+1))
            else
                duL(j) = 0.0d0
                duR(j) = 0.0d0
            endif
        enddo

        dt = CFL*minval(dx(:)/(vweight*(duL(:)+duR(:))+cs(i,:)))
        tstamp = tstamp+dt

        ! EVALUATE GHOST CELLS HERE
        if (BC == 55) then
            ! left, right fixed
            PLghost = PLfix
            PRghost = PRfix
            rhoLghost = rhoLfix
            rhoRghost = rhoRfix
            uLghost = uLfix
            uRghost = uRfix
            ELghost = ELfix
            ERghost = ERfix
        else if (BC == 11) then
            ! Periodic
            PLghost = P(i,ncells-1)
            PRghost = P(i,0)
            rhoLghost = rho(i,ncells-1)
            rhoRghost = rho(i,0)
            uLghost = u(i,ncells-1)
            uRghost = u(i,0)
            ELghost = E(i,ncells-1)
            ERghost = E(i,0)
        else if (BC == 33) then
            ! left, right reflective
            PLghost = P(i,0)
            PRghost = P(i,ncells-1)
            rhoLghost = rho(i,0)
            rhoRghost = rho(i,ncells-1)
            uLghost = -u(i,0)
            uRghost = -u(i,ncells-1)
            ELghost = E(i,0)
            ERghost = E(i,ncells-1)
        else if (BC == 03) then
            ! left fixed, right reflective
            PLghost = PLfix
            PRghost = P(i,ncells-1)
            rhoLghost = rhoLfix
            rhoRghost = rho(i,ncells-1)
            uLghost = uLfix
            uRghost = -u(i,ncells-1)
            ELghost = ELfix
            ERghost = E(i,ncells-1)
        else if (BC == 30) then
            ! left reflective, right fixed
            PLghost = P(i,0)
            PRghost = PRfix
            rhoLghost = rho(i,0)
            rhoRghost = rhoRfix
            uLghost = -u(i,0)
            uRghost = uRfix
            ELghost = E(i,0)
            ERghost = ERfix
        else if (BC == 34) then
            ! left reflective, right vaccuum
            PLghost = P(i,0)
            PRghost = 0.0d0
            rhoLghost = rho(i,0)
            rhoRghost = 0.0d0
            uLghost = -u(i,0)
            uRghost = 0.0d0
            ELghost = E(i,0)
            ERghost = 0.0d0
        endif

        do j = 0,ncells-1
            if (P(i,j).lt.0.0d0) then
                print*, "negative pressure at coord j=",j
                stop
            endif
            if (rho(i,j).lt.0.0d0) then
                print*, "negative density at coord j=",j
                stop
            endif
        enddo

       if (do_ppm == 1) then
            do j = 0,ncells-1
                vol(i,j)=((x_r(i,j)**(alpha+1.0d0))-(x_l(i,j)**(alpha+1.0d0)))/( &
                                    alpha+1.0d0)
                A(i,j) = ((x_l(i+1,j)**(alpha+1.0d0))-(x_l(i,j)**(alpha+1.0d0)))/( &
                                (alpha+1.0d0)*(x_l(i+1,j)-x_l(i,j)))   
                
                if ((j .le. 1) .or. (j.ge.(ncells-5))) then
                    ppm_PL(j) = P(i,j)
                    ppm_PR(j) = P(i,j)
                    ppm_uL(j) = u(i,j)
                    ppm_uR(j) = u(i,j)
                    ppm_rhoL(j) = rho(i,j)
                    ppm_rhoR(j) = rho(i,j)
                else
                    call PPM(alpha,gamma,dt,j,dm,P(i,:),ncells,P(i,:),rho(i,:), &
                            vol(i,:),ap,am,a_int)
                    ppm_PL(j) = ap
                    ppm_PR(j) = am
                    P_interp(j) = a_int
                    call PPM(alpha,gamma,dt,j,dm,rho(i,:),ncells,P(i,:),rho(i,:), &
                            vol(i,:),ap,am,a_int)
                    ppm_rhoL(j) = ap
                    ppm_rhoR(j) = am
                    rho_interp(j) = a_int
                    call PPM(alpha,gamma,dt,j,dm,u(i,:),ncells,P(i,:),rho(i,:), &
                            vol(i,:),ap,am,a_int)
                    ppm_uL(j) = ap
                    ppm_uR(j) = am
                    null = a_int

                endif
            enddo

            do j = 0,ncells-1
                if (BC == 55 .or. BC == 33 .or. BC == 03 .or. BC == 30 .or. BC == 34) then
                    if ((j .gt. 0) .and. (j .lt. ncells-5) ) then
                        call Riemann(gamma,ppm_PL(j),ppm_PR(j),ppm_rhoL(j),ppm_rhoR(j), &
                                ppm_uL(j),ppm_uR(j),Pstar,ustar,success)
                        P_star(j) = Pstar
                        u_star(j) = ustar
                        !div_v = (u(i,j-2)-u(i,j))/(x_l(i,j-1)-x_l(i,j))
                        !if (div_v .lt. 0) then
                        !    print*,"postshock at j=",j
                        !    call Riemann(gamma,P(i,j-1),P(i,j),rho(i,j-1),rho(i,j),u(i,j-1),u(i,j),Pstar,ustar,success)
                        !    P_star(j) = Pstar
                        !    u_star(j) = ustar
                        !    print*,"ppm off"
                        !    !P(i,j) = P_interp(j)
                        !    !rho(i,j) = rho_interp(j)
                        !endif
                    else if (j == 0) then
                        call Riemann(gamma,PLghost,P(i,j),rhoLghost,rho(i,j),uLghost,u(i,j),Pstar,ustar,success)
                        P_star(j) = Pstar
                        u_star(j) = ustar
                    else if ((j .ge. (ncells-5)) .and. (j .lt. (ncells-1))) then
                        call Riemann(gamma,P(i,j-1),P(i,j),rho(i,j-1),rho(i,j),u(i,j-1),u(i,j),Pstar,ustar,success)
                        P_star(j) = Pstar
                        u_star(j) = ustar
                    else if (j == (ncells-1)) then
                        call Riemann(gamma,P(i,j),PRghost,rho(i,j),rhoRghost,u(i,j),uRghost,Pstar,ustar,success)
                        P_star(j) = Pstar
                        u_star(j) = ustar
                    endif

                else if (BC == 11) then
                    if (j == 0 .or. j .ge. ncells-5) then
                        call Riemann(gamma,P(i,j-1),P(i,j),rho(i,j-1),rho(i,j),u(i,j-1),u(i,j),Pstar,ustar,success)
                        P_star(j) = Pstar
                        u_star(j) = ustar
                    else
                        call Riemann(gamma,ppm_PL(j-1),ppm_PR(j-1),ppm_rhoL(j-1),ppm_rhoR(j-1), &
                                ppm_uL(j-1),ppm_uR(j-1),Pstar,ustar,success)
                        P_star(j) = Pstar
                        u_star(j) = ustar
                    endif
                endif
            enddo
        endif

        do j = 0,ncells-1
            vol(i,j)=((x_r(i,j)**(alpha+1.0d0))-(x_l(i,j)**(alpha+1.0d0)))/( &
                    alpha+1.0d0)
            A(i,j) = ((x_l(i+1,j)**(alpha+1.0d0))-(x_l(i,j)**(alpha+1.0d0)))/( &
                        (alpha+1.0d0)*(x_l(i+1,j)-x_l(i,j))) 
        enddo

        if (do_ppm == 0) then 
            if (BC == 55 .or. BC == 33 .or. BC == 03 .or. BC == 30 .or. BC == 34) then
                do j = 0,ncells
                    if (j > 0 .and. j < ncells) then
                    	call Riemann(gamma,P(i,j-1),P(i,j),rho(i,j-1),rho(i,j),u(i,j-1),u(i,j),P_star(j),u_star(j),success)
                    else if (j == 0) then
                        call Riemann(gamma,PLghost,P(i,j),rhoLghost,rho(i,j),uLghost,u(i,j),P_star(j),u_star(j),success)
                    else if (j == ncells) then
                        call Riemann(gamma,P(i,j-1),PRghost,rho(i,j-1),rhoRghost,u(i,j-1),uRghost,P_star(j),u_star(j),success)
                    endif
                enddo


            else if (BC == 11) then
                do j = 0,ncells
                    if (j > 0 .and. j < ncells) then
                        call Riemann(gamma,P(i,j-1),P(i,j),rho(i,j-1),rho(i,j),u(i,j-1),u(i,j),P_star(j),u_star(j),success)
                    else if (j == 0) then
                        call Riemann(gamma,PLghost,P(i,j),rhoLghost,rho(i,j),uLghost,u(i,j),P_star(j),u_star(j),success)
                    else if (j == ncells) then
                        call Riemann(gamma,P(i,j-1),PRghost,rho(i,j-1),rhoRghost,u(i,j-1),uRghost,P_star(j),u_star(j),success)
                    endif
                enddo
            endif
        endif

        do j = 0,ncells
            if (j < ncells) then
                x_l(i+1,j) = x_l(i,j)+u_star(j)*dt      
                if (u_star(j).ne.0.0d0) then
                    Astar(j) = ((x_l(i+1,j)**(alpha+1.0d0))-(x_l(i,j)**(alpha+1.0d0)))/( &
                                    (alpha+1.0d0)*u_star(j)*dt)
                else if (u_star(j).eq.0.0d0) then
                    Astar(j) = 0.0d0
                endif
            else 
                x_r(i+1,j-1) = x_r(i,j-1)+u_star(j)*dt 
                if (u_star(j).ne.0.0d0) then
                    Astar(j) = ((x_r(i+1,j-1)**(alpha+1.0d0))-(x_r(i,j-1)**(alpha+1.0d0)))/( &
                                    (alpha+1.0d0)*u_star(j)*dt)
                else if (u_star(j).eq.0.0d0) then
                    Astar(j) = 0.0d0
                endif
            endif
        enddo

        do j = 0,ncells-1
            x_r(i+1,j) = x_l(i+1,j+1)
            if (j == (ncells-1)) then
                x_r(i+1,j) = x_r(i,j) + u_star(j+1)*dt
            endif
        enddo

        do j = 0,ncells
            vstar(j) = u_star(j)*P_star(j)
        enddo

        do j = 0,ncells-1
            x_c(i+1,j) = 0.5d0*(x_r(i+1,j)+x_l(i+1,j))
            vol(i+1,j) = ((x_r(i+1,j)**(alpha+1.0d0))-(x_l(i+1,j)**(alpha+1.0d0)))/( &
                                    alpha+1.0d0)
        enddo
        
        do j = 0,ncells-1
            if (dm(j) .ne. 0.0d0) then
                dtdm(j) = dt/dm(j)
                rho(i+1,j) = dm(j)/vol(i+1,j)
                
                if (BC == 55 .or. BC == 33 .or. BC == 03 .or. BC == 30 .or. BC == 34) then
                    u(i+1,j) = u(i,j)-0.5d0*(Astar(j+1)+Astar(j))*dtdm(j)* &
                                        (P_star(j+1)-P_star(j))+dt*g(i,j)
                    E(i+1,j) = E(i,j)-dtdm(j)*(Astar(j+1)*vstar(j+1)- &
                                    Astar(j)*vstar(j))+0.5d0*dt*(u(i,j)+u(i+1,j))*g(i,j)
                    !print*, u(i,j)

                else if (BC == 11) then
                    ! Periodic (wrap-around) BCs
                    if (j == ncells-1) then
                        u(i+1,j) = u(i,j)-0.5d0*(Astar(1)+Astar(j))*dtdm(j)* &
                                    (P_star(1)-P_star(j))+dt*g(i,j)
                        E(i+1,j) = E(i,j)-dtdm(j)*(Astar(1)*vstar(1)- &
                                    Astar(j)*vstar(j))+0.5d0*dt*(u(i,j)+u(i+1,j))*g(i,j)
                    else if (j < ncells-1) then
                        u(i+1,j) = u(i,j)-0.5d0*(Astar(j)+Astar(j-1))*dtdm(j)* &
                                        (P_star(j+1)-P_star(j))+dt*g(i,j)
                        E(i+1,j) = E(i,j)-dtdm(j)*(Astar(j+1)*vstar(j+1)- &
                                    Astar(j)*vstar(j))+0.5d0*dt*(u(i,j)+u(i+1,j))*g(i,j)
                    endif
                endif
            else
                u(i+1,j) = u(i,j)
                E(i+1,j) = E(i,j)
                rho(i+1,j) = dm(j)/vol(i+1,j)
            endif

            if (E(i+1,j).ne.E(i+1,j)) then
                print*, "energy nan or inf at j=",j
                print*,"E=",E(i+1,j)
                print*,"P*=",P_star(j)
                print*,"u*=",u_star(j)
                print*,"P*j+1=",P_star(j+1)
                print*,"u*j+1=",u_star(j+1)
                print*,"P=",P(i,j)
                print*,"rho=",rho(i,j)
                stop
            else if (u(i+1,j).ne.u(i+1,j)) then
                print*, "velocity nan or inf at j=",j
                stop
            endif
        enddo

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !=====     GET SHOCK POSITION HERE...                                  =====!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do j = 1,ncells-2
            if (-((u(i,j+1)-u(i,j-1))/(x_c(i,j+1)-x_c(i,j-1))) .gt. divlist(j)) then
                divlist(j) = -((u(i,j+1)-u(i,j-1))/(x_c(i,j+1)-x_c(i,j-1)))
                xlist(j) = x_c(i,j)
                tlist(j) = tstamp
                !print*, "new max at j=",j
                !print*, "velocities ", u(i,j+1),u(i,j-1)
            endif
        enddo


        do j = 1,ncells-2
            cellsqueeze(j) = -(u(i,j+1)-u(i,j-1))/(x_c(i,j+1)-x_c(i,j-1))
        enddo
        maxsqueeze = maxval(cellsqueeze(:))
        shocki = getindex(maxsqueeze,cellsqueeze,ncells)
        shockloc(i) = x_c(i,shocki)
        write(2,*) dt, shockloc(i)

        do j = 0,ncells-1
            if (dm(j) .ne. 0.0d0) then
                e_int(i+1,j) = (E(i+1,j)-0.5*u(i+1,j)**2)
            else
                e_int(i+1,j) = 0.0d0
            endif
            if (e_int(i+1,j).lt.0.0d0) then
                print*, "energy <0 at j=",j
                stop
            endif
            P(i+1,j) = gamm1*rho(i+1,j)*e_int(i+1,j)
            S(i+1,j) = log(P(i+1,j)/rho(i+1,j)**gamma)
            cs(i+1,j) = (gamma*P(i+1,j)/rho(i+1,j))**0.5
            g(i+1,j) = 0.0d0
        enddo


        i = i+1

    enddo

    t_out = i

    do j = 0,ncells-1
        write(3,*) tlist(j), xlist(j), divlist(j)
    enddo

    do j = 0,ncells-2
        shockv(j) = (xlist(j+1)-xlist(j))/(tlist(j+1)-tlist(j))
        write(5,*) xc_i(j),rho_i(j),xlist(j),shockv(j)
    enddo

    return
    close(2)
    close(3)
    close(5)
    end

    !======================================================================

    real*16 function getindex(value,array,n)
        ! This function returns the index at which "value" occurs in "array"
        integer     n,k
        real*16     value,array(n)
        
        do k = 0,n-1
            if ((abs(array(k)-value) .lt. 1.0e-20)) then
                getindex = k
            endif
        enddo
        return 
        end function getindex

    !=======================================================================