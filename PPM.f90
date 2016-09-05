
subroutine PPM(symm,gamma,dt,j,x_array,y_array,len,pressure,density,volume,ap,am,a_int,Astar,x_r,do_A)

    implicit none

    ! 2.6 and 1.12 in CW 1984
    ! note: p, m appended here refer to whether the terms will end up in
    ! apj+0.5 or amj+0.5

    integer         j,len,symm,do_A
    real*16         x_array(len),y_array(len),pressure(len),&
                    density(len),volume(len),Astar(len),x_r(len)
    real*16         dt,Cjm,Ajm,xm,ym,Cjp,Ajp,xp,yp,&
                    aj,ajp1,aL,aLp1,aR,aRp1,deltaaj,a_6j,&
                    deltaajp1,a_6jp1,ap,am,a_int,gamma
    
    Cjm = (gamma*pressure(j+1)*density(j+1))**0.5
    Ajm = 1.0d0
    ym = dt*Cjm*Ajm
    xm = ym/x_array(j+1)

    Cjp = (gamma*pressure(j)*density(j))**0.5
    Ajp = 1.0d0
    yp = dt*Cjp*Ajp
    xp = yp/x_array(j)

    call aRL(j,x_array,y_array,len,aj,ajp1,aL,aLp1,aR,aRp1)

    deltaaj = aR-aL
    a_6j = 6.0d0*(aj-0.5d0*(aR+aL))
    deltaajp1 = aRp1-aLp1
    a_6jp1 = 6.0d0*(ajp1-0.5d0*(aRp1+aLp1))

    ap = aR-0.5d0*xp*(deltaaj-(1-(2.0d0/3.0d0)*xp)*a_6j)

    am = aLp1+0.5d0*xm*(deltaajp1+(1-(2.0d0/3.0d0)*xm)*a_6jp1)
    a_int = 0.5d0*(aR+aLp1)

    return
    end




    real*16 function d2aj(j,x_array,y_array,len)
    !eqn 1.17
    integer     j,len
    real*16     x_array(len),y_array(len)
    real*16     m,mm1,mp1,mp2,aj,ajm1,ajp1,ajp2,&
                A,B,C

    m = x_array(j)
    mm1 = x_array(j-1)
    mp1 = x_array(j+1)
    mp2 = x_array(j+2)

    aj = y_array(j)
    ajm1 = y_array(j-1)
    ajp1 = y_array(j+1)
    ajp2 = y_array(j+2)

    A = (1.0d0/(mm1+m+mp1))
    B = (mp1+m)
    C = (m+mm1)

    d2aj = A*((ajp1-aj)/B - (aj-ajm1)/C)

    return
    end


    real*16 function d_a(j,x_array,y_array,len)
    !This function calculates the variable d_a in CW 1983, eq'n 1.7   
    integer     j,len
    real*16     x_array(len),y_array(len)
    real*16     m,mm1,mp1,mp2,aj,ajm1,ajp1,ajp2,&
                A,B,C,d_aR,d_aL

    m = x_array(j)
    mm1 = x_array(j-1)
    mp1 = x_array(j+1)
    mp2 = x_array(j+2)

    aj = y_array(j)
    ajm1 = y_array(j-1)
    ajp1 = y_array(j+1)
    ajp2 = y_array(j+2)

    ! condensed terms and factors
    ! Note A,B,C are defined in the order these
    ! factors occur in CW 1983 1.7
    ! A,B,C unitless
    A = (m/(mm1+m+mp1))
    B = (2.0d0*mm1+m)/(mp1+m)
    C = (2.0d0*mp1+m)/(mm1+m)
    d_aR = (ajp1-aj)
    d_aL = (aj-ajm1)
    d_a = (A*(B*d_aR+C*d_aL))
    return
    end


    real*16 function dmin_a(j,x_array,y_array,len)
    ! CW 1984 Eq'n 1.8 - slope limiter
    integer     j,len
    real*16     x_array(len),y_array(len)
    real*16     da_array,aj,ajm1,ajp1,ajp2,d_a,getsign

    da_array = d_a(j,x_array,y_array,len)

    aj = y_array(j)
    ajm1 = y_array(j-1)
    ajp1 = y_array(j+1)
    ajp2 = y_array(j+2)

    if ((ajp1-aj)*(aj-ajm1).gt.0.0d0) then
        dmin_a = min(abs(da_array),2*abs(aj-ajm1))*getsign(da_array)
        return
    else
        dmin_a = 0.0d0
        return
    endif
    end



    real*16 function ahalf(j,x_array,y_array,len)
    ! This function calculates values for a(j+0.5), as introduced in CW 1983, eq'n 1.6
    integer     j,len
    real*16     x_array(len),y_array(len)
    real*16     m,mm1,mp1,mp2,aj,ajm1,ajp1,ajp2,&
                A,B,C,D,Ec,F,G,d_aR,d_aL,d_aj,d_ap1,dmin_a

    ! Define m values
    ! m here is Xi in eq'n 1.6 of CW 1983
    m = x_array(j)
    mm1 = x_array(j-1)
    mp1 = x_array(j+1)
    mp2 = x_array(j+2)

    aj = y_array(j)
    ajm1 = y_array(j-1)
    ajp1 = y_array(j+1)
    ajp2 = y_array(j+2)

    ! Note that these A,B,C,D's etc are different
    ! from those specified in d_a(...)
    ! They are defined in order of appearance in CW 1983 1.6
    A = (m/(m+mp1))
    ! B has units 1/[x_array]
    B = (1.0d0/(mm1+m+mp1+mp2))
    ! C has units [x_array]
    C = (2.0d0*mp1*m)/(m+mp1)
    D = (mm1+m)/(2.0d0*m+mp1)
    Ec = (mp2+mp1)/(2.0d0*mp1+m)
    ! F, G have units of [x_array]
    F = m*(mm1+m)/(2.0d0*m+mp1)
    G = mp1*(mp1+mp2)/(m+2.0d0*mp1)

    d_aR = (ajp1-aj)
    d_aL = (aj-ajm1)

    ! Call and evaluate delta_a functions
    d_aj = dmin_a(j,x_array,y_array,len)
    d_ap1 = dmin_a(j+1,x_array,y_array,len)
    ahalf = (aj+A*d_aR+B*(C*(D-Ec)*d_aR-F*d_ap1+G*d_aj))    

    !if (aj.lt.0 .or. aL .lt. 0 .or. aR .lt. 0) then
    !    print*, "interpolation less than 0..."
    !    print*,"j=",j
    !    print*,"aj=",aj
    !    print*,"aL=",aL
    !    print*,"aR=",aR
    !    stop
    !endif

    return
    end


    subroutine aRL(j,x_array,y_array,len,aj,ajp1,aL,aLp1,aR,aRp1)
    ! Here we look at cases and set aR, aL to be used
    ! in the rest of the PPM scheme
    ! aR,j = aj+1/2 --> ahalf(j,...)
    ! aL,j = aj-1/2 --> ahalf(j-1,...)
    integer     j,len
    real*16     y_array(len),x_array(len)
    real*16     aj,aL1,aR1,eta,ajp1,aLp1,aRp1,&
                adL,adR,aL,aR,ahalf

    ! define aRj, aLj (left and right points of j)
    !aj = y_array(j)
    !aL1 = ahalf(j-1,x_array,y_array,len)
    !aR1 = ahalf(j,x_array,y_array,len)

    !call adLR(j,x_array,y_array,len,adL,adR)

    !eta = etaj(j,x_array,y_array,len)

    !aL = aL1*(1-eta)+adL*eta
    !aR = aR1*(1-eta)+adR*eta

    aj = y_array(j)
    aL = ahalf(j-1,x_array,y_array,len)
    aR = ahalf(j,x_array,y_array,len)

    ! check if aj is an extremum
    if (((aR-aj)*(aj-aL)).le.0.0d0) then
        aL = aj
        aR = aj
    else if ((aR-aL)*(aj-0.5*(aL+aR)).gt. &
         (1.0d0/6.0d0)*(aR-aL)**2) then
        aL = 3.0d0*aj-2.0d0*aR
    else if ((-1.0d0/6.0d0)*(aR-aL)**2 .gt. &
         (aR-aL)*(aj-0.5d0*(aR+aL))) then
        aR = 3.0d0*aj-2.0d0*aL
    endif

    ! define aRj+1,aLj+1 (left and right points of j+1)
    ajp1 = y_array(j+1)
    aLp1 = ahalf(j,x_array,y_array,len)
    aRp1 = ahalf(j+1,x_array,y_array,len)

    if ((aRp1-ajp1)*(ajp1-aLp1).le.0.0d0) then
        aLp1 = ajp1
        aRp1 = ajp1
    else if ((aRp1-aLp1)*(ajp1-0.5d0*(aLp1+aRp1)) .gt. &
         (1.0d0/6.0d0)*(aRp1-aLp1)**2) then
        aLp1 = 3.0d0*ajp1-2.0d0*aRp1

    else if ((-1.0d0/6.0d0)*(aRp1-aLp1)**2 .gt. &
         (aRp1-aLp1)*(ajp1-0.5d0*(aRp1+aLp1))) then
        aRp1 = 3.0d0*ajp1-2.0d0*aLp1
    endif
    return
    end

