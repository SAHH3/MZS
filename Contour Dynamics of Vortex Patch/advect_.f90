MODULE advect_

    USE dyn_com
    USE uv_
    IMPLICIT NONE
    PRIVATE
    PUBLIC advect

CONTAINS

! Advects all of the nodes by a 4th order runga-kutta scheme.

SUBROUTINE advect

    REAL(kr) :: dt2,dt3,dt6
    dt2 = half*dt
    dt3 = dt/three
    dt6 = half *dt3
    xi = x
    yi = y
    xf = x
    yf = y
    
    CALL uv
    x = xi+dt2*u
    y = yi+dt2*v
    xf = xf+dt6*u
    yf = yf+dt6*v
    t = t+dt2
    CALL uv
    x = xi+dt2*u
    y = yi+dt2*v
    xf = xf+dt3*u
    yf = yf+dt3*v
    CALL uv
    x = xi+dt*u
    y = yi+dt*v
    xf = xf+dt3*u
    yf = yf+dt3*v
    t = t+dt2
    CALL uv
    x = xf+dt6*u
    y = yf+dt6*v
    
        do i=1,npt 
     if (x(i)**2+y(i)**2 .lt. rad1**2) then 
        x(i)=rad1*x(i)/sqrt(x(i)**2+y(i)**2) 
        y(i)=rad1*y(i)/sqrt(x(i)**2+y(i)**2) 
     elseif ((x(i)-cent)**2+y(i)**2 .lt. rad2**2) then  
        x(i)=rad2*(x(i)-cent)/sqrt((x(i)-cent)**2+y(i)**2)+cent 
        y(i)=rad2*y(i)/sqrt((x(i)-cent)**2+y(i)**2) 
     endif 
    end do   

END SUBROUTINE advect

END MODULE advect_
