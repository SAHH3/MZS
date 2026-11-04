MODULE cubic_

!   Calculates the interpolation coefficients for every node
!   on every contour.  The interpolation approximately enforces
!   continuity of curvature (except at corners which have
!   effectively infinite curvature).

    USE dyn_com
    IMPLICIT NONE
    PRIVATE
    PUBLIC cubic

CONTAINS

SUBROUTINE cubic (eu,alp,bet,gam,c,e,g)

    REAL(kr), DIMENSION(:)              :: c,e,g,eu,alp,bet,gam
    REAL(kr), DIMENSION(:), ALLOCATABLE :: a,b,d,f

    ALLOCATE (a(npt),b(npt),d(npt),f(npt))

    dx(:npt-1) = x(2:)-x(:npt-1)
    dy(:npt-1) = y(2:)-y(:npt-1)
    dx(mb(:n)+np(:n)) = x(mb(:n)+1)-x(mb(:n)+np(:n))
    dy(mb(:n)+np(:n)) = y(mb(:n)+1)-y(mb(:n)+np(:n))
    b = dx*dx+dy*dy
    e = SQRT(b)
    eu = one/(b+small)
    a(2:) = b(:npt-1)
    f(2:) = dx(:npt-1)
    g(2:) = dy(:npt-1)
    a(mb(:n)+1) = b(mb(:n)+np(:n))
    f(mb(:n)+1) = dx(mb(:n)+np(:n))
    g(mb(:n)+1) = dy(mb(:n)+np(:n))
    c = f*b+dx*a
    d = g*b+dy*a
    b = (f*dy-dx*g)/SQRT(c**2+d**2+small**3)
    
    ! Set curvature to zero at corners
    WHERE (corn) b = zero
    c(:npt-1) = e(:npt-1)*(b(2:)+b(:npt-1))
    d(:npt-1) = e(:npt-1)*(b(2:)-b(:npt-1))
    c(mb(:n)+np(:n)) = e(mb(:n)+np(:n))*(b(mb(:n)+1)+b(mb(:n)+np(:n)))
    d(mb(:n)+np(:n)) = e(mb(:n)+np(:n))*(b(mb(:n)+1)-b(mb(:n)+np(:n)))
    
    ! Cubic interpolation coefficients:
    gam = d/three
    bet = half*(c-d)
    alp = half*(gam-c)

    DEALLOCATE (a,b,d,f)

END SUBROUTINE cubic

END MODULE cubic_
