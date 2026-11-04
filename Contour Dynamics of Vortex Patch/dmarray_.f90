MODULE dmarray_
    USE dyn_com
    IMPLICIT NONE
    PRIVATE
    PUBLIC dimarrays

CONTAINS

!   Allocates the variable length arrays

SUBROUTINE dimarrays

    DEALLOCATE (corn)
    DEALLOCATE (x,y,o,u,v,dx,dy)
    DEALLOCATE (xi,yi,xf,yf)
    deallocate(ugap,vgap,xzi,eta)
    deallocate(z,zeta,ucomp)
    
!    DEALLOCATE (o,eu,alp,bet,gam)

    ALLOCATE (corn(npt))
    ALLOCATE (x(npt),y(npt),o(npt),u(npt),v(npt),dx(npt),dy(npt))
    ALLOCATE (xi(npt),yi(npt),xf(npt),yf(npt))
    allocate (ugap(npt),vgap(npt),xzi(npt),eta(npt))
    allocate(z(npt),zeta(npt),ucomp(npt))   
!    ALLOCATE (o(npt),eu(npt),alp(npt),bet(npt),gam(npt))

    x = xw(:npt)
    y = yw(:npt)
    corn = cornw(:npt)

    DO j = 1,n
        o(mb(j)+1:mb(j)+np(j)) = om(j)
    END DO

END SUBROUTINE dimarrays

END MODULE dmarray_
