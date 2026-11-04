MODULE renode_

! Re-nodes each contour while preserving corner locations.
! Uses square-root dependence on non-local curvature

    USE dyn_com
    USE cubic_
    USE dmarray_
    PRIVATE
    PUBLIC renode

CONTAINS

SUBROUTINE renode

    REAL(kr), DIMENSION(:), ALLOCATABLE :: a,b,c,dis,e,g,h,eu,alp,bet,gam
    REAL(kr) :: eta,p,q,fac,dsum
    INTEGER  :: im,lm,l,ia,ib,i1,i2,mbo
    LOGICAL last

    ALLOCATE (a(npt),b(npt),c(npt),dis(npt),e(npt),g(npt),h(npt),   &
              eu(npt),alp(npt),bet(npt),gam(npt))
    
    ! Find corners:
    DO j = 1,n
        DO i = 1,np(j)
            im = mb(j)+i
            ia = mb(j)+MOD(i,np(j))+1
            ib = mb(j)+MOD(i-2+np(j),np(j))+1
            corn(im) = (x(ia)-x(im))*(x(ib)-x(im)) +  &
                       (y(ia)-y(im))*(y(ib)-y(im)) > zero
        END DO
    END DO
    ! Get the updated cubic interpolation coefficients:
    CALL cubic (eu,alp,bet,gam,c,e,g)
    
    ! Determine the fractional number of nodes to be placed
    ! between each pair of nodes using the node density fn:
    h = SQRT((e/ell)**2+c*c)
    ! Use the spherical curvature expression to ensure solely
    ! that the minimum node density for each segment is positive :
    xi = x+half*dx
    yi = y+half*dy
    DO i = 1,npt
        dis = (x(i)-xi)**2 + (y(i)-yi)**2
        a(i) = SUM (h*abs(o)/(dis+dm**2), MASK = dis < rdsq)
        b(i) = SUM (e*abs(o)/(dis+dm**2), MASK = dis < rdsq)
    ! original code simply had g=a/b:replace with following
    IF  ((a(i) .EQ. 0.d0) .AND. (a(i) .EQ. 0.d0)) THEN
    g(i) = 0.d0
    ELSE
    g(i)=a(i)/b(i)
    ENDIF
    END DO

    h(:npt-1) = half*(g(:npt-1)+g(2:))
    h(mb(:n)+np(:n)) = half*(g(mb(:n)+np(:n))+g(mb(:n)+1))
    o = e * MIN (one/dm, densf*SQRT(h)+half*h )

    ! npt tracks the total number of points on all renoded contours:
    npt = 0   

    DO j=1,n
        ! Renode between corners, if they exist:
        ! l tracks the total number of points on a renoded contour
        l=1
        mbo = mb(j)
        mb(j) = npt
        i1 = mbo+1
r1:     DO
            cornw(npt+l) = corn(i1)
            xw(npt+l) = x(i1)
            yw(npt+l) = y(i1)
            dsum = zero
            i = i1
            DO
                dsum = dsum+o(i)
                i = i+1
                last = i > mbo+np(j)
                IF (last) THEN
                    IF (dsum < small) THEN
                        l = l-1
                        EXIT r1
                    END IF
                    i2 = i-1
                    EXIT
                ELSE IF (corn(i)) THEN
                    IF (dsum < small) THEN
                        i1 = i
                        CYCLE r1
                    END IF
                    i2 = i-1
                    EXIT
                END IF
            END DO

            ! m-1 is the number of nodes to be placed between i1 and i2.
            m = NINT(dsum)+1

            fac = m/dsum
            o(i1:i2) = fac*o(i1:i2)

            ! The sum of the o's is equal to m.
            ! The first node along a contour (segment) is fixed.
            ! Find the adjusted nodal positions:
            q = zero
            i = i1-1
            DO lm = npt+l+1,npt+l+m-1
                IF (q < one) THEN
                    DO
                        i = i+1
                        q = q+o(i)
                        IF (q >= one) EXIT
                    END DO
                END IF
                q = q-one
                p = one-q/o(i)
                eta = p*(alp(i)+p*(bet(i)+p*gam(i)))
                cornw(lm) = .FALSE.
                xw(lm) = x(i)+p*dx(i)-eta*dy(i)
                yw(lm) = y(i)+p*dy(i)+eta*dx(i)
            END DO
            IF (last) THEN
                l = l+m-1
                EXIT r1
            ELSE
                l = l+m
                i1 = i2+1
            END IF
        END DO r1
        np(j) = l
        npt = npt+l
    END DO

    ! Dimension and switch arrays around again:
    CALL dimarrays

    DEALLOCATE (a,b,c,dis,e,g,h,eu,alp,bet,gam)

END SUBROUTINE renode

END MODULE renode_

