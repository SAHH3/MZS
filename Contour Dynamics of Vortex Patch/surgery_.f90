MODULE surgery_

    USE dyn_com
    USE dmarray_
    IMPLICIT NONE
    PRIVATE
    PUBLIC surgery

CONTAINS

SUBROUTINE surgery

    ! Vortex splitting and merger routine

    REAL(kr), DIMENSION(:), ALLOCATABLE :: b,e
    REAL(kr) :: om1,om2,dis,p,ff,gg
    INTEGER  :: ind1,ind2,np1,im,ih,j1,j2,is,js2,js1,js,np2,ii
    LOGICAL  :: keep1,keep2

    ALLOCATE (b(npt),e(npt))

    xw(:npt) = x
    yw(:npt) = y

    ! First the case when a single contour is involved.
    ! Here the contour could break into two, either in the
    ! form of two separate vortices or an annular vortex. 
    ! Both cases are treated similarly.

    j = 0
s1: DO
        j = j+1

        ! Contour "j"

        IF (j > n) EXIT
        CALL prop(j,j1,j2)
       
        DO i = 1,j2-1

            DO is = i+1,j2

            ! Check whether node i gets within a distance dm of the
            ! line segment between is and is+1.
            ff = xw(j1+i)-xw(j1+is)
            gg = yw(j1+i)-yw(j1+is)
            p = (ff*dx(is)+gg*dy(is))/b(is)
            IF ( p*(one-p) < small) CYCLE

            ! If so, a perpendicular can be dropped onto the
            ! line segment from is to is+1.  
            ! Check the distance to the line segment:
            dis = ABS((ff*dy(is)-gg*dx(is)))/e(is)
            IF (dis >= dm) CYCLE

            ! Break the contour into two at node i.  
            ! Temporarily store the two pieces in the (xi,yi) and (xf,yf) arrays.
            ! Get rid of contour j, and then introduce the two new contours.
            ! First piece:
            xw(j1+i) = xw(j1+i)+half*(-ff+p*dx(is))
            yw(j1+i) = yw(j1+i)+half*(-gg+p*dy(is))
            np1 = is-i+1
            keep1 = np1 > 4
            IF (keep1) THEN
                om1 = om(j)
                ind1 = ind(j)
                im = j1+i
                xi(:np1) = xw(im:im+np1-1)
                yi(:np1) = yw(im:im+np1-1)
            END IF

            ! Second piece:
            np2 = j2-is+i
            keep2 = np2 > 4
            IF (keep2) THEN
                om2 = om(j)
                ind2 = ind(j)
                xf(1) = xw(j1+i)
                yf(1) = yw(j1+i)
                ii = is-2
                DO ih = 2,np2
                    im = j1+MOD(ii+ih,j2)+1
                    xf(ih) = xw(im)
                    yf(ih) = yw(im)
                END DO
            END IF
!
!get rid of contour j
!
        CALL purge(j)
!
!reassign temporary arrays:
!
!--  1st contour:
        IF (keep1) THEN
          n = n+1
          np(n) = np1
          om(n) = om1
          ind(n) = ind1
          mb(n) = npt
          xw(npt+1:npt+np1) = xi(:np1)
          yw(npt+1:npt+np1) = yi(:np1)
          npt = npt+np1
        END IF
!--  2nd contour:
        IF (keep2) THEN
          n = n+1
          np(n) = np2
          om(n) = om2
          ind(n) = ind2
          mb(n) = npt
          xw(npt+1:npt+np2) = xf(:np2)
          yw(npt+1:npt+np2) = yf(:np2)
          npt = npt+np2
        END IF
        j = j-1
        CYCLE s1
      END DO
    END DO
!
!      next consider the merger of two distinct vortices (or the
!breakup of an annular vortex).   this process replaces two
!contours with one.
!
    js = j
    DO
      js = js+1
      IF (js > n) CYCLE s1
!
!this is the mergeability condition.
!
      IF (ind(j) /=  ind(js)) CYCLE
      CALL prop(js,js1,js2)
      DO i = 1,j2
        DO is = 1,js2

          ! See whether the node i gets within a distance dm of the
          ! line segment between is and is+1.

          ff = xw(j1+i)-xw(js1+is)
          gg = yw(j1+i)-yw(js1+is)
          p = (ff*dx(is)+gg*dy(is))/b(is)
          IF ( p*(one-p) <small ) CYCLE
!
!passing this condition, a perpendicular can be dropped onto the
!line segment from is to is+1.  check distance to line segment:
!
          dis = ABS((ff*dy(is)-gg*dx(is)))/e(is)
          IF (dis >= dm) CYCLE
!
!now we paste two contours together at node i:
!
          xw(j1+i) = xw(j1+i)+half*(-ff+p*dx(is))
          yw(j1+i) = yw(j1+i)+half*(-gg+p*dy(is))

         ! Use the xi,yi arrays as temporary storage of the
         ! combined contour, get rid of contours j and js, and finally
         ! introduce a new contour, the combined contour

          ii = i-2
          DO ih = 1,j2+1
            im = j1+MOD(ii+ih,j2)+1
            xi(ih) = xw(im)
            yi(ih) = yw(im)
          END DO
          ii = is-1
          DO ih = 1,js2
            im = js1+MOD(ii+ih,js2)+1
            xf(ih) = xw(im)
            yf(ih) = yw(im)
          END DO
          np1 = js2+j2+1
          om1 = om(j)
          ind1 = ind(j)
          CALL purge(js)
          CALL purge(j)
          n = n+1
          om(n) = om1
          ind(n) = ind1
          np(n) = np1
          mb(n) = npt
          xw(npt+1:npt+j2+1) = xi(:j2+1)
          yw(npt+1:npt+j2+1) = yi(:j2+1)
          xw(npt+j2+2:npt+np1) = xf(:js2)
          yw(npt+j2+2:npt+np1) = yf(:js2)
          npt = npt+np1
          j = j-1
                    CYCLE s1
                END DO
            END DO
        END DO
    END DO s1

    DEALLOCATE (b,e)

    CALL dimarrays  

CONTAINS

SUBROUTINE prop(ij,i1,i2)

! Get a single contour for surgery

    INTEGER, intent(in)  :: ij
    INTEGER, intent(out) :: i1,i2

    i1 = mb(ij)
    i2 = np(ij)

    DEALLOCATE (dx,dy,b,e)
    ALLOCATE   (dx(i2),dy(i2),b(i2),e(i2))

    dx(:i2-1) = xw(i1+2:i1+i2)-xw(i1+1:i1+i2-1)
    dy(:i2-1) = yw(i1+2:i1+i2)-yw(i1+1:i1+i2-1)
    dx(i2) = xw(i1+1)-xw(i1+i2)
    dy(i2) = yw(i1+1)-yw(i1+i2)
    b = dx**2+dy**2
    e = SQRT(b)

END SUBROUTINE prop

END SUBROUTINE surgery

SUBROUTINE purge(ij)

! Purge contour j from the set of contours

    INTEGER, INTENT(IN) :: ij
    INTEGER             :: npo,npto

    n = n-1
    npto = npt
    npo = np(ij)
    npt = npt-npo
    IF (ij > n) RETURN
    x(:npt-mb(ij)) = xw(mb(ij+1)+1:npto)
    y(:npt-mb(ij)) = yw(mb(ij+1)+1:npto)
    xw(mb(ij)+1:npt) = x(:npt-mb(ij))
    yw(mb(ij)+1:npt) = y(:npt-mb(ij))
    np(ij:n) = np(ij+1:n+1)
    mb(ij:n) = mb(ij+1:n+1)-npo
    ind(ij:n) = ind(ij+1:n+1)
    om(ij:n) = om(ij+1:n+1)

END SUBROUTINE purge

END MODULE surgery_
