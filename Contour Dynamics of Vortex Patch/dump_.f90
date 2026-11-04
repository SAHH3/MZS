MODULE dump_

    USE dyn_com
    IMPLICIT NONE
    PRIVATE
    PUBLIC dump

CONTAINS

! Dump info to tape11

SUBROUTINE dump
    call findcentroid

    OPEN (11,FILE='fort.11',STATUS='OLD',POSITION='APPEND')

    WRITE (11,'(2I8,F20.10)') n,npt,t
    PRINT *, n,npt,t
    WRITE (11,'(3I8,F20.10)') (np(i),mb(i),ind(i),om(i), i=1,n)
    WRITE (11,'(2F20.10,L8)') xloc,yloc	
    WRITE (11,'(2F20.10,L8)')    (x(i),y(i),corn(i), i=1,npt)

    CLOSE (11)

END SUBROUTINE dump

 SUBROUTINE findcentroid

! Calculates the centroid of the largest remaining patch.

    IMPLICIT none
    REAL(kr) :: maxarea,thisarea
    INTEGER  :: maxnum
    
    real,dimension(npt) :: xave,yave
    
    dx(:npt-1) = x(2:)-x(:npt-1)
    dy(:npt-1) = y(2:)-y(:npt-1)
    xave(:npt-1)=(x(2:)+x(:npt-1))/two
    yave(:npt-1)=(y(2:)+y(:npt-1))/two
    
    DO j=1,n
    
        dx(mb(j)+np(j)) = x(mb(j)+1)-x(mb(j)+np(j))
        dy(mb(j)+np(j)) = y(mb(j)+1)-y(mb(j)+np(j))
	xave(mb(j)+np(j)) = (x(mb(j)+1)+x(mb(j)+np(j)))/two
        yave(mb(j)+np(j)) = (y(mb(j)+1)+y(mb(j)+np(j)))/two 
	 

    END DO    
        
    maxarea=zero
    maxnum=0

    DO i=1,n 
       !The contour must have the same index as original patch
       IF (ind(i).ne.1) CYCLE 
       
   
       
       !calculate the area by Green's theorem.
       thisarea=zero
       DO j=mb(i)+1,mb(i)+np(i) 
         thisarea=thisarea+half*(xave(j)*dy(j) - yave(j)*dx(j))
       ENDDO
       
       thisarea=abs(thisarea)
       
       IF (thisarea.gt.maxarea) THEN 
          maxarea=thisarea
          maxnum=i
       ENDIF
    ENDDO

    !print*,"Area =",maxarea
    !print*,"on contour ",maxnum

    !find the centroid co-ordinates
    xloc=zero
    yloc=zero
    DO j=mb(maxnum)+1,mb(maxnum)+np(maxnum)
       xloc=xloc+half*half*(xave(j)*xave(j)*dy(j))-half*xave(j)*yave(j)*dx(j)
       yloc=yloc-half*half*(yave(j)*yave(j)*dx(j))+half*xave(j)*yave(j)*dy(j)
    ENDDO
    xloc=xloc/maxarea
    yloc=yloc/maxarea
    
END SUBROUTINE findcentroid   

END MODULE dump_
