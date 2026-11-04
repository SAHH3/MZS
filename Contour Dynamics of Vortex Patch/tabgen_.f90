MODULE tabgen_

! Generate the table coefficients for the derivative of the g function.

    USE dyn_com
    IMPLICIT NONE
    PRIVATE
    PUBLIC tabgen

CONTAINS

SUBROUTINE tabgen

! Since bessk1 is declared below in this module 
! the following statement would be erroneous:
!   REAL(kr) :: bessk1

    REAL(kr)           :: r

    drdi = rdi/dr
    rdsq = one/rdi**2
    bdsq = (barg/rdi)**2

    agf(0) = zero
    DO j = 1,ntab
        r = j*dr
        agf(j) = r*bessk1(r) - one
        ! bessk1 is the Bessel function K_1 from Numerical Recipes.
    END DO

    DO j = 0,ntab-1
        bgf(j) = half*(agf(j+1)-agf(j))
        agf(j) = half*agf(j)-j*bgf(j)
    END DO

    bgf(ntab) = zero
    agf(ntab) = -half
      
END SUBROUTINE tabgen

REAL(kr) FUNCTION bessk1(x)

! Modified Bessel function K_1(x) for real x>0, using bessi1.f.


! Since bessi1 is declared below in this module 
! the following statement would be erroneous:
!   REAL(kr) :: bessi1

    REAL(kr)      ::  x,y,p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
    SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7
    DATA p1,p2,p3,p4,p5,p6,p7/1.0D0,0.15443144D0,-0.67278579D0,         &
           -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
    DATA q1,q2,q3,q4,q5,q6,q7/1.25331414D0,0.23498619D0,-0.3655620D-1,  &
            0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/

    IF (x.le.two) THEN
        y = (x/two)**2
        bessk1 = (LOG(x/two)*bessi1(x))+(one/x)*(p1+y*(p2+              &
                                    y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
    ELSE
        y = two/x
        bessk1 = (EXP(-x)/SQRT(x))*(q1+y*(q2+y*(q3+                     &
                                          y*(q4+y*(q5+y*(q6+y*q7))))))
    END IF

END FUNCTION bessk1


REAL(kr) FUNCTION bessi1(x)

! Modified Bessel function I_1(x) for any real x

    REAL(kr) :: x,ax,y,p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
    SAVE        p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
    REAL(kr), PARAMETER :: three75 = 3.75_kr
    DATA p1,p2,p3,p4,p5,p6,p7/0.5D0,0.87890594D0,0.51498869D0,    &
          0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
    DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228D0,-0.3988024D-1,   &
          -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,    &
          -0.2895312D-1,0.1787654D-1,-0.420059D-2/

    IF ( ABS(x) < three75 ) THEN
        y = (x/three75)**2
        bessi1 = x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))

    ELSE
        ax = ABS(x)
        y = three75/ax
        bessi1 = (EXP(ax)/SQRT(ax))*(q1+y*(q2+y*(q3+y*(q4+        &
                           y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
    END IF

END FUNCTION bessi1

END MODULE tabgen_




