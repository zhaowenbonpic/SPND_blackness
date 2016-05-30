SUBROUTINE GAUSCP(MMAX,SEC,ST1,VSI,P,Y,X)
! st1 inverse sig total
! vsi volume*sigma
DIMENSION SEC(mmax),ST1(mmax),VSI(mmax),P(mmax,mmax),Y(mmax),X(mmax,mmax)
!
! RHS AND MATRIX, USING THE ORIGINAL SYMMETRY
!
DO I = 1,MMAX
  DO J = 1,MMAX
    X(I,J)=ST1(J)*P(J,I)
    P(J,I)=-SEC(I)*P(J,I)
  ENDDO
  P(I,I) = P(I,I)+VSI(I)
ENDDO
!
! SWEEPING LOWER DIAGONAL MATRIX CLEAN
!
DO I = 1,MMAX
  F= 1.0/P(I,I)
  Y(I)=F*Y(I)
  DO J=1,MMAX
    X(I,J)=F*X(I,J)
  ENDDO
  IF(I.EQ.MMAX) GOTO 70
  I1 = I+1
  DO J = I1,MMAX
    P(I,J)=F*P(I,J)
  ENDDO
  !
  ! THE ACTUAL SWEEPING
  DO J = I1,MMAX
    Y(J)=Y(J)-Y(I)*P(J,I)
    DO K=1,MMAX
      X(J,K) = X(J,K)-X(I,K)*P(J,I)
    ENDDO
    DO K = I1,MMAX
      P(J,K)=P(J,K)-p(I,K)*P(J,I)
    ENDDO
  ENDDO
ENDDO
!
! SOLVE
70 IMAX = MMAX-1
IF(IMAX.EQ.0) GOTO 90
DO I=1,IMAX
  J=MMAX-i
  J1=J+1
  DO K = J1,MMAX
    Y(J)=Y(J)-Y(K)*P(J,K)
    DO KK = 1,MMAX
      X(J,KK)=X(J,KK)-X(K,KK)*P(J,K)
    ENDDO
  ENDDO
ENDDO
90 RETURN
END



