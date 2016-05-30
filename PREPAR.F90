subroutine prepar(rad,vol,sig,sigs,sigr,sig1,sec,sigv)
! sigv = total sig times vol
! sig1 = the inverse of total sig
! SEC  = number of secondaries
dimension rad(mmax),vol(mmax),sig(mmax),sigs(mmax),sigr(mmax),sig1(mmax),sec(mmax),sigv(mmax)
common /CPR001/ nin,NOT,nx1(2)
common /CPR002/ mmax,ig,nx203,VTOT,SB,NX206
logical swerr,swhite
DATA PI / 3.1415927/
data small /1e-8/

1010 format(6e12.4)
2010 format('number of regions',i4,/,'number of gauss points', i4,/)
2020 format ('   REGION   RADIUS      VOLUME      SIGTOT      SIGSCAT    SIGREMO   SECONDARIES',//)
2030 FORMAT(I5,1X,1P6E12.4)
2040 FORMAT(' TOTAL VOLUME =',1PE11.4,'OUTER SURFACE = ',E11.4)

READ(NiN,*) (RAD(M),M=1,MMAX)
READ(NiN,*) (SIG(M),M=1,MMAX)
READ(NiN,*) (SIGS(M),M=1,MMAX)

SWERR = .FALSE.
VA = 0.0
DO M = 1,MMAX
  R2 = RAD(M)*RAD(M)
  VOL(M) = PI*(R2-VA)
  VA = R2
  IF(VOL(M).LE.0.0) SWERR = .TRUE.
  SIGR(M) = SIG(M)-SIGS(M)
  SIG(M) = SIG(M) + SMALL
  SIG1(M) = 1.0/SIG(M)
  SEC(M) = SIGS(M)*SIG1(M)
  SIGV(M) = SIG(M)*VOL(M)
ENDDO

SB = 2.0*PI*RAD(MMAX)
VTOT = 0.5*SB*RAD(MMAX)

WRITE(NOT,2010) MMAX,IG
WRITE(NOT,2020)
DO M = 1,MMAX
  WRITE(NOT,2030) M,RAD(M),VOL(M),SIG(M),SIGS(M),SIGR(M),SEC(M)
ENDDO
WRITE(NOT,2040) VTOT,SB
IF(SWERR) WRITE(NOT,"('WRONG INPUT')")
RETURN
END




