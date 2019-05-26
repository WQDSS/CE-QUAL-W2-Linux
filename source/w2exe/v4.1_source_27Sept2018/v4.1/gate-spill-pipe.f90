
!***********************************************************************************************************************************
!**                                          S U B R O U T I N E   G A T E  F L O W                                               **
!***********************************************************************************************************************************

SUBROUTINE GATE_FLOW
  USE STRUCTURES; USE GLOBAL; USE GEOMC
  IMPLICIT NONE
    INTEGER :: JG,ISUB,IGT
    REAL(R8)    :: ELIU,ELID,HTAIL,HENERGY,DLEL


  DO JG=1,NGT
   IF(DYNGTC(JG) == '    FLOW')THEN    
    QGT(JG) = BGT(JG)
   ELSE

 !   ELIU =  ELWS(IUGT(JG))       !EL(KTWB(JWUGT(JG)),IUGT(JG))-Z(IUGT(JG))*COSA(BS(JWUGT(JG)))
    IF (LATERAL_GATE(JG)) THEN
      ELIU =  ELWS(IUGT(JG))      !EL(KTWB(JWUGT(JG)),IUGT(JG))-Z(IUGT(JG))*COSA(BS(JWUGT(JG)))
    ELSE
      !ELIU=ELWS(IUGT(JG))-SINA(JBUGT(JG))*DLX(IUGT(JG))*0.5     !EL(KTWB(JWUGT(JG)),IUGT(JG))-Z(IUGT(JG))*COSA(BS(JWUGT(JG)))-SINA(JBUGT(JG))*DLX(IUGT(JG))*0.5
      ELIU= ELWS(IUGT(JG)) + (ELWS(IUGT(JG))-ELWS(IUGT(JG)-1))/(0.5*(DLX(IUGT(JG))+DLX(IUGT(JG)-1)))*DLX(IUGT(JG))*0.5     ! LINEAR INTERPOLATION OF THE WATER LEVEL TO THE EDGE
    END IF
    IF (IDGT(JG) /= 0)THEN
      IF (US(JBDGT(JG)) /= IDGT(JG)) THEN
        ELID = ELWS(IDGT(JG))     !EL(KTWB(JWDGT(JG)),IDGT(JG))-Z(IDGT(JG))*COSA(BS(JWDGT(JG)))
      ELSE
        !ELID = ELWS(IDGT(JG))+SINA(JBDGT(JG))*DLX(IDGT(JG))*0.5          !EL(KTWB(JWDGT(JG)),IDGT(JG))-Z(IDGT(JG))*COSA(BS(JWDGT(JG)))+SINA(JBDGT(JG))*DLX(IDGT(JG))*0.5
        ELID = ELWS(IDGT(JG)) - (ELWS(IDGT(JG)+1) - ELWS(IDGT(JG)))/(0.5*(DLX(IDGT(JG))+DLX(IDGT(JG)+1)))*DLX(IDGT(JG))*0.5  
      END IF
    ELSE
      ELID = -100.0
    END IF
    IF (BGT(JG) /= 0.0) THEN
      IF (ELID > EGT(JG) .OR. ELIU > EGT(JG)) THEN
        ISUB = 0
        IF (A2GT(JG) /= 0.0 .AND. IDGT(JG) /= 0) THEN                                    ! SW 8/21/2013
          HTAIL =  ELID-EGT(JG)                                                          ! SW 5/10/05
          IF (HTAIL > 0) THEN
            HENERGY = (U(KTWB(JWUGT(JG)),IUGT(JG))**2)/(2.0*G)+ELIU-EGT(JG)              ! SW 5/10/05
            IF (HTAIL/HENERGY > 0.67) ISUB = 1
          END IF
        END IF
        IGT = 0
        IF (BGT(JG) >= 0.8*(ELIU-EGT(JG)) .AND. GTA1(JG) /= 0.0) IGT = 1
        IF (IGT == 0) THEN
          IF (ISUB == 0) THEN
            DLEL = ELIU-EGT(JG)
            IF (A2GT(JG) == 0.0 .AND. G2GT(JG) /= 0.0) DLEL = ELIU-G2GT(JG)
            IF (DLEL < 0.0) THEN
              DLEL    = -DLEL
              QGT(JG) = -A1GT(JG)*(DLEL**B1GT(JG))*BGT(JG)**G1GT(JG)
            ELSE
              QGT(JG) =  A1GT(JG)*(DLEL**B1GT(JG))*BGT(JG)**G1GT(JG)
            END IF
          ELSE IF (ELID > ELIU) THEN
            DLEL    =  ELID-ELIU
            QGT(JG) = -A2GT(JG)*DLEL**B2GT(JG)*BGT(JG)**G2GT(JG)
          ELSE
            DLEL    =  ELIU-ELID
            QGT(JG) =  A2GT(JG)*DLEL**B2GT(JG)*BGT(JG)**G2GT(JG)
          END IF
        ELSE IF (ISUB == 0) THEN
          DLEL = ELIU-EGT(JG)
          IF (ELID > EGT(JG)) DLEL = ELIU-ELID
          IF (DLEL < 0.0) THEN
            DLEL    = -DLEL
            QGT(JG) = -GTA1(JG)*DLEL**GTB1(JG)
          ELSE
            QGT(JG) =  GTA1(JG)*DLEL**GTB1(JG)
          END IF
        ELSE IF (ELID > ELIU) THEN
          DLEL    =  ELID-ELIU
          QGT(JG) = -GTA2(JG)*DLEL**GTB2(JG)
        ELSE
          DLEL    =  ELIU-ELID
          QGT(JG) =  GTA2(JG)*DLEL**GTB2(JG)
        END IF
      ELSE
        QGT(JG) = 0.0
      END IF
    ELSE
      QGT(JG) = 0.0
    END IF
   endif
  END DO
END SUBROUTINE GATE_FLOW

!***********************************************************************************************************************************
!**                                         S U B R O U T I N E   S P I L L W A Y  F L O W                                        **
!***********************************************************************************************************************************

SUBROUTINE SPILLWAY_FLOW
  USE STRUCTURES; USE GLOBAL; USE GEOMC
    INTEGER :: JS,ISUB
    REAL(R8):: ELIU,ELID,HTAIL,HENERGY,DLEL

  DO JS=1,NSP
    IF (LATERAL_SPILLWAY(JS)) THEN
       ELIU =  ELWS(IUSP(JS))                                        !EL(KTWB(JWUSP(JS)),IUSP(JS))-Z(IUSP(JS))*COSA(BS(JWUSP(JS)))
    ELSE
      ! ELIU =  ELWS(IUSP(JS))-SINA(JBUSP(JS))*DLX(IUSP(JS))*0.5      !EL(KTWB(JWUSP(JS)),IUSP(JS))-Z(IUSP(JS))*COSA(BS(JWUSP(JS)))-SINA(JBUSP(JS))*DLX(IUSP(JS))*0.5
        ELIU= ELWS(IUSP(JS)) + (ELWS(IUSP(JS))-ELWS(IUSP(JS)-1))/(0.5*(DLX(IUSP(JS))+DLX(IUSP(JS)-1)))*DLX(IUSP(JS))*0.5     ! LINEAR INTERPOLATION OF THE WATER LEVEL TO THE EDGE
    END IF
    IF (IDSP(JS) /= 0) THEN
      IF (US(JBDSP(JS)) /= IDSP(JS)) THEN
         ELID = ELWS(IDSP(JS))                                       !EL(KTWB(JWDSP(JS)),IDSP(JS))-Z(IDSP(JS))*COSA(BS(JWDSP(JS)))
      ELSE
      !   ELID = ELWS(IDSP(JS))+SINA(JBDSP(JS))*DLX(IDSP(JS))*0.5     !EL(KTWB(JWDSP(JS)),IDSP(JS))-Z(IDSP(JS))*COSA(BS(JWDSP(JS)))+SINA(JBDSP(JS))*DLX(IDSP(JS))*0.5
          ELID = ELWS(IDSP(JS)) - (ELWS(IDSP(JS)+1) - ELWS(IDSP(JS)))/(0.5*(DLX(IDSP(JS))+DLX(IDSP(JS)+1)))*DLX(IDSP(JS))*0.5  
      END IF
    ELSE
      ELID = -1.0
    END IF
    IF (ELID >= ESP(JS) .OR. ELIU >= ESP(JS)) THEN
      ISUB = 0
      IF (A2SP(JS) /= 0.0 .AND. IDSP(JS) /= 0) THEN
        HTAIL   =  ELID-ESP(JS)                                                 ! SW 5/10/05
        IF (HTAIL > 0) THEN
          HENERGY = (U(KTWB(JWUSP(JS)),IUSP(JS))**2)/(2.0*G)+ELIU-ESP(JS)       ! SW 5/10/05
          IF (HTAIL/HENERGY > 0.67) ISUB = 1
        END IF
      END IF
      IF (ISUB == 0) THEN
        DLEL = ELIU-ESP(JS)
        IF (DLEL < 0.0) THEN
          DLEL    = -DLEL
          QSP(JS) = -A1SP(JS)*DLEL**B1SP(JS)
        ELSE
          QSP(JS) =  A1SP(JS)*DLEL**B1SP(JS)
        END IF
      ELSE IF (ELID > ELIU) THEN
        DLEL    =  ELID-ELIU
        QSP(JS) = -A2SP(JS)*DLEL**B2SP(JS)
      ELSE
        DLEL    =  ELIU-ELID
        QSP(JS) =  A2SP(JS)*DLEL**B2SP(JS)
      END IF
    ELSE
      QSP(JS) = 0.0
    END IF
  END DO
END SUBROUTINE SPILLWAY_FLOW

!***********************************************************************************************************************************
!**                                            S U B R O U T I N E   P I P E  F L O W                                             **
!***********************************************************************************************************************************

SUBROUTINE PIPE_FLOW_INITIALIZE
  USE GLOBAL; USE GEOMC; USE STRUCTURES; USE SCREENC, ONLY: NIT
  REAL(R8) :: DTQ,DLTX,EL1,EL2,HIE,EPS,DCHECK,D1,D2,DTEST,DCRIT,DEPTHCRIT,VTOT,TOTT
  INTEGER  :: JP     !,NIT

  SAVE

  ALLOCATE (BEGIN(NPI), WLFLAG(NPI), VMAX(NPI))
  QOLD   =  0.01;  VMAX   =  0.01
  BEGIN  = .TRUE.; WLFLAG = .TRUE.
RETURN

ENTRY PIPE_FLOW      !(NIT)
  DTQ = DLT/10.0
  DO JP=1,NPI
    DIA   = WPI(JP)
    CLEN  = DLXPI(JP)
    FMAN  = FPI(JP)
    CLOSS = FMINPI(JP)
    UPIE  = EUPI(JP)
    DNIE  = EDPI(JP)
    DLTX  = CLEN/(REAL(NC-1)*0.5)
    IF (LATERAL_PIPE(JP)) THEN
      EL1   = ELWS(IUPI(JP))                                                                 !EL(KTWB(JWUPI(JP)),IUPI(JP))-Z(IUPI(JP))*COSA(JBUPI(JP))
    ELSE
     ! EL1   = ELWS(IUPI(JP))-SINA(JBDPI(JP))*DLX(IUPI(JP))*0.5                               !EL(KTWB(JWUPI(JP)),IUPI(JP))-Z(IUPI(JP))*COSA(JBUPI(JP))-SINA(JBDPI(JP))*DLX(IUPI(JP))*0.5
       EL1 = ELWS(IUPI(JP)) + (ELWS(IUPI(JP))-ELWS(IUPI(JP)-1))/(0.5*(DLX(IUPI(JP))+DLX(IUPI(JP)-1)))*DLX(IUPI(JP))*0.5     ! LINEAR INTERPOLATION OF THE WATER LEVEL TO THE EDGE
    END IF
    IF (IDPI(JP) /= 0) THEN
      IF (US(JBDPI(JP)) /= IDPI(JP)) THEN
        EL2   = ELWS(IDPI(JP))                                                               !EL(KTWB(JWDPI(JP)),IDPI(JP))-Z(IDPI(JP))*COSA(JBDPI(JP))
      ELSE
       ! EL2   = ELWS(IDPI(JP))+SINA(JBDPI(JP))*DLX(IDPI(JP))*0.5                             !EL(KTWB(JWDPI(JP)),IDPI(JP))-Z(IDPI(JP))*COSA(JBDPI(JP))+SINA(JBDPI(JP))*DLX(IDPI(JP))*0.5
         EL2 = ELWS(IDPI(JP)) - (ELWS(IDPI(JP)+1) - ELWS(IDPI(JP)))/(0.5*(DLX(IDPI(JP))+DLX(IDPI(JP)+1)))*DLX(IDPI(JP))*0.5  
      END IF
    ELSE
      EL2 = -1.0
    END IF
    HIE = MAX(UPIE,DNIE)
    IF (DIA == 0.0) THEN
      QPI(JP)    =  0.0
      WLFLAG(JP) = .TRUE.
      GO TO 140
    END IF
    EPS = 0.001
    IF ((HIE+EPS) >= EL1 .AND. (HIE+EPS) >= EL2) THEN
      QPI(JP)    =  0.0
      WLFLAG(JP) = .TRUE.
      GO TO 140
    END IF
    IF (EL1 > EL2) THEN
      DCHECK = EL1-UPIE
    ELSE
      DCHECK = EL2-DNIE
    END IF
    IF (DCHECK < 0.02) THEN
      QPI(JP)    =  0.0
      WLFLAG(JP) = .TRUE.
      GO TO 140
    END IF
    IF (ABS(QOLD(JP)) < 0.001) QOLD(JP) = 0.001
    IF (EL1 >= (UPIE+DIA) .AND. EL2 >= (DNIE+DIA)) THEN
      D1 = EL1
      D2 = EL2
      GO TO 120
    END IF
    IF (EL1 > EL2) THEN
      DTEST = EL2-DNIE
    ELSE
      DTEST = EL1-UPIE
    END IF
    DCRIT = DEPTHCRIT(ABS(QOLD(JP)))
    IF (DTEST <= DCRIT) THEN
      IF (EL1 <= EL2) THEN
        D1 = UPIE+DCRIT
        D2 = EL2
      ELSE
        D1 = EL1
        D2 = DNIE+DCRIT
      END IF
      VTOT = 0.0
      TOTT = 0.0
110   CONTINUE
      IF (NIT /= 0) THEN
        DTQ = OMEGA*DLTX/VMAX(JP)
        IF (DTQ > (DLT-TOTT)) THEN
          DTQ = DLT-TOTT
        ELSE IF ((2.0*DTQ) > (DLT-TOTT)) THEN
          DTQ = (DLT-TOTT)*0.5
        END IF
      END IF
      CALL OPEN_CHANNEL (D1,D2,QPI(JP),JP,DTQ)
      DCRIT = DEPTHCRIT(ABS(QPI(JP)))
      IF (EL1 <= EL2) THEN
        D1 = UPIE+DCRIT
      ELSE
        D2 = DNIE+DCRIT
      END IF
      VTOT = VTOT+DTQ*QPI(JP)
      TOTT = DTQ+TOTT
      IF (TOTT < (DLT-EPS2)) GO TO 110
      QPI(JP) = VTOT/DLT
      GO TO 140
    END IF
    D1 = EL1
    D2 = EL2
120 CONTINUE
    TOTT = 0.0
    VTOT = 0.0
130 CONTINUE
    IF (NIT /= 0) THEN
      DTQ = OMEGA*DLTX/VMAX(JP)
      IF (DTQ > (DLT-TOTT)) THEN
        DTQ = DLT-TOTT
      ELSE IF ((2.0*DTQ) > (DLT-TOTT)) THEN
        DTQ = (DLT-TOTT)*0.5
      END IF
    END IF
    CALL OPEN_CHANNEL (D1,D2,QPI(JP),JP,DTQ)
    VTOT = VTOT+DTQ*QPI(JP)
    TOTT = DTQ+TOTT
    IF (TOTT < (DLT-EPS2)) GO TO 130
    QPI(JP) = VTOT/DLT
140 CONTINUE
    QOLD(JP) = QPI(JP)
    IF (QPI(JP) == 0.0) WLFLAG(JP) = .TRUE.
  END DO
RETURN
ENTRY DEALLOCATE_PIPE_FLOW
  DEALLOCATE (BEGIN, WLFLAG, VMAX)
RETURN
END SUBROUTINE PIPE_FLOW_INITIALIZE

!***********************************************************************************************************************************
!**                                        S U B R O U T I N E   O P E N  C H A N N E L                                           **
!***********************************************************************************************************************************

SUBROUTINE OPEN_CHANNEL_INITIALIZE
  USE GLOBAL; USE STRUCTURES
  REAL(R8), PARAMETER :: THETA=0.55
  REAL(R8)    :: DLTX,PHI,VTOT,SLOPE,DIST,BEPR1,BEPR2,BC1,EL1,BC2,EL2,WLSLOPE,DLTX2,BAR1,BAREA,RAD1,RAD2
  REAL(R8)    :: TWIDTH,VAVG,QSUM,QAVG,QOUT,WETPER,BAR2


! Type declarations

  INTEGER                              :: J,IC,N,NP,NQCNT
  REAL(R8)                                 :: DT,D
  REAL(R8),    ALLOCATABLE, DIMENSION(:)   :: Y,   B,   V,   CAREA, TOPW,  BELEV, Q, VOLD, YOLD     ! CB 10/4/07
  REAL(R8),    ALLOCATABLE, DIMENSION(:)   :: YT,  VT, VPR, YPR, TAREA, TOPWT, RT
  REAL(R8),    ALLOCATABLE, DIMENSION(:,:) :: DAA, AL
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: INDX
  LOGICAL                              :: SMOOTH_WATER_LEVELS      !, OPENWRN
  SAVE

! Allocation declarations

  ALLOCATE (Y(NN),    V(NN),     CAREA(NN),  TOPW(NN),   BELEV(NN),  Q(NN),     VOLD(NN), YOLD(NN), B(NN))         ! CB 10/4/07
  ALLOCATE (YT(NN),   VT(NN),    VPR(NN),    YPR(NN),    TAREA(NN),  TOPWT(NN), RT(NN),   INDX(NN))
  ALLOCATE (AL(NN,2), DAA(NN,NN))
RETURN

ENTRY OPEN_CHANNEL (EL1,EL2,QOUT,IC,DT)

! Variable initializtion

  B     = 0.0; Y     = 0.0; V = 0.0; VT = 0.0; YT = 0.0; RT = 0.0; DAA = 0.0; YPR = 0.0; VPR = 0.0; TOPW = 0.0; TOPWT = 0.0
  CAREA = 0.0; TAREA = 0.0
  BELEV(1)  = UPIE
  BELEV(NC) = DNIE
  PHI       = ASIN((UPIE-DNIE)/CLEN)
  DLTX      = CLEN/(REAL(NC-1)*0.5)
  DO J=2,NC-1
    DLTX2    =  DLTX*0.5
    SLOPE    = (UPIE-DNIE)/CLEN
    DIST     = (REAL(J-1)*DLTX2)
    BELEV(J) =  UPIE-SLOPE*DIST
  END DO
  BEPR1 =  UPIE+SLOPE*DLTX2
  BEPR2 =  DNIE-SLOPE*DLTX2
  BC1   = (EL1-BEPR1)*COS(PHI)
  IF (BC1 <= 0.0) BC1 = EL1-UPIE
  BC2 = (EL2-BEPR2)*COS(PHI)
  IF (BC2 <= 0.0) BC2 = EL2-DNIE
  IF (.NOT. BEGIN(IC)) THEN
    IF (WLFLAG(IC)) THEN
      DO J=2,NC-1,2
        WLSLOPE = ((BC1-BC2)/(CLEN+DLTX))*DCOS(PHI)
        DIST    = (REAL(J-1)*0.5*DLTX)+DLTX2
        Y(J)    =  BC1-WLSLOPE*DIST
        YT(J)   =  Y(J)
        DTP(IC) =  DT
      END DO
    ELSE
      DO I=2,NC-1,2
        Y(I)  = YS(I,IC)
        YT(I) = YST(I,IC)
      END DO
    END IF
  END IF
  DO I=1,NC,2
    V(I)  = VS(I,IC)
    VT(I) = VST(I,IC)
  END DO
  IF (BEGIN(IC)) THEN
    BEGIN(IC) = .FALSE.
    DO J=2,NC-1,2
      WLSLOPE = ((BC1-BC2)/(CLEN+DLTX))*DCOS(PHI)
      DIST    = (REAL(J-1)*0.5*DLTX)+DLTX2
      Y(J)    =  BC1-WLSLOPE*DIST
      YT(J)   =  Y(J)
      DTP(IC) =  DT
    END DO
    DO J=1,NC,2
      V(J)  = 0.0
      VT(J) = V(J)
    END DO
!    OPENWRN = .TRUE.
  END IF
  SMOOTH_WATER_LEVELS = .FALSE.
  DO N=1,NC,2
    IF (N == NC) THEN
      BAR1 = BAREA(BC2,DIA)
      RAD1 = BAR1/WETPER(BC2,DIA)
    ELSE
      BAR1 = BAREA(Y(N+1),DIA)
      RAD1 = BAR1/WETPER(Y(N+1),DIA)
    END IF
    IF (N == 1) THEN
      BAR2 = BAREA(BC1,DIA)
      RAD2 = BAR2/WETPER(BC1,DIA)
    ELSE
      BAR2 = BAREA(Y(N-1),DIA)
      RAD2 = BAR2/WETPER(Y(N-1),DIA)
    END IF
    RT(N) = (RAD1+RAD2)*0.5
  END DO
  DO N=2,NC-1,2
    TAREA(N) = BAREA(Y(N),DIA)
    TOPWT(N) = TWIDTH(Y(N),DIA)
    CAREA(N) = BAREA(Y(N),DIA)
  END DO

! Projected water levels and velocities

  DO J=1,NC,2
    VPR(J) = V(J)+DT*(V(J)-VT(J))/DTP(IC)
  END DO
  DO J=2,NC-1,2
    YPR(J) = Y(J)+DT*(Y(J)-YT(J))/DTP(IC)
  END DO

! Matrix setup

  VTOT = 0.0
  DO J=1,NC,2
    VTOT = VTOT+V(J)
  END DO
  VAVG = VTOT/(REAL(NC-1)*0.5)

! Continuity

  DO N=2,NC-1,2
    VPR(N) = (VPR(N-1)+VPR(N+1))*0.5D0
    V(N)   = (V(N-1)+V(N+1))*0.5D0
    IF (N /= 2) DAA(N,N-2) = -THETA*(DT/DLTX)*(VPR(N)*0.5)
    DAA(N,N-1) = -THETA*(DT/DLTX)*(TAREA(N)/TOPWT(N))
    DAA(N,N)   =  1.0D0
    DAA(N,N+1) =  THETA*(DT/DLTX)*(TAREA(N)/TOPWT(N))
    IF (N /= NC-1) DAA(N,N+2) = THETA*(DT/DLTX)*(VPR(N)*0.5D0)
    IF (N == 2) THEN
      B(N) = Y(N)-(1.0D0-THETA)*(DT/DLTX)*(TAREA(N)/TOPWT(N))*(V(N+1)-V(N-1))-(1.0D0-THETA)*(DT/DLTX)*(V(N)*0.5D0)*(Y(N+2)-BC1)          &
             +THETA*(DT/DLTX)*(VPR(N)*0.5D0)*BC1
    ELSE IF (N == NC-1) THEN
      B(N) = Y(N)-(1.0D0-THETA)*(DT/DLTX)*(TAREA(N)/TOPWT(N))*(V(N+1)-V(N-1))-(1.0D0-THETA)*(DT/DLTX)*(V(N)*0.5D0)*(BC2-Y(N-2))          &
             -THETA*(DT/DLTX)*(VPR(N)*0.5D0)*BC2
    ELSE
      B(N) = Y(N)-(1.0D0-THETA)*(DT/DLTX)*(TAREA(N)/TOPWT(N))*(V(N+1)-V(N-1))-(1.0D0-THETA)*(DT/DLTX)*(V(N)*0.5D0)*(Y(N+2)-Y(N-2))
    END IF
  END DO
  IF (VAVG > 0.0 .OR. (VAVG == 0.0 .AND. EL1 > EL2)) THEN

!** Momentum

    DO N=1,NC,2
      IF (N /= 1) THEN
        DAA(N,N-2) = -THETA*(DT/DLTX)*VPR(N)
        DAA(N,N-1) = -THETA*(DT/DLTX)*G*DCOS(PHI)
      END IF
      DAA(N,N) = 1.0+THETA*DT*G*(FMAN**2)*DABS(VPR(N))/(RT(N)**(4.0/3.0))+THETA*(DT/DLTX)*VPR(N)+THETA*(CLOSS*0.5D0)*(DT/CLEN)        &
                 *DABS(VPR(N))
      IF (N /= NC) DAA(N,N+1) = THETA*(DT/DLTX)*G*DCOS(PHI)
      IF (N == 1) THEN
        B(N) = V(N)-(1.0D0-THETA)*(DT/DLTX)*G*(Y(N+1)-BC1)*DCOS(PHI)-(1.0D0-THETA)*V(N)*(DT/DLTX)*V(N)-(1.0D0-THETA)*DT*G*(FMAN**2)       &
               /(RT(N)**(4.0/3.0))*V(N)*DABS(V(N))+DT*G*DSIN(PHI)-(1.0D0-THETA)*(DT/CLEN)*(CLOSS*0.5D0)*V(N)*DABS(V(N))+THETA*(DT/DLTX)   &
               *G*DCOS(PHI)*BC1
      ELSE IF (N == NC) THEN
        B(N) = V(N)-(1.0D0-THETA)*(DT/DLTX)*G*(BC2-Y(N-1))*DCOS(PHI)-(1.0D0-THETA)*V(N)*(DT/DLTX)*(V(N)-V(N-2))-(1.0D0-THETA)             &
               *DT*G*(FMAN**2)/(RT(N)**(4.0/3.0))*V(N)*DABS(V(N))+DT*G*DSIN(PHI)-(1.0D0-THETA)*(DT/CLEN)*(CLOSS*0.5D0)*V(N)*DABS(V(N))    &
               -THETA*(DT/DLTX)*G*DCOS(PHI)*BC2
      ELSE
        B(N) = V(N)-(1.0D0-THETA)*(DT/DLTX)*G*(Y(N+1)-Y(N-1))*COS(PHI)-(1.0D0-THETA)*V(N)*(DT/DLTX)*(V(N)-V(N-2))-(1.0D0-THETA)          &
               *DT*G*(FMAN**2)/(RT(N)**(4.0/3.0))*V(N)*DABS(V(N))+DT*G*DSIN(PHI)-(1.0D0-THETA)*(DT/CLEN)*(CLOSS*0.5D0)*V(N)*DABS(V(N))
      END IF
    END DO
  ELSE
    DO N=1,NC,2
      IF (N /= NC) THEN
        DAA(N,N+2) = THETA*(DT/DLTX)*VPR(N)
        DAA(N,N+1) = THETA*(DT/DLTX)*G*DCOS(PHI)
      END IF
      DAA(N,N) = 1.0+THETA*DT*G*(FMAN**2)*DABS(VPR(N))/(RT(N)**(4.0/3.0))-THETA*(DT/DLTX)*VPR(N)+THETA*(CLOSS*0.5D0)*(DT/CLEN)        &
                 *DABS(VPR(N))
      IF (N /= 1) DAA(N,N-1) = -THETA*(DT/DLTX)*G*DCOS(PHI)
      IF (N == NC) THEN
        B(N) = V(N)-(1.0D0-THETA)*(DT/DLTX)*G*(BC2-Y(N-1))*DCOS(PHI)-(1.0-THETA)*V(N)*(DT/DLTX)*(-V(N))-(1.0D0-THETA)*DT*G*(FMAN**2)    &
               /(RT(N)**(4.0/3.0))*V(N)*DABS(V(N))+DT*G*DSIN(PHI)-(1.0-THETA)*(DT/CLEN)*(CLOSS*0.5)*V(N)*DABS(V(N))-THETA*(DT/DLTX)   &
               *G*DCOS(PHI)*BC2
      ELSE IF (N == 1) THEN
        B(N) = V(N)-(1.0D0-THETA)*(DT/DLTX)*G*(Y(N+1)-BC1)*DCOS(PHI)-(1.0-THETA)*V(N)*(DT/DLTX)*(V(N+2)-V(N))-(1.0D0-THETA)             &
               *DT*G*(FMAN**2)/(RT(N)**(4.0/3.0))*V(N)*ABS(V(N))+DT*G*SIN(PHI)-(1.0-THETA)*(DT/CLEN)*(CLOSS*0.5D0)*V(N)*DABS(V(N))    &
               +THETA*(DT/DLTX)*G*DCOS(PHI)*BC1
      ELSE
        B(N) = V(N)-(1.0D0-THETA)*(DT/DLTX)*G*(Y(N+1)-Y(N-1))*DCOS(PHI)-(1.0D0-THETA)*V(N)*(DT/DLTX)*(V(N+2)-V(N))-(1.0D0-THETA)          &
               *DT*G*(FMAN**2)/(RT(N)**(4.0/3.0))*V(N)*DABS(V(N))+DT*G*DSIN(PHI)-(1.0D0-THETA)*(DT/CLEN)*(CLOSS*0.5D0)*V(N)*DABS(V(N))
      END IF
    END DO
  END IF
  NP = NN
  CALL LUDCMP (DAA,NC,NP,INDX,D)
  CALL LUBKSB (DAA,NC,NP,INDX,B)
  DO I=2,NC-1,2
    YOLD(I)   = Y(I)
    YST(I,IC) = Y(I)
  END DO
  DO I=2,NC-1,2
    Y(I) = B(I)
  END DO

! Smooth water levels

  DO I=2,NC-1,2
    IF (Y(I) <= 0.0) THEN
 !     IF (OPENWRN) THEN
 !       OPEN (391,FILE='culvert.wrn',STATUS='unknown')
 !       OPENWRN = .FALSE.
 !     END IF
      SMOOTH_WATER_LEVELS = .TRUE.
    END IF
  END DO
  IF (SMOOTH_WATER_LEVELS) THEN
    DO J=2,NC-1,2
      WLSLOPE = ((BC1-BC2)/(CLEN+DLTX))*DCOS(PHI)
      DIST    = (REAL(J-1)*0.5D0*DLTX)+DLTX2
      Y(J)    =  BC1-WLSLOPE*DIST
    END DO
  !  WRITE (391,10010) IC, JDAY
    SMOOTH_WATER_LEVELS = .FALSE.
  END IF

! Flows

  NQCNT = 0
  QSUM  = 0.0
  DO I=1,NC,2
    VOLD(I)   = V(I)
    VST(I,IC) = V(I)
    V(I)      = B(I)
    IF (I == NC) THEN
      BAR1 = BAREA(BC2,DIA)
    ELSE
      BAR1 = BAREA(Y(I+1),DIA)
    END IF
    IF (I == 1) THEN
      BAR2 = BAREA(BC1,DIA)
    ELSE
      BAR2 = BAREA(Y(I-1),DIA)
    END IF
    CAREA(I) = (BAR1+BAR2)*0.5D0
    Q(I)     =  V(I)*CAREA(I)
    NQCNT    =  NQCNT+1
    QSUM     =  QSUM+Q(I)
  END DO
  QAVG = QSUM/REAL(NQCNT)
  DO I=2,NC-1,2
    YS(I,IC) = Y(I)
  END DO
  VMAX(IC) = 0.0
  DO I=1,NC,2
    VS(I,IC) = V(I)
    VMAX(IC) = MAX(ABS(V(I)),VMAX(IC))
  END DO
  DTP(IC)    =  DT
  QOUT       =  QAVG
  QOLD(IC)   =  QOUT
  WLFLAG(IC) = .FALSE.
10010 FORMAT ('water levels for culvert ',I3,' on Julian Day ',F10.3,' are <= 0 - predictions have been smoothed')
RETURN
ENTRY DEALLOCATE_OPEN_CHANNEL
  DEALLOCATE (Y, V, CAREA, TOPW, BELEV, Q, VOLD, YOLD, B, YT, VT, VPR, YPR, TAREA, TOPWT, RT, INDX, AL, DAA)      ! CB 10/4/07
RETURN
END SUBROUTINE OPEN_CHANNEL_INITIALIZE

!***********************************************************************************************************************************
!**                                             S U B R O U T I N E   G R I D  A R E A 1                                          **
!***********************************************************************************************************************************

SUBROUTINE GRID_AREA1 (EL1,EL2,DIFF,BTOP)
  USE GLOBAL; USE GEOMC; USE PREC

  REAL(R8) :: BAREA1,BAREA2,DIST,DIST1,SLPE,DIST2
  INTEGER  :: K,K1,K2
  REAL(R8) :: EL1,EL2,DIFF,BTOP

! Difference in areas for trapezoidal geometry

  DO K=2,KB(I)
    IF (EL(K,I) <= EL1) THEN
      K1 = K
      EXIT
    END IF
  END DO
  DO K=2,KB(I)
    IF (EL(K,I) <= EL2) THEN
      K2 = K
      EXIT
    END IF
  END DO
  BAREA1 = 0.0
  BAREA2 = 0.0
  DO K=KB(I),K1,-1
    BAREA1 = BAREA1+BH(K,I)
  END DO
  DIST = EL1-EL(K1,I)
  IF (H(K1-1,JW)/2.0 < DIST) THEN
    DIST1  =  H(K1-1,JW)*0.5
    SLPE   = (B(K1-1,I)-BB(K1-1,I))/(0.5*H(K1-1,JW))
    BAREA1 =  BAREA1+BB(K1-1,I)*DIST1+0.5*SLPE*DIST1*DIST1
    DIST2  =  DIST-H(K1-1,JW)*0.5
    SLPE   = (BB(K1-2,I)-B(K1-1,I))/(0.5*H(K1-1,JW))
    BAREA1 =  BAREA1+B(K1-1,I)*DIST2+0.5*SLPE*DIST2*DIST2
    BTOP   =  B(K1-1,I)+DIST2*SLPE
  ELSE
    SLPE   = (B(K1-1,I)-BB(K1-1,I))/(0.5*H(K1-1,JW))
    BAREA1 =  BAREA1+BB(K1-1,I)*DIST+0.5*SLPE*DIST*DIST
    BTOP   =  BB(K1-1,I)+DIST*SLPE
  END IF
  DO K=KB(I),K2,-1
    BAREA2 = BAREA2+BH(K,I)
  END DO
  DIST = EL2-EL(K2,I)
  IF (H(K2-1,JW)/2. < DIST) THEN
    DIST1  =  H(K2-1,JW)*0.5
    SLPE   = (B(K2-1,I)-BB(K2-1,I))/(0.5*H(K2-1,JW))
    BAREA2 =  BAREA2+BB(K2-1,I)*DIST1+0.5*SLPE*DIST1*DIST1
    DIST2  =  DIST-H(K2-1,JW)*0.5
    SLPE   = (BB(K2-2,I)-B(K2-1,I))/(0.5*H(K2-1,JW))
    BAREA2 =  BAREA2+B(K2-1,I)*DIST2+0.5*SLPE*DIST2*DIST2
  ELSE
    SLPE   = (B(K2-1,I)-BB(K2-1,I))/(0.5*H(K2-1,JW))
    BAREA2 =  BAREA2+BB(K2-1,I)*DIST+0.5*SLPE*DIST*DIST
  END IF
  DIFF = BAREA1-BAREA2
  RETURN
END SUBROUTINE

!***********************************************************************************************************************************
!**                                             S U B R O U T I N E   G R I D  A R E A 2                                          **
!***********************************************************************************************************************************

SUBROUTINE GRID_AREA2
  USE GLOBAL; USE GEOMC; USE RSTART
  INTEGER   :: K
  REAL(R8)  :: AREA,SL,A_COEF,B_COEF,C_COEF

  AREA   = (EL(KT,I)-SZ(I)-(EL(KT,I)-Z(I)))*BI(KT,I)
  SL     = (B(KT,I)-BB(KT,I))/(0.5*H(KT,JW))
  A_COEF = -1.0
  B_COEF =  SZ(I)*2.+BI(KT,I)/(0.5*SL)
  C_COEF = -AREA/(0.5*SL)-SZ(I)**2-BI(KT,I)*2.*SZ(I)/SL
  Z(I)   = (-B_COEF+SQRT(B_COEF**2-4.*A_COEF*C_COEF))/(2.0*A_COEF)
  KTI(I) = 2
  DO K=2,KB(I)
    IF (EL(K,I) <= EL(KT,I)-Z(I)) THEN
      KTI(I) = K-1
      EXIT
    END IF
  END DO
  RETURN
END SUBROUTINE

!***********************************************************************************************************************************
!**                                              S U B R O U T I N E   L U D C M P                                                **
!***********************************************************************************************************************************

SUBROUTINE LUDCMP (A,N,NP,INDX,D)
USE PREC
  INTEGER         :: I,J,K,IMAX,NP,N
  INTEGER         :: INDX(NP)
  REAL(R8)        :: A(NP,NP),D
  REAL(R8)        :: VV(500),AAMAX,SUM,DUM
  REAL, PARAMETER :: TINY=1.0E-20


  D = 1.0
  DO I=1,N
    AAMAX = 0.0
    DO J=1,N
      IF (ABS(A(I,J)) > AAMAX) AAMAX = ABS(A(I,J))
    END DO
    VV(I) = 1.0/AAMAX
  END DO
  DO J=1,N
    DO I=1,J-1
      SUM = A(I,J)
      DO K=1,I-1
        SUM = SUM-A(I,K)*A(K,J)
      END DO
      A(I,J) = SUM
    END DO
    AAMAX = 0.0
    DO I=J,N
      SUM = A(I,J)
      DO K=1,J-1
        SUM = SUM-A(I,K)*A(K,J)
      END DO
      A(I,J) = SUM
      DUM = VV(I)*ABS(SUM)
      IF (DUM >= AAMAX) THEN
        IMAX  = I
        AAMAX = DUM
      END IF
    END DO
    IF (J /= IMAX) THEN
      DO K=1,N
        DUM       = A(IMAX,K)
        A(IMAX,K) = A(J,K)
        A(J,K)    = DUM
      END DO
      D        = -D
      VV(IMAX) =  VV(J)
    END IF
    INDX(J) = IMAX
    IF (A(J,J) == 0.0) A(J,J) = TINY
    IF (J /= N) THEN
      DUM = 1.0/A(J,J)
      DO I=J+1,N
        A(I,J) = A(I,J)*DUM
      END DO
    END IF
  END DO
END SUBROUTINE LUDCMP

!***********************************************************************************************************************************
!**                                              S U B R O U T I N E   L U B K S B                                                **
!***********************************************************************************************************************************

SUBROUTINE LUBKSB (A,N,NP,INDX,B)
USE PREC
  INTEGER :: N, NP
  REAL(R8):: A(NP,NP), B(N),SUM
  INTEGER :: INDX(NP)
  INTEGER :: I, II, J, LL

  II = 0
  DO I=1,N
    LL    = INDX(I)
    SUM   = B(LL)
    B(LL) = B(I)
    IF (II /= 0) THEN
      DO J=II,I-1
        SUM = SUM-A(I,J)*B(J)
      END DO
    ELSE IF (SUM /= 0.0) THEN
      II = I
    END IF
    B(I) = SUM
  END DO
  DO I=N,1,-1
    SUM = B(I)
    DO J=I+1,N
      SUM = SUM-A(I,J)*B(J)
    END DO
    B(I) = SUM/A(I,I)
  END DO
END SUBROUTINE LUBKSB

!***********************************************************************************************************************************
!**                                                  F U N C T I O N   B A R E A                                                  **
!***********************************************************************************************************************************

REAL (R8) FUNCTION BAREA (DEPTH,DIA)
  USE PREC
  REAL(R8), PARAMETER ::PI=3.14159265359D0
  REAL(R8) :: DEPTH,DIA
  IF (DEPTH < DIA) THEN
    BAREA = (DEPTH-DIA*0.5D0)*DSQRT(DEPTH*DIA-DEPTH**2)+(DIA**2*0.25D0)*DASIN((2.0D0/DIA)*(DEPTH-DIA*0.5D0))+(PI*DIA**2)/8.0D0
  ELSE
    BAREA = (PI*DIA**2)*0.25D0
  END IF
END FUNCTION BAREA

!***********************************************************************************************************************************
!**                                                  F U N C T I O N   T W I D T H                                                **
!***********************************************************************************************************************************

REAL (R8) FUNCTION TWIDTH (DEPTH,DIA)
USE PREC
REAL(R8) :: DEPTH,DIA
  IF (DEPTH < DIA) THEN
    TWIDTH = 2.0D0*DSQRT((DIA*DEPTH)-DEPTH**2)
  ELSE
    TWIDTH = 0.005D0*DIA
  END IF
END FUNCTION TWIDTH

!***********************************************************************************************************************************
!**                                                  F U N C T I O N   W E T P E R                                                **
!***********************************************************************************************************************************

REAL (R8) FUNCTION WETPER (DEPTH,DIA)
USE PREC
REAL(R8) :: DEPTH,DIA
  REAL(R8), PARAMETER :: PI=3.14159265359D0
  IF (DEPTH < DIA) THEN
    WETPER = DIA*(DASIN((2.0D0/DIA)*(DEPTH-DIA*0.5D0))+PI*0.5D0)
  ELSE
    WETPER = PI*DIA
  END IF
END FUNCTION WETPER

!***********************************************************************************************************************************
!**                                                F U N C T I O N   D E P T H C R I T                                            **
!***********************************************************************************************************************************

REAL (R8) FUNCTION DEPTHCRIT (FLOW)
  USE STRUCTURES
  REAL(R8) :: FLOW,X1,X2,TOL,ZBRENT1
  X1        = DIA/1.0D7
  X2        = DIA
  TOL       = 0.001
  DEPTHCRIT = ZBRENT1(X1,X2,TOL,FLOW)
END FUNCTION DEPTHCRIT

!***********************************************************************************************************************************
!**                                                  F U N C T I O N   C D F U N C                                                **
!***********************************************************************************************************************************

REAL (R8) FUNCTION CDFUNC (DEPTH,FLOW)
  USE STRUCTURES
  REAL(R8)  :: DEPTH,FLOW,BAREA,TWIDTH
  CDFUNC = (FLOW**2*TWIDTH(DEPTH,DIA))/(BAREA(DEPTH,DIA)**3*9.81D0)-1.0D0
END FUNCTION CDFUNC

!***********************************************************************************************************************************
!**                                                  F U N C T I O N   Z B R E N T                                                **
!***********************************************************************************************************************************

REAL (R8) FUNCTION ZBRENT1 (X1,X2,TOL,BARG)
USE PREC
  EXTERNAL   CDFUNC
  REAL, PARAMETER :: FACTOR=0.1,NTRY=50,ITMAX=100,EPS=3.E-8
  INTEGER         :: I,J,ITER
  REAL(R8)            :: F1,F2,X1,X2,TOL,BARG,BA,B,FA,FB,FC,CDFUNC
  REAL(R8)            :: C,D,E,TOL1,XM,S,P,Q,R

  F1 = CDFUNC(X1,BARG)
  F2 = CDFUNC(X2,BARG)
  IF (F1 <= 0.0) THEN
    DO I=1,40
      X1 = X1/10.0
      F1 = CDFUNC(X1,BARG)
      IF (F1 > 0.0) EXIT
    END DO
  END IF
  DO J=1,NTRY
    IF (F1*F2 < 0.0) EXIT
    IF (ABS(F1) < ABS(F2)) THEN
      X1 = X1+FACTOR*(X1-X2)
      F1 = CDFUNC(X1,BARG)
    ELSE
      X2 = X2+FACTOR*(X2-X1)
      F2 = CDFUNC(X2,BARG)
    END IF
  END DO
  BA = X1
  B  = X2
  FA = CDFUNC(BA,BARG)
  FB = CDFUNC(B,BARG)
  FC = FB
  DO ITER=1,ITMAX
    IF (FB*FC > 0.0) THEN
      C  = BA
      FC = FA
      D  = B-BA
      E  = D
    END IF
    IF (ABS(FC) < ABS(FB)) THEN
      BA = B
      B  = C
      C  = BA
      FA = FB
      FB = FC
      FC = FA
    END IF
    TOL1 = 2.0*EPS*ABS(B)+0.5*TOL
    XM   = 0.5*(C-B)
    IF (ABS(XM) <= TOL1 .OR. FB == 0.0) THEN
      ZBRENT1 = B; EXIT
    END IF
    IF (ABS(E) >= TOL1 .AND. ABS(FA) > ABS(FB)) THEN
      S = FB/FA
      IF (BA == C) THEN
        P = 2.0*XM*S
        Q = 1.0-S
      ELSE
        Q =  FA/FC
        R =  FB/FC
        P =  S*(2.*XM*Q*(Q-R)-(B-BA)*(R-1.0))
        Q = (Q-1.0)*(R-1.0)*(S-1.0)
      END IF
      IF (P > 0.0) Q = -Q
      P = ABS(P)
      IF (2.0*P < MIN(3.0*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
        E = D
        D = P/Q
      ELSE
        D = XM
        E = D
      END IF
    ELSE
      D = XM
      E = D
    END IF
    BA = B
    FA = FB
    IF (ABS(D) > TOL1) THEN
      B = B+D
    ELSE
      B = B+SIGN(TOL1,XM)
    END IF
    FB = CDFUNC(B,BARG)
  END DO
  ZBRENT1 = B
END FUNCTION ZBRENT1
