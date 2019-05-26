
!***********************************************************************************************************************************
!**                                            S U B R O U T I N E   K I N E T I C S                                              **
!***********************************************************************************************************************************

SUBROUTINE KINETICS
  USE SCREENC; USE GLOBAL; USE KINETIC; USE GEOMC; USE TVDC; USE LOGICC; USE SURFHE
  USE MACROPHYTEC; USE ZOOPLANKTONC; USE MAIN, ONLY:NGCTDG, EPIPHYTON_CALC, BOD_CALC, ALG_CALC, BOD_CALCN, BOD_CALCP, PO4_CALC, N_CALC, DSI_CALC, SEDCOMP_EXIST, JG_AGE, WATER_AGE_ACTIVE
  !SP CEMA
  Use CEMAVars
  !End SP CEMA

! Type declarations
  IMPLICIT NONE
  
  REAL                                :: LAM1,   LAM2,   NH4PR,  NO3PR,  LIMIT,  LIGHT,  L, L0, L1, EA, N2SAT  ! SW 10/17/15
  REAL                                :: KW,     INCR,   OH,     K1,     K2, bicart
  REAL                                :: CART,ALKT,T1K,S2,SQRS2,DH1,DH2,H2CO3T,CO3T,PHT,F,HION,HCO3T
  REAL                                :: LTCOEFM, LAVG,  MACEXT, TMAC,MACEXT1         ! CB 4/20/11
  REAL                                :: FETCH, U2, COEF1,COEF2,COEF3,COEF4,HS,TS,COEF,UORB,TAU
  REAL                                :: EPSILON, CBODSET, DOSAT,O2EX,CO2EX,SEDSI,SEDEM, SEDSO,SEDSIP
  REAL                                :: SEDSOP,SEDSON,SEDSOC,SEDSIC,SEDSIDK,SEDSUM,SEDSUMK,XDUM
  REAL                                :: BLIM, SEDSIN, COLB,COLDEP,BMASS,BMASSTEST,CVOL
  REAL                                :: ALGEX, SSEXT, TOTMAC, ZOOEXT, TOTSS0, FDPO4, ZMINFAC, SSR
  REAL                                :: ZGZTOT,CBODCT,CBODNT,CBODPT,BODTOT  ! CB 6/6/10
  REAL                                :: ALGP,ALGN,ZOOP,ZOON,TPSS,XX   ! SW 4/5/09
  ! enhanced pH buffering start
  real                                :: ammt,phost,omct,dh3,dhh,po4t,ht,hpo4t,h2po4t,oht
  real                                :: nh4t,nh3t,h3po4t,kamm,kp1,kp2,kp3
  ! enhanced pH buffering end
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: OMTRM,  SODTRM, NH4TRM, NO3TRM, BIBH2
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: DOM,    POM,    PO4BOD, NH4BOD, TICBOD
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: LAM2M  
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ATRM,   ATRMR,  ATRMF
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ETRM,   ETRMR,  ETRMF
  INTEGER                             :: K, JA, JE, M, JS, JT, JJ, JJZ, JG, JCB, JBOD, LLM,J,JD
  INTEGER                             :: MI,JAF,N,ITER,IBOD
! CEMA start
  integer                             :: jcg
  real                                :: h2sex,ch4ex,ticch4
  real                                :: sdalgc,sdepc,sdbodc,sdalgn,sdepn,sdbodn,sdalgp,sdepp,sdbodp
  logical                             :: FeMn
! CEMA end
  SAVE

! Allocation declarations

  ALLOCATE (OMTRM(KMX,IMX),    SODTRM(KMX,IMX),    NH4TRM(KMX,IMX),    NO3TRM(KMX,IMX), DOM(KMX,IMX), POM(KMX,IMX))
  ALLOCATE (PO4BOD(KMX,IMX),   NH4BOD(KMX,IMX),    TICBOD(KMX,IMX))
  ALLOCATE (ATRM(KMX,IMX,NAL), ATRMR(KMX,IMX,NAL), ATRMF(KMX,IMX,NAL))
  ALLOCATE (ETRM(KMX,IMX,NEP), ETRMR(KMX,IMX,NEP), ETRMF(KMX,IMX,NEP))
  ALLOCATE (lam2m(KMX,kmx),    BIBH2(KMX,IMX))
RETURN

!***********************************************************************************************************************************
!**                                      T E M P E R A T U R E  R A T E  M U L T I P L I E R S                                    **
!***********************************************************************************************************************************

ENTRY TEMPERATURE_RATES
  DO I=IU,ID
    DO K=KT,KB(I)
      LAM1        = FR(T1(K,I),NH4T1(JW),NH4T2(JW),NH4K1(JW),NH4K2(JW))
      NH4TRM(K,I) = LAM1/(1.0+LAM1-NH4K1(JW))
      LAM1        = FR(T1(K,I),NO3T1(JW),NO3T2(JW),NO3K1(JW),NO3K2(JW))
      NO3TRM(K,I) = LAM1/(1.0+LAM1-NO3K1(JW))
      LAM1        = FR(T1(K,I),OMT1(JW),OMT2(JW),OMK1(JW),OMK2(JW))
      OMTRM(K,I)  = LAM1/(1.0+LAM1-OMK1(JW))
      LAM1        = FR(T1(K,I),SODT1(JW),SODT2(JW),SODK1(JW),SODK2(JW))
      SODTRM(K,I) = LAM1/(1.0+LAM1-SODK1(JW))
      DO JA=1,NAL
        IF(ALG_CALC(JA))THEN
        LAM1          = FR(T1(K,I),AT1(JA),AT2(JA),AK1(JA),AK2(JA))
        LAM2          = FF(T1(K,I),AT3(JA),AT4(JA),AK3(JA),AK4(JA))
        ATRMR(K,I,JA) = LAM1/(1.0+LAM1-AK1(JA))
        ATRMF(K,I,JA) = LAM2/(1.0+LAM2-AK4(JA))
        ATRM(K,I,JA)  = ATRMR(K,I,JA)*ATRMF(K,I,JA)
        ENDIF
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
        LAM1          = FR(T1(K,I),ET1(JE),ET2(JE),EK1(JE),EK2(JE))
        LAM2          = FF(T1(K,I),ET3(JE),ET4(JE),EK3(JE),EK4(JE))
        ETRMR(K,I,JE) = LAM1/(1.0+LAM1-EK1(JE))
        ETRMF(K,I,JE) = LAM2/(1.0+LAM2-EK4(JE))
        ETRM(K,I,JE)  = ETRMR(K,I,JE)*ETRMF(K,I,JE)
        endif
      END DO
      DO M=1,NMC
      IF(MACROPHYTE_CALC(JW,M))THEN
        LAM1    = FR(T1(K,I),MT1(M),MT2(M),MK1(M),MK2(M))
        LAM2    = FF(T1(K,I),MT3(M),MT4(M),MK3(M),MK4(M))
        MACTRMR(K,I,M) = LAM1/(1.0+LAM1-MK1(M))
        MACTRMF(K,I,M) = LAM2/(1.0+LAM2-MK4(M))
        MACTRM(K,I,M)  = MACTRMR(K,I,M)*MACTRMF(K,I,M)
      endif
      end do
      IF(ZOOPLANKTON_CALC)THEN
	    DO JZ = 1, NZP
          LAM1       = FR(T1(K,I),ZT1(JZ),ZT2(JZ),ZK1(JZ),ZK2(JZ))
          LAM2       = FF(T1(K,I),ZT3(JZ),ZT4(JZ),ZK3(JZ),ZK4(JZ))
          ZOORMR(K,I,JZ)= LAM1/(1.+LAM1-ZK1(JZ))
          ZOORMF(K,I,JZ)= LAM2/(1.+LAM2-ZK4(JZ))
          ZOORM(K,I,JZ) = ZOORMR(K,I,JZ)*ZOORMF(K,I,JZ)
        END DO
	  end if
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                 K I N E T I C   R A T E S                                                     **
!***********************************************************************************************************************************

ENTRY KINETIC_RATES

! Decay rates
!!$OMP PARALLEL DO
  DO I=IU,ID
    DO K=KT,KB(I)
      DO1(K,I)          = O2(K,I)/(O2(K,I)+KDO)                  
      DO2(K,I)          = 1.0 - DO1(K,I)                         !O2(K,I)/(O2(K,I)+KDO)
      DO3(K,I)          = (1.0+SIGN(1.0,O2(K,I)-1.E-10)) *0.5
      SEDD(K,I)         =   SODTRM(K,I) *SDKV(K,I)   *SED(K,I) *DO3(K,I)   !CB 10/22/06
      SEDDP(K,I)         =  SODTRM(K,I) *SDKV(K,I)   *SEDP(K,I) *DO3(K,I)
      SEDDN(K,I)         =  SODTRM(K,I) *SDKV(K,I)   *SEDN(K,I) *DO3(K,I)
      SEDDC(K,I)         =  SODTRM(K,I) *SDKV(K,I)   *SEDC(K,I) *DO3(K,I)
! Amaila start
      IF(SEDCOMP_EXIST)THEN
      SEDD1(K,I)         =   SODTRM(K,I) *SDK1(jw)   *SED1(K,I) *DO3(K,I)
      SEDD2(K,I)         =   SODTRM(K,I) *SDK2(jw)   *SED2(K,I) *DO3(K,I)
      ENDIF
! amaila end
      SEDBR(K,I)         =  SEDB(JW)    *SED(K,I)                           !CB 11/30/06
      SEDBRP(K,I)        =  SEDB(JW)    *SEDP(K,I)                          !CB 11/30/06
      SEDBRN(K,I)        =  SEDB(JW)    *SEDN(K,I)                          !CB 11/30/06
      SEDBRC(K,I)        =  SEDB(JW)    *SEDC(K,I)                          !CB 11/30/06
      NH4D(K,I)         =  NH4TRM(K,I) *NH4DK(JW) *NH4(K,I) *DO1(K,I)
      NO3D(K,I)         =  NO3TRM(K,I) *NO3DK(JW) *NO3(K,I) *DO2(K,I)
      LDOMD(K,I)        =  OMTRM(K,I)  *LDOMDK(JW)*LDOM(K,I)*DO3(K,I)
      RDOMD(K,I)        =  OMTRM(K,I)  *RDOMDK(JW)*RDOM(K,I)*DO3(K,I)
      LPOMD(K,I)        =  OMTRM(K,I)  *LPOMDK(JW)*LPOM(K,I)*DO3(K,I)
      RPOMD(K,I)        =  OMTRM(K,I)  *RPOMDK(JW)*RPOM(K,I)*DO3(K,I)
      LRDOMD(K,I)       =  OMTRM(K,I)  *LRDDK(JW) *LDOM(K,I)*DO3(K,I)
      LRPOMD(K,I)       =  OMTRM(K,I)  *LRPDK(JW) *LPOM(K,I)*DO3(K,I)
      CBODD(K,I,1:NBOD) =  KBOD(1:NBOD)*TBOD(1:NBOD)**(T1(K,I)-20.0)*DO3(K,I)
        IF(K == KB(I))THEN     ! SW 4/18/07
	  SODD(K,I)         =  SOD(I)/BH2(K,I)*SODTRM(K,I)*BI(K,I)
	    ELSE
      SODD(K,I)         =  SOD(I)/BH2(K,I)*SODTRM(K,I)*(BI(K,I)-BI(K+1,I))
	    ENDIF

! Inorganic suspended solids settling rates - P adsorption onto SS and Fe
      FPSS(K,I) = PARTP(JW)         /(PARTP(JW)*TISS(K,I)+PARTP(JW)*FE(K,I)*DO1(K,I)+1.0)
      FPFE(K,I) = PARTP(JW)*FE(K,I) /(PARTP(JW)*TISS(K,I)+PARTP(JW)*FE(K,I)*DO1(K,I)+1.0)
      SSSI(K,I) = SSSO(K-1,I)
      TOTSS0    = 0.0
      DO JS=1,NSS
        TOTSS0 = TOTSS0+SSS(JS)*FPSS(K,I)*SS(K,I,JS)
      END DO
      SSSO(K,I) = (TOTSS0+FES(JW)*FPFE(K,I))*BI(K,I)/BH2(K,I)*DO1(K,I)                ! SW 11/7/07
      FPSS(K,I) =  FPSS(K,I)*TISS(K,I)

! OM stoichiometry
        ORGPLD(K,I)=0.0
        ORGPRD(K,I)=0.0
        ORGPLP(K,I)=0.0
        ORGPRP(K,I)=0.0
        ORGNLD(K,I)=0.0
        ORGNRD(K,I)=0.0
        ORGNLP(K,I)=0.0
        ORGNRP(K,I)=0.0
        IF(CAC(NLDOMP) == '      ON')THEN
          IF(LDOM(K,I).GT.0.0)THEN
          ORGPLD(K,I)=LDOMP(K,I)/LDOM(K,I)
          ELSE
          ORGPLD(K,I)=ORGP(JW)
          ENDIF
        ELSE
          ORGPLD(K,I)=ORGP(JW)
        END IF
        IF(CAC(NRDOMP) == '      ON')THEN
          IF(RDOM(K,I).GT.0.0)THEN
          ORGPRD(K,I)=RDOMP(K,I)/RDOM(K,I)
          ELSE
          ORGPRD(K,I)=ORGP(JW)
          ENDIF
        ELSE
          ORGPRD(K,I)=ORGP(JW)
        END IF
        IF(CAC(NLPOMP) == '      ON')THEN
          IF(LPOM(K,I).GT.0.0)THEN
          ORGPLP(K,I)=LPOMP(K,I)/LPOM(K,I)
          ELSE
          ORGPLP(K,I)=ORGP(JW)
          ENDIF
        ELSE
          ORGPLP(K,I)=ORGP(JW)
        END IF
        IF(CAC(NRPOMP) == '      ON')THEN
          IF(RPOM(K,I).GT.0.0)THEN
          ORGPRP(K,I)=RPOMP(K,I)/RPOM(K,I)
          ELSE
          ORGPRP(K,I)=ORGP(JW)
          ENDIF
        ELSE
          ORGPRP(K,I)=ORGP(JW)
        END IF
        IF(CAC(NLDOMN) == '      ON')THEN
          IF(LDOM(K,I).GT.0.0)THEN
          ORGNLD(K,I)=LDOMN(K,I)/LDOM(K,I)
          ELSE
          ORGNLD(K,I)=ORGN(JW)
          ENDIF
        ELSE
          ORGNLD(K,I)=ORGN(JW)
        END IF
        IF(CAC(NRDOMN) == '      ON')THEN
          IF(RDOM(K,I).GT.0.0)THEN
          ORGNRD(K,I)=RDOMN(K,I)/RDOM(K,I)
          ELSE
          ORGNRD(K,I)=ORGN(JW)
          ENDIF
        ELSE
          ORGNRD(K,I)=ORGN(JW)
        END IF
        IF(CAC(NLPOMN) == '      ON')THEN
          IF(LPOM(K,I).GT.0.0)THEN
          ORGNLP(K,I)=LPOMN(K,I)/LPOM(K,I)
          ELSE
          ORGNLP(K,I)=ORGN(JW)
          ENDIF
        ELSE
          ORGNLP(K,I)=ORGN(JW)
        END IF
        IF(CAC(NRPOMN) == '      ON')THEN     ! SR 8/2/2017
          IF(RPOM(K,I).GT.0.0)THEN
          ORGNRP(K,I)=RPOMN(K,I)/RPOM(K,I)
          ELSE
          ORGNRP(K,I)=ORGN(JW)
          ENDIF
        ELSE
          ORGNRP(K,I)=ORGN(JW)
        END IF

! Light Extinction Coefficient
      IF (.NOT. READ_EXTINCTION(JW)) THEN
      ALGEX = 0.0; SSEXT = 0.0; ZOOEXT = 0.0                                                     ! SW 11/8/07
        DO JA=1,NAL
          IF(ALG_CALC(JA))ALGEX = ALGEX+EXA(JA)*ALG(K,I,JA)
        END DO
        DO JS=1,NSS
          SSEXT = SSEXT+EXSS(JW)*SS(K,I,JS)
        END DO
 !       TOTMAC=0.0                                                                ! SW 4/20/11 Delete this section?
 !       DO M=1,NMC
 !         IF(MACROPHYTE_CALC(JW,M))THEN
 !           JT=KTI(I)
 !           JE=KB(I)
 !           DO JJ=JT,JE
 !             TOTMAC = EXM(M)*MACRM(JJ,K,I,M)+TOTMAC
 !           END DO
 !         END IF
 !       END DO
 !       MACEXT=TOTMAC/(BH2(K,I)*DLX(I))

	    IF(ZOOPLANKTON_CALC)THEN
	        DO JZ = 1,NZP
	        ZOOEXT = ZOOEXT + ZOO(K,I,JZ)*EXZ(JZ)
	        END DO
	    ENDIF
		
		GAMMA(K,I) = EXH2O(JW)+SSEXT+EXOM(JW)*(LPOM(K,I)+RPOM(K,I))+ALGEX+ZOOEXT         ! sw 4/21/11
		
	    IF(NMC>0)THEN    ! cb 4/20/11
	      MACEXT1=0.0    ! cb 4/20/11
          IF(KTICOL(I))THEN
            JT=KTI(I)
          ELSE
            JT=KTI(I)+1
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            TOTMAC=0.0
            DO M=1,NMC
              IF(MACROPHYTE_CALC(JW,M))THEN
                TOTMAC = EXM(M)*MACRM(JJ,K,I,M)+TOTMAC
              END IF
            END DO
            IF(CW(JJ,I).GT.0.0)THEN
              MACEXT=TOTMAC/(CW(JJ,I)*DLX(I)*H2(K,I))
            ELSE
              MACEXT=0.0
            END IF
			GAMMAJ(JJ,K,I) = GAMMA(K,I)+MACEXT       ! SW 4/20/11
            MACEXT1 = MACEXT*CW(JJ,I)+MACEXT1    ! cb 4/20/11
          END DO
          GAMMA(K,I) = GAMMA(K,I) + MACEXT1/B(JT,I)                                      ! SW 4/21/11
        end if
      ELSE
        GAMMA(K,I) = EXH2O(JW)
      END IF

! Zooplankton Rates
   IF(ZOOPLANKTON_CALC)THEN
      DO JZ=1,NZP
        TGRAZE(K,I,JZ)=PREFP(JZ)*LPOM(K,I)
        DO JJZ = 1, NZP
          TGRAZE(K,I,JZ) = TGRAZE(K,I,JZ) + PREFZ(JJZ,JZ)*ZOO(K,I,JJZ)          !CB 5/17/2007
      END DO
        DO JA=1,NAL
          IF(ALG_CALC(JA))TGRAZE(K,I,JZ)=PREFA(JA,JZ)*ALG(K,I,JA)+TGRAZE(K,I,JZ)
        END DO
        ZMINFAC  = (1.0+SIGN(1.0,ZOO(K,I,JZ)-ZOOMIN(JZ)))*0.5
        ZRT(K,I,JZ) =  ZOORMR(K,I,JZ)*ZR(JZ)*ZMINFAC*DO3(K,I)
        IF (TGRAZE(K,I,JZ) <= 0.0 .OR. O2(K,I) < 2.0) THEN
          ZMU(K,I,JZ)       = 0.0
          AGZ(K,I,1:NAL,JZ) = 0.0
		  ZGZ(K,I,JZ,:) = 0.0
          IF (O2(K,I) < 2.0) ZMINFAC = 2*ZMINFAC
        ELSE
          ZMU(K,I,JZ) = MAX(ZOORM(K,I,JZ)*ZG(JZ)*(TGRAZE(K,I,JZ)-ZOOMIN(JZ))/(TGRAZE(K,I,JZ)+ZS2P(JZ)), 0.0)
          DO JA=1,NAL
          IF(ALG_CALC(JA))AGZ(K,I,JA,JZ) = ZMU(K,I,JZ)*ZOO(K,I,JZ)*(ALG(K,I,JA)*PREFA(JA,JZ)/TGRAZE(K,I,JZ))                      !  KV 5/26/2007
          END DO
          DO JJZ = 1,NZP ! OMNIVOROUS ZOOPLANKTON
          ZGZ(K,I,JJZ,JZ)  = ZMU(K,I,JZ)*ZOO(K,I,JZ)*(ZOO(K,I,JJZ)*PREFZ(JJZ,JZ)/TGRAZE(K,I,JZ))         !KV 5/26/2007
          END DO
        END IF
        ZMT(K,I,JZ) = MAX(1.0-ZOORMF(K,I,JZ),0.02)*ZM(JZ)*ZMINFAC
      END DO   ! ZOOP LOOP
   ENDIF

    END DO ! K LOOP
  END DO   ! I LOOP
!!$OMP END PARALLEL DO

! Algal rates
   DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
      do i=iu,id
!**** Limiting factor
      LIGHT = (1.0-BETA(JW))*SRON(JW)*SHADE(I)/ASAT(JA)
      LAM1  =  LIGHT
      LAM2  =  LIGHT
      DO K=KT,KB(I)

!****** Limiting factor
        LAM1           = LAM2
        LAM2           = LAM1*EXP(-GAMMA(K,I)*H2(K,I))
        FDPO4          = 1.0-FPSS(K,I)-FPFE(K,I)
        ALLIM(K,I,JA)  = 2.718282*(EXP(-LAM2)-EXP(-LAM1))/(GAMMA(K,I)*H2(K,I))
        IF (AHSP(JA)  /= 0.0 .and. po4_calc) APLIM(K,I,JA) =  FDPO4*PO4(K,I)/(FDPO4*PO4(K,I)+AHSP(JA)+NONZERO)       ! cb 10/12/11
        IF (AHSN(JA)  /= 0.0 .and. n_calc) ANLIM(K,I,JA) = (NH4(K,I)+NO3(K,I))/(NH4(K,I)+NO3(K,I)+AHSN(JA)+NONZERO)  ! cb 10/12/11
        IF (AHSSI(JA) /= 0.0 .and. DSI_CALC) ASLIM(K,I,JA) =  DSI(K,I)/(DSI(K,I)+AHSSI(JA)+NONZERO)                  ! cb 10/12/11
        LIMIT          = MIN(APLIM(K,I,JA),ANLIM(K,I,JA),ASLIM(K,I,JA),ALLIM(K,I,JA))

!****** Algal rates
        AGR(K,I,JA) =  ATRM(K,I,JA)*AG(JA)*LIMIT
        ARR(K,I,JA) =  ATRM(K,I,JA)*AR(JA)*DO3(K,I)
        AMR(K,I,JA) = (ATRMR(K,I,JA)+1.0-ATRMF(K,I,JA))*AM(JA)
        AER(K,I,JA) =  MIN((1.0-ALLIM(K,I,JA))*AE(JA)*ATRM(K,I,JA),AGR(K,I,JA))
        IF (AS(JA) >= 0.0) THEN
          IF(K == KT)THEN
          ASR(K,I,JA) =  AS(JA)*(-ALG(K,I,JA))*BI(K,I)/BH2(K,I)
          ELSE
          ASR(K,I,JA) =  AS(JA)*(ALG(K-1,I,JA)-ALG(K,I,JA))*BI(K,I)/BH2(K,I)
          ENDIF
        ELSE
          IF(K == KB(I))THEN
            ASR(K,I,JA) = -AS(JA)*(-ALG(K,I,JA)  *BI(K,I)/BH2(K,I))                                           !SW 11/8/07
          ELSEIF(K == KT)THEN
            ASR(K,I,JA) = -AS(JA)* ALG(K+1,I,JA)*BI(K+1,I)*DLX(I)/VOL(K,I)                                   !SW 11/8/07
          ELSE
            ASR(K,I,JA) = -AS(JA)*(ALG(K+1,I,JA)*BI(K+1,I)/BH2(K,I)-ALG(K,I,JA)*BI(K,I)/BH2(K,I))             !SP 8/27/07
          END IF
        END IF
      end do
    end do
    ENDIF
  END DO    ! ALGAE LOOP

! Macrophyte Light/Nutrient Limitation and kinetic rates
  do m=1,nmc
  mGR(:,:,iu:id,m)=0.0; mRR(:,iu:id,m)=0.0; mmR(:,iu:id,m)=0.0  ! cb 3/8/16
  if(macrophyte_calc(jw,m))then
    DO I=IU,ID
      LTCOEFm = (1.0-BETA(jw))*SRON(jw)*SHADE(I)
      if(kticol(i))then
        jt=kti(i)
      else
        jt=kti(i)+1
      end if
      je=kb(i)
      do jj=jt,je
        lam1=ltcoefm
        lam2m(jj,kt)=lam1*exp(-gammaj(jj,kt,i)*h2(kt,i))
        lavg=(lam1-lam2m(jj,kt))/(GAMMAj(jj,kt,i)*H2(kt,i))
        mLLIM(jj,kt,I,m) = lavg/(lavg+msat(m))
        IF (mHSP(m)  /= 0.0.and.psed(m) < 1.0)then
          mPLIM(kt,I,m) =  FDPO4*PO4(kt,I)/(FDPO4*PO4(kt,I)+mHSP(m)+nonzero)
        else
          mPLIM(kt,I,m)=1.0
        end if
        IF (mHSN(m)  /= 0.0.and.nsed(m) < 1.0)then
          mNLIM(kt,I,m) = NH4(kt,I)/(NH4(kt,I)+mHSN(m)+nonzero)
        else
          mNLIM(kt,I,m)=1.0
        end if
        IF (mHSc(m) /= 0.0)then
          mcLIM(kt,i,m) = co2(kt,I)/(co2(kt,I)+mHSc(m)+NONZERO)
        end if
        LIMIT          = MIN(mPLIM(kt,I,m),mNLIM(kt,I,m),mcLIM(kt,I,m),mLLIM(jj,kt,I,m))

!************* sources/sinks

        mGR(jj,Kt,I,m) = macTRM(Kt,I,m)*mG(m)*LIMIT

      end do

      mRR(Kt,I,m) = macTRM(Kt,I,m)*mR(m)*DO3(Kt,I)
      mMR(Kt,I,m) = (macTRMR(Kt,I,m)+1.0-mAcTRMF(Kt,I,m))*mM(m)

      DO K=KT+1,KB(I)
        jt=k
        je=kb(i)
        do jj=jt,je
          lam1=lam2m(jj,k-1)
          lam2m(jj,k)=lam1*exp(-gammaj(jj,k,i)*h2(k,i))
          lavg=(lam1-lam2m(jj,k))/(GAMMAj(jj,k,i)*H2(k,i))
          mLLIM(jj,K,I,m) = lavg/(lavg+msat(m))
          IF (mHSP(m)  /= 0.0.and.psed(m) < 1.0)then
            mPLIM(K,I,m) =  FDPO4*PO4(K,I)/(FDPO4*PO4(K,I)+mHSP(m)+nonzero)
          else
            mPLIM(K,I,m)=1.0
          end if
          IF (mHSN(m)  /= 0.0.and.nsed(m) < 1.0)then
            mNLIM(K,I,m) = NH4(K,I)/(NH4(K,I)+mHSN(m)+nonzero)
          else
             mNLIM(K,I,m)=1.0
          end if
          IF (mHSc(m) /= 0.0)then
            mcLIM(k,i,m) = co2(K,I)/(co2(K,I)+mHSc(m)+NONZERO)
          end if
          LIMIT          = MIN(mPLIM(K,I,m),mNLIM(K,I,m),mcLIM(K,I,m),mLLIM(jj,K,I,m))

!************* sources/sinks

          mGR(jj,K,I,m) = macTRM(K,I,m)*mG(m)*LIMIT

        end do

        mRR(K,I,m) = macTRM(K,I,m)*mR(m)*DO3(K,I)
        mMR(K,I,m) = (macTRMR(K,I,m)+1.0-mAcTRMF(K,I,m))*mM(m)
      end do
    END DO
    ENDIF
  END DO

RETURN

!***********************************************************************************************************************************
!**                                             G E N E R I C   C O N S T I T U E N T                                             **
!***********************************************************************************************************************************

ENTRY GENERIC_CONST (JG)

IF(WATER_AGE_ACTIVE .AND. JG==JG_AGE)THEN    ! SW 7/27/2017  Speed of computation
    DO I=IU,ID
      DO K=KT,KB(I)
               CGSS(K,I,JG) =-CG0DK(JG)
      ENDDO
    ENDDO
ELSE
    
XX=0.0
FeMn= .false.
IF(SD_global)THEN
sdinFeOOH(:,iu:id)=0.0
sdinMnO2(:,iu:id)=0.0
if(jg == ngFe2 .or. jg==ngFeOOH.or.jg == ngMn2 .or. jg==ngMnO2)FeMn=.true.
if(jg == ngh2s)then
  h2sd(:,IU:ID) = 0.0;h2sreaer(:,IU:ID) = 0.0
end if
if(jg == ngch4)then
  ch4d(:,IU:ID) = 0.0;ch4reaer(:,IU:ID) = 0.0
end if
if(jg == ngfe2) fe2d(:,IU:ID) = 0.0
if(jg == ngMn2) Mn2d(:,IU:ID) = 0.0
ENDIF
! CEMA end
DO I=IU,ID
LIGHT=(1.0-BETA(JW))*SRON(JW)*SHADE(I)                  !LCJ 2/26/15
LAM1  =  LIGHT
LAM2  =  LIGHT

      DO K=KT,KB(I)
        LAM1           = LAM2
        LAM2           = LAM1*EXP(-GAMMA(K,I)*H2(K,I))
        LIGHT          = LAM1*(1.-EXP(-GAMMA(K,I)*H2(K,I)))/(GAMMA(K,I)*H2(K,I))     ! SW 10/17/15

! CEMA start
        if(FeMn)then  
          if(jg == ngFeOOH)then
            IF(K == KT)THEN
              xx =  FeSetVel*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I) 
            ELSE
              xx =  FeSetVel*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I) 
            ENDIF
            sdinFeOOH(k,i)=xx
          end if
          if(jg == ngMnO2)then
            IF(K == KT)THEN
              xx =  MnSetVel*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I) 
            ELSE
              xx =  MnSetVel*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I) 
            ENDIF
            sdinMnO2(k,i)=xx
          end if
        else
! CEMA end
         IF (CGS(JG) > 0.0) THEN
          IF(K == KT)THEN
          xx =  CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I)    ! AS(JA)*(-ALG(K,I,JA))*BI(K,I)/BH2(K,I)
          ELSE
          xx =  CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I)     !AS(JA)*(ALG(K-1,I,JA)-ALG(K,I,JA))*BI(K,I)/BH2(K,I)
          ENDIF
         ELSEif(cgs(jg)<0.0)then
          IF(K == KB(I))THEN
            xx = -CGS(JG)*(-CG(K,I,JG))*BI(K,I)/BH2(K,I)    !-AS(JA)*(-ALG(K,I,JA)  *BI(K,I)/BH2(K,I))                                           !SW 11/8/07
          ELSEIF(K == KT)THEN
            xx = -CGS(JG)*CG(K+1,I,JG)*BI(K+1,I)*DLX(I)/VOL(K,I)    !-AS(JA)* ALG(K+1,I,JA)*BI(K+1,I)*DLX(I)/VOL(K,I)                                   !SW 11/8/07
          ELSE
            xx = -CGS(JG)*(CG(K+1,I,JG)*BI(K+1,I)/BH2(K,I)-CG(K,I,JG)*BI(K,I)/BH2(K,I))    !-AS(JA)*(ALG(K+1,I,JA)*BI(K+1,I)/BH2(K,I)-ALG(K,I,JA)*BI(K,I)/BH2(K,I))             !SP 8/27/07
          END IF
         ENDIF        
        end if  ! CEMA

        
! CEMA start
        if(FeMn)then  
          if(jg == ngFeOOH)then
            CGSS(K,I,JG) = kFe_oxid*O2(k,i)*10**(2.0*(pH(k,i)-7.0))*CG(K,I,ngFe2) - kFe_red*(KFeOOH_HalfSat/(o2(k,i)+KFeOOH_HalfSat))*cg(k,i,jg) + xx
          end if
          if(jg == ngFe2)then
            fe2d(k,i)=-kFe_oxid*O2(k,i)*10**(2.0*(pH(k,i)-7.0))*CG(K,I,JG)
            CGSS(K,I,JG) = fe2d(k,i) + kFe_red*(KFeOOH_HalfSat/(o2(k,i)+KFeOOH_HalfSat))*cg(k,i,ngFeOOH) + xx
          end if
          if(jg == ngMnO2)then
            CGSS(K,I,JG) = kMn_oxid*O2(k,i)*10**(2.0*(pH(k,i)-7.0))*CG(K,I,ngMn2) - kMn_red*(KMnO2_HalfSat/(o2(k,i)+KMnO2_HalfSat))*cg(k,i,jg) + xx
          end if
          if(jg == ngMn2)then
            Mn2d(k,i)=-kMn_oxid*O2(k,i)*10**(2.0*(pH(k,i)-7.0))*CG(K,I,JG)
            CGSS(K,I,JG) = Mn2d(k,i) + kMn_red*(KMnO2_HalfSat/(o2(k,i)+KMnO2_HalfSat))*cg(k,i,ngMnO2) + xx
          end if
        else
! CEMA end

         IF (CGQ10(JG) /= 0.0) THEN
!             CGSS(K,I,JG) = -CG0DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)-CG1DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)*CG(K,I,JG)+xx            ! SW 4/5/09 CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I)
! CEMA start
             if(jg == ngh2s)then
               h2sd(k,i)=(-CG0DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)-CG1DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)*CG(K,I,JG))  *DO3(K,I)                              
               CGSS(K,I,JG) = h2sd(k,i)
             else if(jg == ngch4)then
               ch4d(k,i)=(-CG0DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)-CG1DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)*CG(K,I,JG))  *DO3(K,I)
               CGSS(K,I,JG) = ch4d(k,i)
             else
               CGSS(K,I,JG) = -CG0DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)-CG1DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)*CG(K,I,JG)-CGLDK(JG)*LIGHT*CG(K,I,JG)+xx
             end if
! CEMA end
         ELSE
!             CGSS(K,I,JG) = -CG0DK(JG)-CG1DK(JG)*CG(K,I,JG)+xx                                                                ! SW 4/5/09 CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I)  
! CEMA start
             if(jg == ngh2s)then
               h2sd(k,i)=-CG0DK(JG)-CG1DK(JG)*CG(K,I,JG) *DO3(K,I) 
               CGSS(K,I,JG) = h2sd(k,i)+xx
             else if(jg == ngch4)then
               ch4d(k,i)=-CG0DK(JG)-CG1DK(JG)*CG(K,I,JG) *DO3(K,I)
               CGSS(K,I,JG) = ch4d(k,i)+xx
             else
               CGSS(K,I,JG) = -CG0DK(JG)-CG1DK(JG)*CG(K,I,JG)-CGLDK(JG)*LIGHT*CG(K,I,JG)+xx
             end if
! CEMA end
         ENDIF
        end if  ! CEMA
        if(jg == ngso4)cgSS(K,I,jg) = cgSS(K,I,jg) - h2sd(k,i)   ! sulfate production from sulfide decay, negative sign because 'h2sd' is negative
     END DO
! CEMA start
     
     if(jg == ngh2s)then
       IF (.NOT. ICE(I)) THEN
         IF (REAER(I) == 0.0) CALL GAS_TRANSFER
         H2SEX       = REAER(I)*0.984
         H2SREAER(KT,I)=h2sEX*(-cg(KT,I,jg))*BI(KT,I)/BH2(KT,I)
         cgSS(KT,I,jg) = cgSS(KT,I,jg)+h2sREAER(KT,I)
       END IF       
     end if
     if(jg == ngch4)then
       IF (.NOT. ICE(I)) THEN
         IF (REAER(I) == 0.0) CALL GAS_TRANSFER
         ch4EX       = REAER(I)*1.188
         ch4REAER(KT,I)=ch4EX*(-cg(KT,I,jg))*BI(KT,I)/BH2(KT,I)
         cgSS(KT,I,jg) = cgSS(KT,I,jg)+ch4REAER(KT,I)
       END IF       
     end if
     IF(CGKLF(JG) /= 0.0)THEN
         IF (.NOT. ICE(I)) THEN
         IF (REAER(I) == 0.0) CALL GAS_TRANSFER
         IF(CGCS(JG) == -1.0)THEN     ! THIS IS FOR N2 GAS
             EA = DEXP(2.3026D0*(7.5D0*TDEW(JW)/(TDEW(JW)+237.3D0)+0.6609D0))*0.001316   ! in mm Hg   0.0098692atm=7.5006151mmHg   
            ! PN2=0.79*(PALT(I)-EA)   ! atm with water vapor correction since 0.79 atm is for dry air
             N2SAT=1.5568D06*0.79*(PALT(I)-EA)*(1.8816D-5 - 4.116D-7 * T1(KT,I) + 4.6D-9 * T1(KT,I)**2)
     ! N2SS(KT,I) = (N2SAT-N2(KT,I))*REAER(I)*1.304*BI(KT,I)/BH2(KT,I)
         cgSS(KT,I,jg) = cgSS(KT,I,jg)+REAER(I)*CGKLF(JG)*(N2SAT-cg(KT,I,jg))*BI(KT,I)/BH2(KT,I)      ! fixed value of KLN2=1.034*KLO2
                  DO K=KT,KB(I)                  ! NOTE THERE IS 1 TIME STEP LAG WITH TDG SINCE PLACED HERE
                  DOSAT = SATO(T1(K,I),TDS(K,I),PALT(I),SALT_WATER(JW))
                  TDG(K,I)=100.*((0.79*CG(K,I,NGCTDG)/N2SAT)+O2(K,I)/DOSAT*0.21)    !(1.5568D06*0.79*(PALT(I)-EA))*(1.8816D-5 - 4.116D-7 * T1(KT,I) + 4.6D-9 * T1(KT,I)**2))
                  ENDDO
                  
         ELSE
         
         cgSS(KT,I,jg) = cgSS(KT,I,jg)+REAER(I)*CGKLF(JG)*(CGCS(JG)-cg(KT,I,jg))*BI(KT,I)/BH2(KT,I)
        END IF    
         ENDIF    
     ENDIF
     
END DO
ENDIF

 RETURN  



!  IF (CGQ10(JG) /= 0.0) THEN
!    DO I=IU,ID
!      DO K=KT,KB(I)
!        CGSS(K,I,JG) = -CG0DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)-CG1DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)*CG(K,I,JG)+xx            ! SW 4/5/09 CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I)
!      END DO
!    END DO
!  ELSE
!    DO I=IU,ID
!      DO K=KT,KB(I)
!        CGSS(K,I,JG) = -CG0DK(JG)-CG1DK(JG)*CG(K,I,JG)+                                                                   ! SW 4/5/09 CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I)
!      END DO
!    END DO
!  END IF
!RETURN

!***********************************************************************************************************************************
!**                                               S U S P E N D E D   S O L I D S                                                 **
!***********************************************************************************************************************************

ENTRY SUSPENDED_SOLIDS (J)

    !SP CEMA
    if(sediment_diagenesis)then
    If(IncludeBedConsolidation)Then
        !All resuspension done in CEMA code
        SEDIMENT_RESUSPENSION(J) = .FALSE.
    End If
    end if
    !End SP CEMA

  DO I=IU,ID
    SSR = 0.0
    IF (SEDIMENT_RESUSPENSION(J)) THEN
      FETCH = FETCHD(I,JB)
      IF (COS(PHI(JW)-PHI0(I)) < 0.0) FETCH = FETCHU(I,JB)
      FETCH = MAX(FETCH,BI(KT,I),DLX(I))
      U2    = WIND(JW)*WSC(I)*WIND(JW)*WSC(I)+NONZERO
      COEF1 = 0.53  *(G*DEPTHB(KT,I)/U2)**0.75
      COEF2 = 0.0125*(G*FETCH/U2)**0.42
      COEF3 = 0.833* (G*DEPTHB(KT,I)/U2)**0.375
      COEF4 = 0.077* (G*FETCH/U2)**0.25
      HS    = 0.283 *U2/G*0.283*TANH(COEF1)*TANH(COEF2/TANH(COEF1))
      !TS    = 2.0*PI*U2/G*1.2*  TANH(COEF3)*TANH(COEF4/TANH(COEF3))
      TS    = 2.0*PI*sqrt(U2)/G*1.2*  TANH(COEF3)*TANH(COEF4/TANH(COEF3))   ! cb 7/15/14
      L0    = G*TS*TS/(2.0*PI)
        L1 = L0             ! SW 6/28/2018 Allow for resuspension of surface layer
        L  = L0*TANH(2.0*PI*DEPTHB(KT,I)/L1)
        DO WHILE (ABS(L-L1) > 0.001)
          L1 = L
          L  = L0*TANH(2.0*PI*DEPTHB(KT,I)/L1)
        END DO
        COEF = MIN(710.0,2.0*PI*DEPTHB(KT,I)/L)
        UORB = PI*HS/TS*100.0/SINH(COEF)
        TAU  = 0.003*UORB*UORB
        IF (TAU-TAUCR(J) > 0.0) EPSILON = MAX(0.0,0.008/49.0*(TAU-TAUCR(J))**3*10000.0/DLT)
        SSR = EPSILON*DLX(I)*(BI(KT,I)-BI(KT+1,I))/VOL(KT,I)
      END IF
    SSSS(KT,I,J) = -SSS(J)*SS(KT,I,J)*BI(KT,I)/BH2(KT,I)+SSR
 !   DO K=KT-1,KB(I)-1                                             ! SW 4/3/09   KT,KB
     DO K=KT+1,KB(I)-1                 ! cb 9/29/14
      IF (SEDIMENT_RESUSPENSION(J)) THEN
        L1 = L0
        L  = L0*TANH(2.0*PI*DEPTHB(K,I)/L1)
        DO WHILE (ABS(L-L1) > 0.001)
          L1 = L
          L  = L0*TANH(2.0*PI*DEPTHB(K,I)/L1)
        END DO
        COEF = MIN(710.0,2.0*PI*DEPTHB(K,I)/L)
        UORB = PI*HS/TS*100.0/SINH(COEF)
        TAU  = 0.003*UORB*UORB
        IF (TAU-TAUCR(J) > 0.0) EPSILON = MAX(0.0,0.008/49.0*(TAU-TAUCR(J))**3*10000.0/DLT)
		if(k == kb(i))then   ! SW 4/18/07
		SSR = EPSILON*DLX(I)*BI(K,I)/VOL(K,I)
		else
        SSR = EPSILON*DLX(I)*(BI(K,I)-BI(K+1,I))/VOL(K,I)
		endif
      END IF
      SSSS(K,I,J) = SSS(J)*(SS(K-1,I,J)-SS(K,I,J))*BI(K,I)/BH2(K,I)+SSR
    END DO
    IF (SEDIMENT_RESUSPENSION(J)) SSR = EPSILON*DLX(I)*BI(KB(I),I)/VOL(KB(I),I)
!CEMA SP
     if(sediment_diagenesis)then
    If(IncludeFFTLayer .and. FFTActive)Then
        SSSS(KB(I),I,J) = (SSS(J)*SS(KB(I)-1,I,J)-FFTLayerSettVel*SS(KB(I),I,J))/H(KB(I),JW)+SSR
    End If
    If(IncludeFFTLayer .and. .NOT. FFTActive)Then
        SSSS(KB(I),I,J) = 0.d0
    End If
    If(.NOT. IncludeFFTLayer)Then
       SSSS(KB(I),I,J) = SSS(J)*(SS(KB(I)-1,I,J)-SS(KB(I),I,J))/H(KB(I),JW)+SSR
  ! Flocculation              !SR                                                      !New section on flocculation          !SR 04/21/13
    !DO K=KT,KB(I)
    !  SSF = 0.0
    !  IF (J > 1 .AND. SSFLOC(J-1) > 0.0) THEN
    !    IF (FLOCEQN(J-1) == 0) THEN
    !      SSF = MIN(SSFLOC(J-1), SS(K,I,J-1)/DLT)
    !    ELSE IF (FLOCEQN(J-1) == 1) THEN
    !      SSF = SSFLOC(J-1)*SS(K,I,J-1)
    !    ELSE IF (FLOCEQN(J-1) == 2) THEN
    !      SSF = SSFLOC(J-1)*SS(K,I,J-1)*SS(K,I,J-1)
    !    END IF
    !  END IF
    !  IF (J < NSS .AND. SSFLOC(J) > 0.0) THEN
    !    IF (FLOCEQN(J) == 0) THEN
    !      SSF = SSF - MIN(SSFLOC(J), SS(K,I,J)/DLT)
    !    ELSE IF (FLOCEQN(J) == 1) THEN
    !      SSF = SSF - SSFLOC(J)*SS(K,I,J)
    !    ELSE IF (FLOCEQN(J) == 2) THEN
    !      SSF = SSF - SSFLOC(J)*SS(K,I,J)*SS(K,I,J)
    !    END IF
    !  END IF
    !  SSSS(K,I,J) = SSSS(K,I,J) + SSF
    !END DO                                                                        !End new section on flocculation      !SR 04/21/13
    End If
     end if
    !End CEMA SP
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      P H O S P H O R U S                                                      **
!***********************************************************************************************************************************

ENTRY PHOSPHORUS
  PO4AR(:,IU:ID) = 0.0; PO4AG(:,IU:ID) = 0.0; PO4ER(:,IU:ID) = 0.0; PO4EG(:,IU:ID) = 0.0; PO4BOD(:,IU:ID) = 0.0
  PO4MR(:,IU:ID) = 0.0; PO4MG(:,IU:ID) = 0.0; PO4ZR(:,IU:ID)=0.0   

  DO I=IU,ID
    DO K=KT,KB(I)
      DO JCB=1,NBOD
!        IF(BOD_CALC(JCB))PO4BOD(K,I) = PO4BOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODP(JCB)
         IF(BOD_CALCp(JCB))then                                                ! cb 5/19/11
           PO4BOD(K,I) = PO4BOD(K,I)+CBODD(K,I,JCB)*CBODp(K,I,JCB)    
         else
           PO4BOD(K,I) = PO4BOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODP(JCB)
         end if
      END DO
      DO JA=1,NAL
        IF(ALG_CALC(JA))THEN
        PO4AG(K,I) = PO4AG(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*AP(JA)
        PO4AR(K,I) = PO4AR(K,I)+ARR(K,I,JA)*ALG(K,I,JA)*AP(JA)
        ENDIF
      END DO
      DO JE=1,NEP
      IF (EPIPHYTON_CALC(JW,JE))then
        PO4EG(K,I) = PO4EG(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*EP(JE)
        PO4ER(K,I) = PO4ER(K,I)+ERR(K,I,JE)*EPC(K,I,JE)*EP(JE)
      endif
      END DO
      PO4EP(K,I)  = PO4ER(K,I)-PO4EG(K,I)
      PO4AP(K,I)  = PO4AR(K,I)-PO4AG(K,I)
      PO4POM(K,I) = ORGPLP(k,i)*LPOMD(K,I)+orgprp(k,i)*RPOMD(K,I)
      PO4DOM(K,I) = ORGPLD(k,i)*LDOMD(K,I)+orgprd(k,i)*RDOMD(K,I)
      PO4OM(K,I)  = PO4POM(K,I)+PO4DOM(K,I)
            IF(SEDCOMP_EXIST)THEN  ! SW 5/26/15
           ! PO4SD(K,I)  = SEDDp(K,I)+sedd1(k,i)*orgp(jw) + sedd2(k,i)*orgp(jw)   ! Amaila
            PO4SD(K,I)  = SEDDp(K,I)+sedd1(k,i)*pbiom(jw) + sedd2(k,i)*pbiom(jw)   ! Amaila, cb 6/7/17
            ELSE
            PO4SD(K,I)  = SEDDp(K,I)
            ENDIF
      PO4SR(K,I)  = PO4R(JW)*SODD(K,I)*DO2(K,I)
      PO4NS(K,I)  = SSSI(K,I)*PO4(K-1,I)-SSSO(K,I)*PO4(K,I)

      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            PO4MG(K,I)= PO4MG(K,I)+MGR(JJ,K,I,M)*MACRM(JJ,K,I,M)*MP(M)*(1.0-PSED(M))
            PO4MR(K,I)= PO4MR(K,I)+MRR(K,I,M)*MACRM(JJ,K,I,M)*MP(M)
          END DO
        END IF
      END DO
      PO4MR(K,I)=PO4MR(K,I)/(DLX(I)*BH(K,I))
      PO4MG(K,I)=PO4MG(K,I)/(DLX(I)*BH(K,I))
      IF(ZOOPLANKTON_CALC)THEN
      DO JZ = 1,NZP
        PO4ZR(K,I) = PO4ZR(K,I) + ZRT(K,I,JZ)*ZOO(K,I,JZ)*ZP(JZ)
	  END DO
	  ENDIF


      PO4SS(K,I)  = PO4AP(K,I)+PO4EP(K,I)+PO4OM(K,I)+PO4SD(K,I)+PO4SR(K,I)+PO4NS(K,I)+PO4BOD(K,I)  &
                    +PO4MR(K,I)-PO4MG(K,I) +PO4ZR(K,I)    

    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                        A M M O N I U M                                                        **
!***********************************************************************************************************************************

ENTRY AMMONIUM
  NH4AG(:,IU:ID) = 0.0; NH4AR(:,IU:ID) = 0.0; NH4ER(:,IU:ID) = 0.0; NH4EG(:,IU:ID) = 0.0; NH4BOD(:,IU:ID) = 0.0
  NH4MG(:,IU:ID) = 0.0; NH4MR(:,IU:ID) = 0.0; NH4ZR(:,IU:ID)=0.0   
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JCB=1,NBOD
!        IF(BOD_CALC(JCB))NH4BOD(K,I) =  NH4BOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODN(JCB)
         IF(BOD_CALCn(JCB))then                                                ! cb 5/19/11
           NH4BOD(K,I) =  NH4BOD(K,I)+CBODD(K,I,JCB)*CBODn(K,I,JCB)
         else
           NH4BOD(K,I) =  NH4BOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODN(JCB)
         end if
      END DO
      DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
        IF (ANEQN(JA).EQ.2) THEN
        NH4PR      = NH4(K,I)*NO3(K,I)/((ANPR(JA)+NH4(K,I))*(ANPR(JA)+NO3(K,I)))+NH4(K,I)*ANPR(JA)/((NO3(K,I)  &
                                        +NH4(K,I)+NONZERO)*(ANPR(JA)+NO3(K,I)))
        ELSE
        NH4PR = NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO)
        ENDIF
        IF (AHSN(JA) > 0.0) NH4AG(K,I) = NH4AG(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*AN(JA)*NH4PR
        NH4AR(K,I) = NH4AR(K,I)+ARR(K,I,JA)*ALG(K,I,JA)*AN(JA)
      ENDIF
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
        IF (ENEQN(JE) == 2)THEN
        NH4PR = NH4(K,I)*NO3(K,I)/((ENPR(JE)+NH4(K,I))*(ENPR(JE)+NO3(K,I)))+NH4(K,I)*ENPR(JE)/((NO3(K,I)  &
                                        +NH4(K,I)+NONZERO)*(ENPR(JE)+NO3(K,I)))
        ELSE
        NH4PR = NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO)
        ENDIF
        IF (EHSN(JE) > 0.0) NH4EG(K,I) = NH4EG(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*EN(JE)*NH4PR
        NH4ER(K,I) = NH4ER(K,I)+ERR(K,I,JE)*EPC(K,I,JE)*EN(JE)
        endif
      END DO
      NH4EP(K,I)  =  NH4ER(K,I) -NH4EG(K,I)
      NH4AP(K,I)  =  NH4AR(K,I) -NH4AG(K,I)

      NH4DOM(K,I) = LDOMD(K,I)*orgnld(k,i) +RDOMD(K,I)*ORGNrd(k,i)
      NH4POM(K,I) = LPOMD(K,I)*orgnlp(k,i) +RPOMD(K,I)*ORGNrp(k,i)

      NH4OM(K,I)  =  NH4DOM(K,I)+NH4POM(K,I)

            IF(SEDCOMP_EXIST)THEN  ! SW 5/26/15
            !NH4SD(K,I)  =  SEDDn(K,I) +sedd1(k,i)*orgn(jw) + sedd2(k,i)*orgn(jw)   ! Amaila
            NH4SD(K,I)  =  SEDDn(K,I) +sedd1(k,i)*nbiom(jw) + sedd2(k,i)*nbiom(jw)   ! Amaila, cb 6/7/17
            ELSE
            NH4SD(K,I)  =  SEDDn(K,I)
            ENDIF

      NH4SR(K,I)  =  NH4R(JW) *SODD(K,I)*DO2(K,I)

      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            NH4MR(K,I)= NH4MR(K,I)+MRR(K,I,M)*MACRM(JJ,K,I,M)*MN(M)
            NH4MG(K,I)= NH4MG(K,I)+MGR(JJ,K,I,M)*MACRM(JJ,K,I,M)*MN(M)*(1.0-NSED(M))
          END DO
        END IF
      END DO
      NH4MR(K,I)=NH4MR(K,I)/(DLX(I)*BH(K,I))
      NH4MG(K,I)=NH4MG(K,I)/(DLX(I)*BH(K,I))
	  IF(ZOOPLANKTON_CALC)THEN
	  DO JZ = 1,NZP
	    NH4ZR(K,I) = NH4ZR(K,I) + ZRT(K,I,JZ)*ZOO(K,I,JZ)*ZN(JZ) 
	  END DO
	  ENDIF
      NH4SS(K,I)  =  NH4AP(K,I)+NH4EP(K,I)+NH4OM(K,I)+NH4SD(K,I)+NH4SR(K,I)+NH4BOD(K,I)-NH4D(K,I)  &
         +NH4MR(K,I)-NH4MG(K,I) +NH4ZR(K,I)     
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                          N I T R A T E                                                        **
!***********************************************************************************************************************************

ENTRY NITRATE
  NO3AG(:,IU:ID) = 0.0; NO3EG(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
        NO3PR = 1.0-NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO)
        IF (ANEQN(JA).EQ.2)  NO3PR      = 1.0-(NH4(K,I)*NO3(K,I)/((ANPR(JA)+NH4(K,I))*(ANPR(JA)+NO3(K,I)))+NH4(K,I)*ANPR(JA)       &
                                          /((NO3(K,I)+NH4(K,I)+NONZERO)*(ANPR(JA)+NO3(K,I))))
        IF (AHSN(JA).GT.0.0) NO3AG(K,I) = NO3AG(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*NO3PR*AN(JA)
      ENDIF
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
        NO3PR = 1.0-NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO)
        IF (ENEQN(JE).EQ.2)  NO3PR      = 1.0-(NH4(K,I)*NO3(K,I)/((ENPR(JE)+NH4(K,I))*(ENPR(JE)+NO3(K,I)))+NH4(K,I)*ENPR(JE)       &
                                          /((NO3(K,I)+NH4(K,I)+NONZERO)*(ENPR(JE)+NO3(K,I))))
        IF (EHSN(JE).GT.0.0) NO3EG(K,I) = NO3EG(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*NO3PR*EN(JE)
        ENDIF
      END DO
      IF(K == KB(I)) THEN      ! SW 4/18/07
      NO3SED(K,I) = NO3(K,I)*NO3S(JW)*NO3TRM(K,I)*(BI(K,I))/BH2(K,I)
	  ELSE
      NO3SED(K,I) = NO3(K,I)*NO3S(JW)*NO3TRM(K,I)*(BI(K,I)-BI(K+1,I))/BH2(K,I)
	  ENDIF
      NO3SS(K,I)  = NH4D(K,I)-NO3D(K,I)-NO3AG(K,I)-NO3EG(K,I)-NO3SED(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  D I S S O L V E D   S I L I C A                                              **
!***********************************************************************************************************************************

ENTRY DISSOLVED_SILICA
  DSIAG(:,IU:ID) = 0.0; DSIEG(:,IU:ID) = 0.0                          !; DSIBOD = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
        DSIAG(K,I) = DSIAG(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*ASI(JA)
      ENDIF
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))DSIEG(K,I) = DSIEG(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*ESI(JE)
      END DO
      DSID(K,I)  =  PSIDK(JW)*PSI(K,I)
      DSISD(K,I) =  SEDD(K,I)*ORGSI(JW)
      DSISR(K,I) =  DSIR(JW)*SODD(K,I)*DO2(K,I)
      DSIS(K,I)  = (SSSI(K,I)*DSI(K-1,I)-SSSO(K,I)*DSI(K,I))*PARTSI(JW)
      DSISS(K,I) =  DSID(K,I)+DSISD(K,I)+DSISR(K,I)+DSIS(K,I)-DSIAG(K,I)-DSIEG(K,I)    !+DSIBOD
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                P A R T I C U L A T E   S I L I C A                                            **
!***********************************************************************************************************************************

ENTRY PARTICULATE_SILICA
  PSIAM(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
        PSIAM(K,I) = PSIAM(K,I)+AMR(K,I,JA)*ALG(K,I,JA)*ASI(JA)     !   PSI(K,I)   HA-Z  12/2016
      ENDIF
      END DO
      PSID(K,I)  = PSIDK(JW)*PSI(K,I)
      PSINS(K,I) = PSIS(JW)*(PSI(K-1,I)*DO1(K-1,I)-PSI(K,I)*DO1(K,I))*BI(K,I)/BH2(K,I)
      PSISS(K,I) = PSIAM(K,I)-PSID(K,I)+PSINS(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                            I R O N                                                            **
!***********************************************************************************************************************************

ENTRY IRON
  DO I=IU,ID
    DO K=KT,KB(I)
      FENS(K,I) = FES(JW)*(FE(K-1,I)*DO1(K-1,I)-FE(K,I)*DO1(K,I))*BI(K,I)/BH2(K,I)
      FESR(K,I) = FER(JW)*SODD(K,I)*DO2(K,I)
      FESS(K,I) = FESR(K,I)+FENS(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                       L A B I L E   D O M                                                     **
!***********************************************************************************************************************************

ENTRY LABILE_DOM
  LDOMAP(:,IU:ID) = 0.0; LDOMEP(:,IU:ID) = 0.0; LDOMMAC(:,IU:ID)= 0.0  
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LDOMAP(K,I) = LDOMAP(K,I)+(AER(K,I,JA)+(1.0-APOM(JA))*AMR(K,I,JA))*ALG(K,I,JA)
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))LDOMEP(K,I) = LDOMEP(K,I)+(EER(K,I,JE)+(1.0-EPOM(JE))*EMR(K,I,JE))*EPC(K,I,JE)
      END DO

      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            LDOMMAC(K,I)=LDOMMAC(K,I)+(1.0-MPOM(M))*MMR(K,I,M)*MACRM(JJ,K,I,M)
          END DO
        END IF
      END DO
      LDOMMAC(K,I)=LDOMMAC(K,I)/(DLX(I)*BH(K,I))
      LDOMSS(K,I) = LDOMAP(K,I)+LDOMEP(K,I)-LDOMD(K,I)-LRDOMD(K,I)+LDOMMAC(K,I)

    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   D O M                                                  **
!***********************************************************************************************************************************

ENTRY REFRACTORY_DOM
  DO I=IU,ID
    DO K=KT,KB(I)
      RDOMSS(K,I) = LRDOMD(K,I)-RDOMD(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      L A B I L E   P O M                                                      **
!***********************************************************************************************************************************

ENTRY LABILE_POM
  LPOMAP(:,IU:ID) = 0.0; LPOMEP(:,IU:ID) = 0.0;   LPOMMAC(:,IU:ID) = 0.0; LPZOOIN(:,IU:ID)=0.0;LPZOOOUT(:,IU:ID)=0.0   ! cb 5/19/06
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LPOMAP(K,I) = LPOMAP(K,I)+APOM(JA)*(AMR(K,I,JA)*ALG(K,I,JA))
      END DO
      DO JE=1,NEP                                                          ! cb 5/19/06
        IF (EPIPHYTON_CALC(JW,JE))LPOMEP(K,I) = LPOMEP(K,I)+EPOM(JE)*(EMR(K,I,JE)*EPC(K,I,JE))       ! cb 5/19/06
      END DO                                                               ! cb 5/19/06
      LPOMNS(K,I) = POMS(JW)*(LPOM(K-1,I)-LPOM(K,I))*BI(K,I)/BH2(K,I)

      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          JT=K
          JE=KB(I)
          DO JJ=JT,JE
            LPOMMAC(K,I)=LPOMMAC(K,I)+MPOM(M)*LRPMAC(M)*MMR(K,I,M)*MACRM(JJ,K,I,M)
          END DO
        END IF
      END DO
      LPOMMAC(K,I)=LPOMMAC(K,I)/(DLX(I)*BH(K,I))
      IF(ZOOPLANKTON_CALC)THEN
      DO JZ = 1,NZP
        IF(TGRAZE(K,I,JZ) > 0.0)THEN
          LPZOOOUT(K,I)=LPZOOOUT(K,I)+ZOO(K,I,JZ)*(ZMT(K,I,JZ)+(ZMU(K,I,JZ)-(ZMU(K,I,JZ)*ZEFF(JZ))))
          LPZOOIN(K,I)=LPZOOIN(K,I)+ZOO(K,I,JZ)*PREFP(JZ)*ZMU(K,I,JZ)*LPOM(K,I)/TGRAZE(K,I,JZ)
        ELSE
          LPZOOOUT(K,I)=LPZOOOUT(K,I)+ZOO(K,I,JZ)*(ZMT(K,I,JZ)+(ZMU(K,I,JZ)-(ZMU(K,I,JZ)*ZEFF(JZ))))
          LPZOOIN(K,I)=0.0
        END IF
      END DO
      ENDIF
      LPOMSS(K,I) = LPOMAP(K,I)+LPOMEP(K,I)-LPOMD(K,I)+LPOMNS(K,I)-LRPOMD(K,I)+LPOMMAC(K,I)+LPZOOOUT(K,I)-LPZOOIN(K,I)       ! cb 5/19/06
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   P O M                                                  **
!***********************************************************************************************************************************

ENTRY REFRACTORY_POM
  RPOMMAC(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      RPOMNS(K,I) = POMS(JW)*(RPOM(K-1,I)-RPOM(K,I))*BI(K,I)/BH2(K,I)
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          JT=K
          JE=KB(I)
          DO JJ=JT,JE
            RPOMMAC(K,I)=RPOMMAC(K,I)+MPOM(M)*(1.0-LRPMAC(M))*MMR(K,I,M)*MACRM(JJ,K,I,M)
          END DO
        END IF
      END DO
      RPOMMAC(K,I)=RPOMMAC(K,I)/(DLX(I)*BH(K,I))
      RPOMSS(K,I) = LRPOMD(K,I)+RPOMNS(K,I)-RPOMD(K,I)+RPOMMAC(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                         A L G A E                                                             **
!***********************************************************************************************************************************

ENTRY ALGAE (J)
  AGZT(:,IU:ID,J) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      IF(ZOOPLANKTON_CALC)THEN
      DO JZ = 1,NZP
	  AGZT(K,I,J) = AGZT(K,I,J) + AGZ(K,I,J,JZ)                       ! CB 5/26/07
	  END DO
	  ENDIF
      ASS(K,I,J) = ASR(K,I,J)+(AGR(K,I,J)-AER(K,I,J)-AMR(K,I,J)-ARR(K,I,J))*ALG(K,I,J)-AGZT(K,I,J)	
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                            B I O C H E M I C A L   O 2   D E M A N D                                          **
!***********************************************************************************************************************************

ENTRY BIOCHEMICAL_O2_DEMAND(JBOD)
  IF(JBOD == 1)CBODNS(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      CBODSET = CBODS(JBOD)*(CBOD(K-1,I,JBOD)-CBOD(K,I,JBOD))*BI(K,I)/BH2(K,I)
      CBODNS(K,I)=CBODNS(K,I)+CBODSET
      CBODSS(K,I,JBOD) = -CBODD(K,I,JBOD)*CBOD(K,I,JBOD)+CBODSET
    END DO
  END DO
RETURN

! VARIABLE STOCHIOMETRY FOR CBOD SECTION ! CB 6/6/10
!***********************************************************************************************************************************
!**                                            B I O C H E M I C A L   O 2   D E M A N D   P H O S P H O R U S                    **
!***********************************************************************************************************************************

ENTRY BIOCHEMICAL_O2_DEMAND_P(JBOD)
  IF(JBOD == 1)CBODNSP(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      CBODSET = CBODS(JBOD)*(CBODP(K-1,I,JBOD)-CBODP(K,I,JBOD))*BI(K,I)/BH2(K,I)
      CBODNSP(K,I)=CBODNSP(K,I)+CBODSET
      CBODPSS(K,I,JBOD) = -CBODD(K,I,JBOD)*CBODP(K,I,JBOD)+CBODSET
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                            B I O C H E M I C A L   O 2   D E M A N D   N I T R O G E N                        **
!***********************************************************************************************************************************

ENTRY BIOCHEMICAL_O2_DEMAND_N(JBOD)
  IF(JBOD == 1)CBODNSN(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      CBODSET = CBODS(JBOD)*(CBODN(K-1,I,JBOD)-CBODN(K,I,JBOD))*BI(K,I)/BH2(K,I)
      CBODNSN(K,I)=CBODNSN(K,I)+CBODSET
      CBODNSS(K,I,JBOD) = -CBODD(K,I,JBOD)*CBODN(K,I,JBOD)+CBODSET
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                D I S S O L V E D   O X Y G E N                                                **
!***********************************************************************************************************************************

ENTRY DISSOLVED_OXYGEN
  DOAP(:,IU:ID) = 0.0; DOAR(:,IU:ID) = 0.0; DOEP(:,IU:ID) = 0.0; DOER(:,IU:ID) = 0.0; DOBOD(:,IU:ID) = 0.0
  DOMP(:,IU:ID) = 0.0; DOMR(:,IU:ID) = 0.0; DOZR(:,IU:ID)=0.0    
  doh2s(:,IU:ID) = 0.0; doch4(:,IU:ID) = 0.0; dofe2(:,IU:ID) = 0.0; doMn2(:,IU:ID) = 0.0    ! CEMA

  DO I=IU,ID
    DOSS(KT,I) = 0.0
    DO K=KT,KB(I)
      DO JCB=1,NBOD
        IF(BOD_CALC(JCB))DOBOD(K,I) = DOBOD(K,I)+RBOD(JCB)*CBODD(K,I,JCB)*CBOD(K,I,JCB)
      END DO
      DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
        DOAP(K,I) = DOAP(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*O2AG(JA)
        DOAR(K,I) = DOAR(K,I)+ARR(K,I,JA)*ALG(K,I,JA)*O2AR(JA)
      ENDIF
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))THEN
        DOEP(K,I) = DOEP(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*O2EG(JE)
        DOER(K,I) = DOER(K,I)+ERR(K,I,JE)*EPC(K,I,JE)*O2ER(JE)
        ENDIF
      END DO

      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            DOMP(K,I)=DOMP(K,I)+MGR(JJ,K,I,M)*MACRM(JJ,K,I,M)*O2MG(M)
            DOMR(K,I)=DOMR(K,I)+MRR(K,I,M)*MACRM(JJ,K,I,M)*O2MR(M)
          END DO
        END IF
      END DO
      DOMP(K,I)=DOMP(K,I)/(DLX(I)*BH(K,I))
      DOMR(K,I)=DOMR(K,I)/(DLX(I)*BH(K,I))
      DOPOM(K,I) = (LPOMD(K,I)+RPOMD(K,I))*O2OM(JW)
      DODOM(K,I) = (LDOMD(K,I)+RDOMD(K,I))*O2OM(JW)
      DOOM(K,I)  =  DOPOM(K,I)+DODOM(K,I)+DOBOD(K,I)      
      DONIT(K,I) =  NH4D(K,I)*O2NH4(JW)
            IF(SEDCOMP_EXIST)THEN  ! SW 5/26/15
            DOSED(K,I) =  SEDD(K,I)*O2OM(JW) +SEDD1(K,I)*O2OM(JW)+SEDD2(K,I)*O2OM(JW)   !Amaila
            ELSE
            DOSED(K,I) =  SEDD(K,I)*O2OM(JW)
            ENDIF
      DOSOD(K,I) =  SODD(K,I)*DO3(K,I)
! CEMA start
      do jcg=NGCS,ngce
        if(jcg == nch4)DOCH4(k,i)=ch4d(k,i)*o2ch4      
        if(jcg == nh2s)DOh2s(k,i)=h2sd(k,i)*o2h2s
        if(jcg == nfe2)DOfe2(k,i)=fe2d(k,i)*o2fe2
        if(jcg == nMn2)DOMn2(k,i)=Mn2d(k,i)*o2Mn2
      end do
! CEMA end
     IF(ZOOPLANKTON_CALC)THEN
     DO JZ = 1, NZP
      DOZR(K,I)  = DOZR(K,I)+ZRT(K,I,JZ)*ZOO(K,I,JZ)*O2ZR(JZ)
	 END DO
	 ENDIF
    !DOSS(K,I)  =  DOAP(K,I)+DOEP(K,I)-DOAR(K,I)-DOER(K,I)-DOOM(K,I)-DONIT(K,I)-DOSOD(K,I)-DOSED(K,I)  &
    !                +DOMP(K,I)-DOMR(K,I)-DOZR(K,I)
    DOSS(K,I)  =  DOAP(K,I)+DOEP(K,I)-DOAR(K,I)-DOER(K,I)-DOOM(K,I)-DONIT(K,I)-DOSOD(K,I)-DOSED(K,I)+DOMP(K,I)-DOMR(K,I)-DOZR(K,I)+doch4(k,i)+doh2s(k,i)+dofe2(k,i)+doMn2(k,i)  !&     ! CEMA
                    !+DOMP(K,I)-DOMR(K,I)-DOZR(K,I)+doch4(k,i)+doh2s(k,i)   ! doch4, doh2s,dofe2 already negative...
    END DO
    DOSAT = SATO(T1(KT,I),TDS(KT,I),PALT(I),SALT_WATER(JW))
    IF (.NOT. ICE(I)) THEN
      CALL GAS_TRANSFER
      O2EX       =  REAER(I)
      DOAE(KT,I) = (DOSAT-O2(KT,I))*O2EX*BI(KT,I)/BH2(KT,I)
      DOSS(KT,I) =  DOSS(KT,I)+DOAE(KT,I)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                              I N O R G A N I C   C A R B O N                                                  **
!***********************************************************************************************************************************

ENTRY INORGANIC_CARBON
  TICAP(:,IU:ID) = 0.0; TICEP(:,IU:ID) = 0.0; TICBOD(:,IU:ID) = 0.0
  ticmc(:,iu:id) = 0.0; ticzr(:,iu:id)=0.0  !v3.5
  ticch4=0.0  ! CEMA
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JCB=1,NBOD
        IF(BOD_CALC(JCB))TICBOD(K,I) = TICBOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODC(JCB)
      END DO
      DO JA=1,NAL
        IF(ALG_CALC(JA))TICAP(K,I) = TICAP(K,I)+AC(JA)*(ARR(K,I,JA)-AGR(K,I,JA))*ALG(K,I,JA)
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))TICEP(K,I) = TICEP(K,I)+EC(JE)*(ERR(K,I,JE)-EGR(K,I,JE))*EPC(K,I,JE)
      END DO
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            TICMC(K,I)=TICMC(K,I)+(MRR(K,I,M)-MGR(JJ,K,I,M))*MACRM(JJ,K,I,M)*MC(M)
          END DO
        END IF
      END DO
      TICMC(K,I)=TICMC(K,I)/(DLX(I)*BH(K,I))
      IF(ZOOPLANKTON_CALC)THEN
      DO JZ = 1,NZP
        TICZR(K,I)=TICZR(K,I)+ZRT(K,I,JZ)*ZOO(K,I,JZ)*ZC(JZ) !MLM
      END DO
      ! CEMA start
      do jcg=NGCS,ngce
        if(jcg == nch4)ticCH4=ch4d(k,i)
      end do
      ! CEMA end
      ENDIF

            IF(SEDCOMP_EXIST)THEN  ! SW 5/26/15
               !TICSS(K,I) = TICAP(K,I)+TICEP(K,I)+SEDDC(K,I)+ORGC(JW)*(LPOMD(K,I)+RPOMD(K,I)+LDOMD(K,I)+RDOMD(K,I))                          &
               !    +CO2R(JW)*SODD(K,I)*DO3(K,I)+TICBOD(K,I)+TICMC(K,I)+TICZR(K,I)  + ticch4               &
               !    +sedd1(k,i)*orgc(jw) + sedd2(k,i)*orgc(jw)             ! Amaila
               TICSS(K,I) = TICAP(K,I)+TICEP(K,I)+SEDDC(K,I)+ORGC(JW)*(LPOMD(K,I)+RPOMD(K,I)+LDOMD(K,I)+RDOMD(K,I))                          &
                   +CO2R(JW)*SODD(K,I)*DO3(K,I)+TICBOD(K,I)+TICMC(K,I)+TICZR(K,I)  + ticch4               &
                   +sedd1(k,i)*cbiom(jw) + sedd2(k,i)*cbiom(jw)             ! Amaila, cb 6/7/17

            ELSE
                TICSS(K,I) = TICAP(K,I)+TICEP(K,I)+SEDDC(K,I)+ORGC(JW)*(LPOMD(K,I)+RPOMD(K,I)+LDOMD(K,I)+RDOMD(K,I))                          &
                   +CO2R(JW)*SODD(K,I)*DO3(K,I)+TICBOD(K,I)+TICMC(K,I)+TICZR(K,I)  
            ENDIF
      
    END DO
    IF (.NOT. ICE(I)) THEN
      IF (REAER(I) == 0.0) CALL GAS_TRANSFER
      CO2EX       = REAER(I)*0.923
      CO2REAER(KT,I)=CO2EX*(0.286*EXP(-0.0314*(T2(KT,I))*PALT(I))-CO2(KT,I))*BI(KT,I)/BH2(KT,I)
      TICSS(KT,I) = TICSS(KT,I)+CO2REAER(KT,I)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T                                                          **
!***********************************************************************************************************************************

ENTRY SEDIMENT
  SEDAS(:,IU:ID) = 0.0; LPOMEP(:,IU:ID) = 0.0; SEDCB(:,IU:ID) = 0.0
  DO I=IU,ID
    SEDSI=0.0
    DO K=KT,KB(I)
    IF(K == KB(I))THEN
    BIBH2(K,I)=BI(K,I)/BH2(K,I)
    ELSE
    BIBH2(K,I)=BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
    ENDIF
      DO JA=1,NAL
        IF(ALG_CALC(JA))SEDAS(K,I) = SEDAS(K,I)+MAX(AS(JA),0.0)*ALG(K,I,JA)*BIBH2(K,I)                !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      END DO
      SEDEM = 0.0   ! CB 5/19/06
      DO JE=1,NEP
!        LPOMEP(K,I) = LPOMEP(K,I)+EPOM(JE)*(EMR(K,I,JE)*EPC(K,I,JE))
        IF (EPIPHYTON_CALC(JW,JE))SEDEM = SEDEM+EBR(K,I,JE)/H1(K,I)*EPC(K,I,JE)    ! cb 5/19/06
      END DO
      DO JD=1,NBOD
        IF(BOD_CALC(JD))SEDCB(K,I) = SEDCB(K,I)+MAX(CBODS(JD),0.0)*CBOD(K,I,JD)*BIBH2(K,I)/O2OM(JW)           !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      END DO
      SEDOMS(K,I) = pomS(JW)*(LPOM(K,I)+RPOM(K,I))*BIBH2(K,I)                        !cb 10/22/06
      IF(K==KB(I))THEN
      SEDSO       = 0.0
      ELSE
      SEDSO       = SEDS(JW)*SED(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      ENDIF
      SEDNS(K,I)  = SEDSI-SEDSO
      SEDSI       = SEDSO
      if(k < kb(i))then   ! CEMA sediment in kb layer goes to sediment diagenesis model
        SED(K,I)    = MAX(SED(K,I)+(SEDEM+SEDAS(K,I)+SEDCB(K,I)+SEDOMS(K,I)+SEDNS(K,I)-SEDD(K,I)-SEDBR(K,I))*DLT,0.0)   ! cb 11/30/06
      else if(k == kb(i) .and. .not. sediment_diagenesis)then
        SED(K,I)    = MAX(SED(K,I)+(SEDEM+SEDAS(K,I)+SEDCB(K,I)+SEDOMS(K,I)+SEDNS(K,I)-SEDD(K,I)-SEDBR(K,I))*DLT,0.0)
      end if
    END DO
  END DO
RETURN


!***********************************************************************************************************************************
!**                                                      S E D I M E N T   P H O S P H O R U S                                    **
!***********************************************************************************************************************************

ENTRY SEDIMENTP
  SEDASP(:,IU:ID) = 0.0; LPOMEPP(:,IU:ID) = 0.0; SEDCBP(:,IU:ID) = 0.0
  sdinp(:,iu:id)=0.0   ! CEMA
  DO I=IU,ID
    SEDSIP=0.0
    DO K=KT,KB(I)
      sdalgp=0.0
      DO JA=1,NAL
        IF(ALG_CALC(JA))then 
          SEDASP(K,I) = SEDASP(K,I)+MAX(AS(JA),0.0)*AP(JA)*ALG(K,I,JA)*BIBH2(K,I)          !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
          if(k==kb(i))sdalgp=sdalgp+MAX(AS(JA),0.0)*AP(JA)*ALG(K,I,JA)*BIBH2(K,I)   ! CEMA
        end if
      END DO
      sdepp=0.0
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
          LPOMEPP(K,I) = LPOMEPP(K,I)+EPOM(JE)*EP(JE)*(EMR(K,I,JE)*EPC(K,I,JE))
          if(k==kb(i))sdepp=sdepp+EPOM(JE)*EP(JE)*(EMR(K,I,JE)*EPC(K,I,JE))  ! CEMA
        end if
      END DO
      sdbodp=0.0   ! CEMA
      DO JD=1,NBOD
!        IF(BOD_CALC(JD))SEDCBP(K,I)=SEDCBP(K,I)+MAX(CBODS(JD),0.0)*BODP(JD)*CBOD(K,I,JD)*BIBH2(K,I)      !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
        IF(BOD_CALC(JD))then
          SEDCBP(K,I)=SEDCBP(K,I)+MAX(CBODS(JD),0.0)*CBODP(K,I,JD)*BIBH2(K,I)    ! CB 6/6/10
          if(k==kb(i))sdbodp=sdbodp+MAX(CBODS(JD),0.0)*CBODP(K,I,JD)*BIBH2(K,I)   ! CEMA
        end if
      END DO
      SEDOMSP(K,I) = POMS(JW)*(LPOMP(K,I)+RPOMP(K,I))*BIBH2(K,I)                         !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))  !CB 10/22/06
      IF(K == KB(I))THEN
      SEDSOP       = 0.0
      ELSE
      SEDSOP       = SEDS(JW)*SEDP(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      ENDIF
      SEDNSP(K,I)  = SEDSIP-SEDSOP
      SEDSIP       = SEDSOP
! CEMA start
      SEDPINFLUX(K,I)=(LPOMEPP(K,I)+SEDASP(K,I)+SEDOMSP(K,I)+SEDCBP(K,I))*DLT
      if(k < kb(i))then
        !SEDP(K,I)    = MAX(SEDP(K,I)+(LPOMEPP(K,I)+SEDASP(K,I)+SEDOMSP(K,I)+SEDCBP(K,I)+SEDNSP(K,I)-SEDDP(K,I)   &   ! SW 4/8/16
         SEDP(K,I)    = MAX(SEDP(K,I)+SEDPINFLUX(K,I)+(SEDNSP(K,I)-SEDDP(K,I)   &
                     -SEDBRP(K,I))*DLT,0.0)                                                                 !cb 11/30/06
      else if(k == kb(i) .and. .not. sediment_diagenesis)then
!       SEDP(K,I)    = MAX(SEDP(K,I)+(LPOMEPP(K,I)+SEDASP(K,I)+SEDOMSP(K,I)+SEDCBP(K,I)+SEDNSP(K,I)-SEDDP(K,I)   &
        SEDP(K,I)    = MAX(SEDP(K,I)+SEDPINFLUX(K,I)+(SEDNSP(K,I)-SEDDP(K,I)   &
                     -SEDBRP(K,I))*DLT,0.0)                                                                 !cb 11/30/06
      else if(k == kb(i) .and. sediment_diagenesis)then
        sdinp(k,i)  = sdepp+sdalgp+sdbodp+SEDOMSP(K,I)+SEDNSP(K,I)  ! CEMA calculating P flux to sediment diagnesis model   
      end if
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T   N I T R O G E N                                        **
!***********************************************************************************************************************************

ENTRY SEDIMENTN
  SEDASN(:,IU:ID) = 0.0; LPOMEPN(:,IU:ID) = 0.0; SEDCBN(:,IU:ID) = 0.0
  sdinn(:,IU:ID) = 0.0           ! CEMA
  DO I=IU,ID
    SEDSIN=0.0
    DO K=KT,KB(I)
      sdalgn=0.0
      DO JA=1,NAL
        IF(ALG_CALC(JA))then
          SEDASN(K,I) = SEDASN(K,I)+MAX(AS(JA),0.0)*AN(JA)*ALG(K,I,JA)*BIBH2(K,I)            !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
          if(k==kb(i))sdalgn=sdalgn+MAX(AS(JA),0.0)*AN(JA)*ALG(K,I,JA)*BIBH2(K,I)
        end if
      END DO
      sdepn=0.0
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
          LPOMEPN(K,I) = LPOMEPN(K,I)+EPOM(JE)*EN(JE)*(EMR(K,I,JE)*EPC(K,I,JE))
          if(k==kb(i))sdepn=sdepn+EPOM(JE)*EN(JE)*(EMR(K,I,JE)*EPC(K,I,JE))
        end if
      END DO
      sdbodn=0.0
      DO JD=1,NBOD
!        IF(BOD_CALC(JD))SEDCBN(K,I)=SEDCBN(K,I)+MAX(CBODS(JD),0.0)*BODN(JD)*CBOD(K,I,JD)*BIBH2(K,I)        !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
        IF(BOD_CALC(JD))then
          SEDCBN(K,I)=SEDCBN(K,I)+MAX(CBODS(JD),0.0)*CBODN(K,I,JD)*BIBH2(K,I)    ! CB 6/6/10
          if(k==kb(i))sdbodn=sdbodn+MAX(CBODS(JD),0.0)*CBODN(K,I,JD)*BIBH2(K,I)   ! CEMA
        end if
      END DO
      SEDOMSN(K,I) = POMS(JW)*(LPOMN(K,I)+RPOMN(K,I))*BIBH2(K,I)                           !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))  !CB 10/22/06
      
      IF(K == KB(I)) THEN      ! SW 12/16/07
      !SEDNO3(K,I)  = FNO3SED(JW)*NO3(K,I)*NO3S(JW)*NO3TRM(K,I)*(BI(K,I))/BH2(K,I)      ! CEMA - KB layer N goes to sediment diagensis
      SEDSON       = 0.0
	  ELSE
      SEDNO3(K,I)  = FNO3SED(JW)*NO3(K,I)*NO3S(JW)*NO3TRM(K,I)*(BI(K,I)-BI(K+1,I))/BH2(K,I)
      SEDSON       = SEDS(JW)*SEDN(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
	  ENDIF
      SEDNSN(K,I)  = SEDSIN-SEDSON
      SEDSIN       = SEDSON
! CEMA start
      SEDNINFLUX(K,I)=(LPOMEPN(K,I)+SEDASN(K,I)+SEDOMSN(K,I)+SEDCBN(K,I))*DLT
      if(k < kb(i))then      
!      SEDN(K,I)    = MAX(SEDN(K,I)+(LPOMEPN(K,I)+SEDASN(K,I)+SEDOMSN(K,I)+SEDCBN(K,I)+SEDNSN(K,I)+SEDNO3(K,I)   &
      SEDN(K,I)    = MAX(SEDN(K,I)+SEDNINFLUX(K,I)+(SEDNSN(K,I)+SEDNO3(K,I)   &          
                     -SEDDN(K,I)-SEDBRN(K,I))*DLT,0.0)  !CB 11/30/06                    
      else if(k == kb(i) .and. .not. sediment_diagenesis)then
!       SEDN(K,I)    = MAX(SEDN(K,I)+(LPOMEPN(K,I)+SEDASN(K,I)+SEDOMSN(K,I)+SEDCBN(K,I)+SEDNSN(K,I)+SEDNO3(K,I)   &
        SEDN(K,I)    = MAX(SEDN(K,I)+SEDNINFLUX(K,I)+(SEDNSN(K,I)+SEDNO3(K,I)   &
                     -SEDDN(K,I)-SEDBRN(K,I))*DLT,0.0)  !CB 11/30/06 
      else if(k == kb(i) .and. sediment_diagenesis)then
        sdinn(k,i)  = sdepn+sdalgn+sdbodn++SEDOMSN(K,I)+SEDNSN(K,I)
      end if
! CEMA end
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T   C A R B O N                                            **
!***********************************************************************************************************************************

ENTRY SEDIMENTC
  SEDASC(:,IU:ID) = 0.0; LPOMEPC(:,IU:ID) = 0.0; SEDCBC(:,IU:ID) = 0.0
  sdinc(:,iu:id)=0.0   ! CEMA
  DO I=IU,ID
    DO K=KT,KB(I)
      SEDSIP=0.0
      sdalgc=0.0
      DO JA=1,NAL
        IF(ALG_CALC(JA))then
          SEDASC(K,I) = SEDASC(K,I)+MAX(AS(JA),0.0)*AC(JA)*ALG(K,I,JA)*BIBH2(K,I)             !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
          if(k==kb(i))sdalgc=sdalgc+MAX(AS(JA),0.0)*AC(JA)*ALG(K,I,JA)*BIBH2(K,I)   ! CEMA
        end if
      END DO
      sdepc=0.0
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
          LPOMEPC(K,I) = LPOMEPC(K,I)+EPOM(JE)*EC(JE)*(EMR(K,I,JE)*EPC(K,I,JE))
          if(k==kb(i))sdepc=sdepc+EPOM(JE)*EC(JE)*(EMR(K,I,JE)*EPC(K,I,JE))  ! CEMA
        end if
      END DO
      sdbodc=0.0    ! CEMA
      DO JD=1,NBOD
        IF(BOD_CALC(JD))then
          SEDCBC(K,I)=SEDCBC(K,I)+MAX(CBODS(JD),0.0)*BODC(JD)*CBOD(K,I,JD)*BIBH2(K,I)         !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
          if(k==kb(i))sdbodc=sdbodc+MAX(CBODS(JD),0.0)*BODC(JD)*CBOD(K,I,JD)*BIBH2(K,I)   ! CEMA
        end if
      END DO
      SEDOMSC(K,I) = POMS(JW)*ORGC(JW)*(LPOM(K,I)+RPOM(K,I))*BIBH2(K,I)                     !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))   !CB 10/22/06
      IF(K == KB(I))THEN
      SEDSOC       = 0.0
      ELSE
      SEDSOC       = SEDS(JW)*SEDC(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      ENDIF
      SEDNSC(K,I)  = SEDSIC-SEDSOC
      SEDSIC       = SEDSOC
! CEMA start
      if(k < kb(i))then
      SEDC(K,I)    = MAX(SEDC(K,I)+(LPOMEPC(K,I)+SEDASC(K,I)+SEDOMSC(K,I)+SEDCBC(K,I)+SEDNSC(K,I)-SEDDC(K,I)    &
                     -SEDBRC(K,I))*DLT,0.0)           
      else if(k == kb(i) .and. .not. sediment_diagenesis)then
        SEDC(K,I)    = MAX(SEDC(K,I)+(LPOMEPC(K,I)+SEDASC(K,I)+SEDOMSC(K,I)+SEDCBC(K,I)+SEDNSC(K,I)-SEDDC(K,I)    &
                     -SEDBRC(K,I))*DLT,0.0)
      else if(k == kb(i) .and. sediment_diagenesis)then
        sdinc(k,i)  = sdepc+sdalgc+sdbodc+SEDOMSC(K,I)+SEDNSC(K,I)  ! CEMA calculating C flux to sediment diagnesis model      
      end if
! CEMA end
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T   D E C A Y    R A T E                                   **
!***********************************************************************************************************************************

ENTRY SEDIMENT_DECAY_RATE
  DO I=IU,ID
    SEDSIDK=0.0
    DO K=KT,KB(I)
      SEDSUM=0.0
      SEDSUMK=0.0
      
      DO JA=1,NAL
        IF(ALG_CALC(JA))THEN
        XDUM=MAX(AS(JA),0.0)*ALG(K,I,JA)*BIBH2(K,I)
        SEDSUMK = SEDSUMK + XDUM * LPOMDK(JW)    
        SEDSUM  = SEDSUM  + XDUM
        ENDIF
      END DO
      
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))THEN
        XDUM=EPOM(JE)*(EMR(K,I,JE)*EPC(K,I,JE))
        SEDSUMK = SEDSUMK + XDUM * LPOMDK(JW)
        SEDSUM  = SEDSUM  + XDUM
        ENDIF
      END DO
      
      DO JD=1,NBOD
        IF(BOD_CALC(JD))THEN
        XDUM=MAX(CBODS(JD),0.0)*CBOD(K,I,JD)*BIBH2(K,I)*RBOD(JD)/O2OM(JW)
        SEDSUMK = SEDSUMK+XDUM*CBODD(K,I,JD)               
        SEDSUM  = SEDSUM + XDUM
        ENDIF
      END DO
      
      SEDSUMK = SEDSUMK + POMS(JW)*(LPOM(K,I)*LPOMDK(JW)+RPOM(K,I)*RPOMDK(JW))*BIBH2(K,I)        !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))  ! CB 10/22/06
      SEDSUM  = SEDSUM  + POMS(JW)*(LPOM(K,I)+RPOM(K,I))*BIBH2(K,I)
      
      SEDSUMK = SEDSUMK*DLT
      SEDSUM  = SEDSUM*DLT  
    
      IF((SEDSUM+SED(K,I)) > 0.0)THEN
      SDKV(K,I)    = (SEDSUMK+SED(K,I) * SDKV(K,I))/(SEDSUM+ SED(K,I))
      ELSE
      SDKV(K,I)=0.0
      ENDIF
            
    END DO
  END DO
RETURN

! Amaila start
! additional sediment compartments simulate slow and fast decaying OM left in standing trees
!***********************************************************************************************************************************
!**                                                      S E D I M E N T  1                                                       **
!***********************************************************************************************************************************

ENTRY SEDIMENT1  
  DO I=IU,ID    
    DO K=KT,KB(I)    
      SED1(K,I)    = MAX(SED1(K,I)+(-SEDD1(K,I))*DLT,0.0)      
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T  2                                                       **
!***********************************************************************************************************************************

ENTRY SEDIMENT2  
  DO I=IU,ID    
    DO K=KT,KB(I)    
      SED2(K,I)    = MAX(SED2(K,I)+(-SEDD2(K,I))*DLT,0.0)      
    END DO
  END DO
RETURN

! Amaila end

!***********************************************************************************************************************************
!*                                                         E P I P H Y T O N                                                      **
!***********************************************************************************************************************************

ENTRY EPIPHYTON (J)
  DO I=IU,ID

!** Limiting factor

    LIGHT = (1.0-BETA(JW))*SRON(JW)*SHADE(I)/ESAT(J)
    LAM2  =  LIGHT
    LAM1  =  LIGHT
    DO K=KT,KB(I)

!**** Limiting factor

      LAM1          = LAM2
      LAM2          = LAM1*EXP(-GAMMA(K,I)*H1(K,I))
      FDPO4         = 1.0-FPSS(K,I)-FPFE(K,I)
      ELLIM(K,I,J)  = 2.718282*(EXP(-LAM2)-EXP(-LAM1))/(GAMMA(K,I)*H1(K,I))
      IF (EHSP(J)  /= 0.0) EPLIM(K,I,J) =  FDPO4*PO4(K,I)/(FDPO4*PO4(K,I)+EHSP(J)+NONZERO)
      IF (EHSN(J)  /= 0.0) ENLIM(K,I,J) = (NH4(K,I)+NO3(K,I))/(NH4(K,I)+NO3(K,I)+EHSN(J)+NONZERO)
      IF (EHSSI(J) /= 0.0) ESLIM(K,I,J) =  DSI(K,I)/(DSI(K,I)+EHSSI(J)+NONZERO)
      LIMIT         =  MIN(EPLIM(K,I,J),ENLIM(K,I,J),ESLIM(K,I,J),ELLIM(K,I,J))
      BLIM          =  1.0-(EPD(K,I,J)/(EPD(K,I,J)+EHS(J)))

!**** Sources/sinks

      EGR(K,I,J) =  MIN(ETRM(K,I,J)*EG(J)*LIMIT*BLIM,PO4(K,I)/(EP(J)*DLT*EPD(K,I,J)/H1(KT,I)+NONZERO),(NH4(K,I)+NO3(K,I))/(EN(J)   &
                    *DLT*EPD(K,I,J)/H1(K,I)+NONZERO))
      ERR(K,I,J) =  ETRM(K,I,J)*ER(J)*DO3(K,I)
      EMR(K,I,J) = (ETRMR(K,I,J)+1.0-ETRMF(K,I,J))*EM(J)
      EER(K,I,J) =  MIN((1.0-ELLIM(K,I,J))*EE(J)*ETRM(K,I,J),EGR(K,I,J))
!      EPD(K,I,J) =  MAX(EPD(K,I,J)+EPD(K,I,J)*(EGR(K,I,J)-ERR(K,I,J)-EMR(K,I,J)-EER(K,I,J)-EBR(K,I,J)/(H1(K,I)*0.0025))*DLT,0.0)
      EPD(K,I,J) =  MAX(EPD(K,I,J)+EPD(K,I,J)*(EGR(K,I,J)-ERR(K,I,J)-EMR(K,I,J)-EER(K,I,J)-EBR(K,I,J)/H1(K,I))*DLT,0.00)   ! cb 5/18/06
      if(k == kb(i)) then      ! SW 12/16/07
      EPM(K,I,J) =  EPD(K,I,J)*(BI(K,I)+2.0*H1(K,I))*DLX(I)
	  else
      EPM(K,I,J) =  EPD(K,I,J)*(BI(K,I)-BI(K+1,I)+2.0*H1(K,I))*DLX(I)
	  endif
      EPC(K,I,J) =  EPM(K,I,J)/VOL(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                       L A B I L E   D O M   P H O S P H O R U S                               **
!***********************************************************************************************************************************

ENTRY LABILE_DOM_P
  LDOMPAP(:,IU:ID) = 0.0; LDOMPEP(:,IU:ID) = 0.0; LDOMPMP(:,IU:ID)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LDOMPAP(K,I) = LDOMPAP(K,I)+(AER(K,I,JA)+(1.0-APOM(JA))*AMR(K,I,JA))*ALG(K,I,JA)*AP(JA)
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))LDOMPEP(K,I) = LDOMPEP(K,I)+(EER(K,I,JE)+(1.0-EPOM(JE))*EMR(K,I,JE))*EPC(K,I,JE)*EP(JE)
      END DO
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            LDOMPMP(K,I)=LDOMPMP(K,I)+(1.0-MPOM(M))*MMR(K,I,M)*MACRM(JJ,K,I,M)*MP(M)
          END DO
        END IF
      END DO
      LDOMPMP(K,I)=LDOMPMP(K,I)/(DLX(I)*BH(K,I))
      LDOMPSS(K,I) = LDOMPAP(K,I)+LDOMPEP(K,I)+LDOMPMP(K,I)-(LDOMD(K,I)+LRDOMD(K,I))*ORGPLD(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   D O M   P H O S P H O R U S                            **
!***********************************************************************************************************************************

ENTRY REFRACTORY_DOM_P
  DO I=IU,ID
    DO K=KT,KB(I)
      RDOMPSS(K,I) = LRDOMD(K,I)*ORGPLD(K,I)-RDOMD(K,I)*ORGPRD(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      L A B I L E   P O M   P H O S P H O R U S                                **
!***********************************************************************************************************************************

ENTRY LABILE_POM_P
  LPOMPAP(:,IU:ID) = 0.0;LPOMPMP(:,IU:ID)=0.0;LPZOOINP(:,IU:ID)=0.0; LPZOOOUTP(:,IU:ID)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LPOMPAP(K,I) = LPOMPAP(K,I)+APOM(JA)*(AMR(K,I,JA)*ALG(K,I,JA))*AP(JA)
      END DO
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          JT=K
          JE=KB(I)
          DO JJ=JT,JE
            LPOMPMP(K,I)=LPOMPMP(K,I)+MPOM(M)*LRPMAC(M)*MMR(K,I,M)*MACRM(JJ,K,I,M)*MP(M)
          END DO
        END IF
      END DO
      LPOMPMP(K,I)=LPOMPMP(K,I)/(DLX(I)*BH(K,I))
	IF(ZOOPLANKTON_CALC)THEN
	DO JZ = 1,NZP
      IF(TGRAZE(K,I,JZ) > 0.0)THEN
        LPZOOOUTP(K,I)=LPZOOOUTP(K,I) + ZOO(K,I,JZ)*(ZMT(K,I,JZ)+(ZMU(K,I,JZ)-(ZMU(K,I,JZ)*ZEFF(JZ))))*ZP(JZ)
        LPZOOINP(K,I)=LPZOOINP(K,I) + ZOO(K,I,JZ)*ZMU(K,I,JZ)*PREFP(JZ)*LPOM(K,I)/TGRAZE(K,I,JZ)*ZP(JZ)
      ELSE
        LPZOOOUTP(K,I)=LPZOOOUTP(K,I)+ZOO(K,I,JZ)*(ZMT(K,I,JZ)+(ZMU(K,I,JZ)-(ZMU(K,I,JZ)*ZEFF(JZ))))*ZP(JZ)
        LPZOOINP(K,I)=0.0
      END IF
    END DO
    ENDIF
      LPOMPNS(K,I) = POMS(JW)*(LPOM(K-1,I)*ORGPLP(K-1,I)-LPOM(K,I)*ORGPLP(K,I))*BI(K,I)/BH2(K,I)
      LPOMPSS(K,I) = LPOMPAP(K,I)+LPOMPMP(K,I)-LPOMD(K,I)*ORGPLP(K,I)+LPOMPNS(K,I)-LRPOMD(K,I)*ORGPLP(K,I)
	  IF(ZOOPLANKTON_CALC)THEN
	!  DO JZ = 1,NZP                                           ! KV 4/24/12
	   LPOMPSS(K,I) =LPOMPSS(K,I) + LPZOOOUTP(K,I)-LPZOOINP(K,I)
	!  END DO                                                  ! KV 4/24/12
	  ENDIF

	END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   P O M   P H O S P H O R U S                            **
!***********************************************************************************************************************************

ENTRY REFRACTORY_POM_P
  RPOMPMP(:,IU:ID)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          JT=K
          JE=KB(I)
          DO JJ=JT,JE
            RPOMPMP(K,I)=RPOMPMP(K,I)+MPOM(M)*(1.0-LRPMAC(M))*MMR(K,I,M)*MACRM(JJ,K,I,M)*MP(M)
          END DO
        END IF
      END DO
      RPOMPMP(K,I)=RPOMPMP(K,I)/(DLX(I)*BH(K,I))
      RPOMPNS(K,I) = POMS(JW)*(RPOM(K-1,I)*ORGPRP(K-1,I)-RPOM(K,I)*ORGPRP(K,I))*BI(K,I)/BH2(K,I)
      RPOMPSS(K,I) = LRPOMD(K,I)*ORGPLP(K,I)+RPOMPNS(K,I)-RPOMD(K,I)*ORGPRP(K,I)+RPOMPMP(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                       L A B I L E   D O M   N I T R O G E N                                   **
!***********************************************************************************************************************************

ENTRY LABILE_DOM_N
  LDOMNAP(:,IU:ID) = 0.0; LDOMNEP(:,IU:ID) = 0.0; LDOMNMP(:,IU:ID)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LDOMNAP(K,I) = LDOMNAP(K,I)+(AER(K,I,JA)+(1.0-APOM(JA))*AMR(K,I,JA))*ALG(K,I,JA)*AN(JA)
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))LDOMNEP(K,I) = LDOMNEP(K,I)+(EER(K,I,JE)+(1.0-EPOM(JE))*EMR(K,I,JE))*EPC(K,I,JE)*EN(JE)
      END DO
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            LDOMNMP(K,I)=LDOMNMP(K,I)+(1.0-MPOM(M))*MMR(K,I,M)*MACRM(JJ,K,I,M)*MN(M)
          END DO
        END IF
      END DO
      LDOMNMP(K,I)=LDOMNMP(K,I)/(DLX(I)*BH(K,I))
      LDOMNSS(K,I) = LDOMNAP(K,I)+LDOMNEP(K,I)+LDOMNMP(K,I)-(LDOMD(K,I)+LRDOMD(K,I))*ORGNLD(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   D O M   N I T R O G E N                                **
!***********************************************************************************************************************************

ENTRY REFRACTORY_DOM_N
  DO I=IU,ID
    DO K=KT,KB(I)
      RDOMNSS(K,I) = LRDOMD(K,I)*ORGNLD(K,I)-RDOMD(K,I)*ORGNRD(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      L A B I L E   P O M   N I T R O G E N                                    **
!***********************************************************************************************************************************

ENTRY LABILE_POM_N
  LPOMNAP(:,IU:ID) = 0.0;LPOMNMP(:,IU:ID)=0.0;LPZOOINN(:,IU:ID)=0.0; LPZOOOUTN(:,IU:ID)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LPOMNAP(K,I) = LPOMNAP(K,I)+APOM(JA)*(AMR(K,I,JA)*ALG(K,I,JA))*AN(JA)
      END DO
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          JT=K
          JE=KB(I)
          DO JJ=JT,JE
            LPOMNMP(K,I)=LPOMNMP(K,I)+MPOM(M)*LRPMAC(M)*MMR(K,I,M)*MACRM(JJ,K,I,M)*MN(M)
          END DO
        END IF
      END DO
      LPOMNMP(K,I)=LPOMNMP(K,I)/(DLX(I)*BH(K,I))
	IF(ZOOPLANKTON_CALC)THEN
	DO JZ = 1,NZP
      IF(TGRAZE(K,I,JZ) > 0.0)THEN
        LPZOOOUTN(K,I)=LPZOOOUTN(K,I)+ZOO(K,I,JZ)*(ZMT(K,I,JZ)+(ZMU(K,I,JZ)-(ZMU(K,I,JZ)*ZEFF(JZ))))*ZN(JZ)
        LPZOOINN(K,I)=LPZOOINN(K,I)+ZOO(K,I,JZ)*PREFP(JZ)*ZMU(K,I,JZ)*LPOM(K,I)/TGRAZE(K,I,JZ)*ZN(JZ)
      ELSE
        LPZOOOUTN(K,I)=LPZOOOUTN(K,I)+ZOO(K,I,JZ)*(ZMT(K,I,JZ)+(ZMU(K,I,JZ)-(ZMU(K,I,JZ)*ZEFF(JZ))))*ZN(JZ)
        LPZOOINN(K,I)=0.0
      END IF
	END DO
	ENDIF
      LPOMNNS(K,I) = POMS(JW)*(LPOM(K-1,I)*ORGNLP(K-1,I)-LPOM(K,I)*ORGNLP(K,I))*BI(K,I)/BH2(K,I)
      LPOMNSS(K,I) = LPOMNAP(K,I)+LPOMNMP(K,I)-LPOMD(K,I)*ORGNLP(K,I)+LPOMNNS(K,I)-LRPOMD(K,I)*ORGNLP(K,I) &
            + LPZOOOUTN(K,I)-LPZOOINN(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   P O M   N I T R O G E N                                **
!***********************************************************************************************************************************

ENTRY REFRACTORY_POM_N
  RPOMNMP(:,IU:ID)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          JT=K
          JE=KB(I)
          DO JJ=JT,JE
            RPOMNMP(K,I)=RPOMNMP(K,I)+MPOM(M)*(1.0-LRPMAC(M))*MMR(K,I,M)*MACRM(JJ,K,I,M)*MN(M)
          END DO
        END IF
      END DO
      RPOMNMP(K,I)=RPOMNMP(K,I)/(DLX(I)*BH(K,I))
      RPOMNNS(K,I) = POMS(JW)*(RPOM(K-1,I)*ORGNRP(K-1,I)-RPOM(K,I)*ORGNRP(K,I))*BI(K,I)/BH2(K,I)
      RPOMNSS(K,I) = LRPOMD(K,I)*ORGNLP(K,I)+RPOMNNS(K,I)-RPOMD(K,I)*ORGNRP(K,I)+RPOMNMP(K,I)
    END DO
  END DO
RETURN


!************************************************************************
!**                          M A C R O P H Y T E                       **
!************************************************************************

ENTRY MACROPHYTE(LLM)
  M=LLM
  DO I=IU,ID
    IF(KTICOL(I))THEN
      JT=KTI(I)
    ELSE
      JT=KTI(I)+1
    END IF
    JE=KB(I)
    DO JJ=JT,JE
      IF(JJ.LT.KT)THEN
        COLB=EL(JJ+1,I)
      ELSE
        COLB=EL(KT+1,I)
      END IF
      !COLDEP=ELWS(I)-COLB
       coldep=EL(KT,i)-Z(i)*COSA(JB)-colb  ! cb 3/7/16
      IF(MACRC(JJ,KT,I,M).GT.MMAX(M))THEN
        MGR(JJ,KT,I,M)=0.0
      END IF
      MACSS(JJ,KT,I,M) = (MGR(JJ,KT,I,M)-MMR(KT,I,M)-MRR(KT,I,M))*MACRC(JJ,KT,I,M)
      MACRM(JJ,KT,I,M)   = MACRM(JJ,KT,I,M)+MACSS(JJ,KT,I,M)*DLT*COLDEP*CW(JJ,I)*DLX(I)
    END DO

    DO K=KT+1,KB(I)
      JT=K
      JE=KB(I)
      DO JJ=JT,JE
        IF(MACRC(JJ,K,I,M).GT.MMAX(M))THEN
          MGR(JJ,K,I,M)=0.0
        END IF
        MACSS(JJ,K,I,M) = (MGR(JJ,K,I,M)-MMR(K,I,M)-MRR(K,I,M))*MACRC(JJ,K,I,M)
        IF(MACT(JJ,K,I).GT.MBMP(M).AND.MACT(JJ,K-1,I).LT.MBMP(M).AND.MACSS(JJ,K,I,M).GT.0.0)THEN
          IF(K-1.EQ.KT)THEN
            BMASS=MACSS(JJ,K,I,M)*DLT*H2(K,I)*CW(JJ,I)*DLX(I)
            MACRM(JJ,K-1,I,M)=MACRM(JJ,K-1,I,M)+BMASS
            COLB=EL(KT+1,I)
            !COLDEP=ELWS(I)-COLB
             coldep=EL(KT,i)-Z(i)*COSA(JB)-colb  ! cb 3/7/16
            MACSS(JJ,K-1,I,M)=BMASS/DLT/(COLDEP*CW(JJ,I)*DLX(I)) + MACSS(JJ,K-1,I,M)
          ELSE
            BMASS=MACSS(JJ,K,I,M)*DLT*H2(K,I)*CW(JJ,I)*DLX(I)
            MACRM(JJ,K-1,I,M)=MACRM(JJ,K-1,I,M)+BMASS
            MACSS(JJ,K-1,I,M)=BMASS/DLT/(H2(K-1,I)*CW(JJ,I)*DLX(I))+ MACSS(JJ,K-1,I,M)
          END IF
          MACSS(JJ,K,I,M)=0.0
        ELSE
          BMASSTEST=MACRM(JJ,K,I,M)+MACSS(JJ,K,I,M)*DLT*H2(K,I)*CW(JJ,I)*DLX(I)
          IF(BMASSTEST.GE.0.0)THEN
            MACRM(JJ,K,I,M)   = BMASSTEST
          ELSE
            MACSS(JJ,K,I,M)=-MACRM(JJ,K,I,M)/DLT/(H2(K,I)*CW(JJ,I)*DLX(I))
            MACRM(JJ,K,I,M)=0.0
          END IF
        END IF
      END DO
    END DO
  END DO
  DO I=IU,ID
    TMAC=0.0
    CVOL=0.0
    IF(KTICOL(I))THEN
      JT=KTI(I)
    ELSE
      JT=KTI(I)+1
    END IF
    JE=KB(I)

    DO JJ=JT,JE
      IF(JJ.LT.KT)THEN
        COLB=EL(JJ+1,I)
      ELSE
        COLB=EL(KT+1,I)
      END IF
      !COLDEP=ELWS(I)-COLB
       coldep=EL(KT,i)-Z(i)*COSA(JB)-colb  ! cb 3/7/16
      IF(CW(JJ,I).GT.0.0)THEN
        MACRC(JJ,KT,I,M)=MACRM(JJ,KT,I,M)/(CW(JJ,I)*COLDEP*DLX(I))
      ELSE
        MACRC(JJ,KT,I,M)=0.0
      END IF
      TMAC=TMAC+MACRM(JJ,KT,I,M)
      CVOL=CVOL+CW(JJ,I)*COLDEP*DLX(I)
    END DO

    MAC(KT,I,M)=TMAC/CVOL

    DO K=KT+1,KB(I)
      JT=K
      JE=KB(I)
      TMAC=0.0
      CVOL=0.0
      DO JJ=JT,JE
        IF(CW(JJ,I).GT.0.0)THEN
          MACRC(JJ,K,I,M)=MACRM(JJ,K,I,M)/(CW(JJ,I)*H2(K,I)*DLX(I))
        ELSE
          MACRC(JJ,K,I,M)=0.0
        END IF
        TMAC=TMAC+MACRM(JJ,K,I,M)
        CVOL=CVOL+CW(JJ,I)*H2(K,I)*DLX(I)
      END DO
      MAC(K,I,M)=TMAC/CVOL
    END DO
  END DO

  DO I=IU,ID
    TMAC=0.0
    CVOL=0.0
    DO K=KT,KB(I)
      IF(K.EQ.KT)THEN
        JT=KTI(I)
      ELSE
        JT=K
      END IF
      JE=KB(I)
      DO JJ=JT,JE
        MACT(JJ,K,I)=0.0
        DO MI=1,NMC
          IF(MACROPHYTE_CALC(JW,MI))THEN
            MACT(JJ,K,I)=MACRC(JJ,K,I,MI)+MACT(JJ,K,I)
          END IF
        END DO
      END DO
    END DO
  END DO
  RETURN

!***********************************************************************************************************************************
!*                                                  K I N E T I C   F L U X E S                                                   **
!***********************************************************************************************************************************

ENTRY KINETIC_FLUXES
  DO JAF=1,NAF(JW)
    DO JB=BS(JW),BE(JW)                ! SW 3/9/16
    DO I=CUS(JB),DS(JB)
      DO K=KT,KB(I)
        KFS(K,I,KFCN(JAF,JW)) = KFS(K,I,KFCN(JAF,JW))+KF(K,I,KFCN(JAF,JW))*VOL(K,I)*DLT      ! KF IN G/M3/S x VOL M3 x DT S == G
      END DO
    END DO
    END DO
  ENDDO
  
RETURN

!***********************************************************************************************************************************
!**                                                       p H   C O 2                                                             **
!***********************************************************************************************************************************

ENTRY PH_CO2

! pH and carbonate species

  DO I=IU,ID
    DO K=KT,KB(I)
      CART = TIC(K,I)/12000.0                ! CART=equivalents/liter of C    TIC=mg/l C (MW=12g/mole)
      ALKT = ALK(K,I)/5.0E+04                ! ALK=mg/l as CaCO3 (MW=50 g/mole; EQ=50g/eq))      ALKT=equivalents/l
      T1K  = T1(K,I)+273.15

!**** Ionic strength

      IF (FRESH_WATER(JW)) S2 = 2.5E-05*TDS(K,I)
      IF (SALT_WATER(JW))  S2 = 1.47E-3+1.9885E-2*TDS(K,I)+3.8E-5*TDS(K,I)*TDS(K,I)

!**** Debye-Huckel terms and activity coefficients

      SQRS2  =  SQRT(S2)
      DH1    = -0.5085*SQRS2/(1.0+1.3124*SQRS2)+4.745694E-03+4.160762E-02*S2-9.284843E-03*S2*S2
      DH2    = -2.0340*SQRS2/(1.0+1.4765*SQRS2)+1.205665E-02+9.715745E-02*S2-2.067746E-02*S2*S2
      H2CO3T =  10.0**(0.0755*S2)
      HCO3T  =  10.0**DH1
      CO3T   =  10.0**DH2
      OH     =  HCO3T

!**** Temperature adjustment

      KW = 10.0**(-283.971-0.05069842*T1K+13323.0/T1K+102.24447*LOG10(T1K)-1119669.0/(T1K*T1K))/OH
      K1 = 10.0**(-3404.71/T1K+14.8435-0.032786*T1K)*H2CO3T/HCO3T
      K2 = 10.0**(-2902.39/T1K+ 6.4980-0.023790*T1K)*HCO3T/CO3T

!**** pH evaluation

      PHT = -PH(K,I)-2.1
      IF (PH(K,I) <= 0.0) PHT = -14.0
      INCR = 10.0
      DO N=1,3
        F    = 1.0
        INCR = INCR/10.0
        ITER = 0
        DO WHILE (F > 0.0 .AND. ITER < 12)
          PHT    = PHT+INCR
          HION   = 10.0**PHT
          BICART = CART*K1*HION/(K1*HION+K1*K2+HION*HION)
          F      = BICART*(HION+2.0*K2)/HION+KW/HION-ALKT-HION/OH
          ITER   = ITER+1
        END DO
        PHT = PHT-INCR
      END DO

!**** pH, carbon dioxide, bicarbonate, and carbonate concentrations

      HION      =  10.0**PHT
      PH(K,I)   = -PHT
      CO2(K,I)  =  TIC(K,I)/(1.0+K1/HION+K1*K2/(HION*HION))          ! mg/l as C
      HCO3(K,I) =  TIC(K,I)/(1.0+HION/K1+K2/HION)                    ! mg/l as C
      CO3(K,I)  =  TIC(K,I)/((HION*HION)/(K1*K2)+HION/K2+1.0)        ! mg/l as C
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                       p H   C O 2   N E W                                                     **
!***********************************************************************************************************************************

ENTRY PH_CO2_NEW ! Enhancements added for buffering by ammonia, phosphate, and OM ! SR 01/01/12
! pH and carbonate species
 DO I=IU,ID
   DO K=KT,KB(I)
     T1K = T1(K,I)+273.15
     CART = TIC(K,I)/12011. ! SR 01/01/12
     ALKT = ALK(K,I)/50044. ! SR 01/01/12
     AMMT = NH4(K,I)/14006.74 ! SR 01/01/12
     PHOST = PO4(K,I)/30973.762 ! SR 01/01/12
     OMCT = (LDOM(K,I)+RDOM(K,I))*ORGC(JW)/12011. ! moles carbon per liter from DOM ! SR 01/01/12
     IF (POM_BUFFERING) OMCT = OMCT + (LPOM(K,I)+RPOM(K,I))*ORGC(JW)/12011. ! SR 01/01/12
 !**** Ionic strength
     IF (FRESH_WATER(JW)) S2 = 2.5E-05*TDS(K,I)
     IF (SALT_WATER(JW)) S2 = 1.47E-3+1.9885E-2*TDS(K,I)+3.8E-5*TDS(K,I)*TDS(K,I)
!**** Debye-Huckel terms and activity coefficients
     SQRS2 = SQRT(S2)
     DH1 = -0.5085*SQRS2/(1.0+1.3124*SQRS2)+4.745694E-03+4.160762E-02*S2-9.284843E-03*S2*S2
     DH2 = -2.0340*SQRS2/(1.0+1.4765*SQRS2)+1.205665E-02+9.715745E-02*S2-2.067746E-02*S2*S2
     DH3 = -4.5765*SQRS2/(1.0+1.3124*SQRS2) ! extended Debye-Huckel for PO4 ! SR 01/01/12
     DHH = -0.5085*SQRS2/(1.0+2.9529*SQRS2) ! extended Debye-Huckel for H+ ion ! SR 01/01/12
     H2CO3T = 10.0**(0.0755*S2)
     HCO3T = 10.0**DH1
     CO3T = 10.0**DH2
     PO4T = 10.0**DH3 ! SR 01/01/12
     HT = 10.0**DHH ! activity coefficient for H+ ! SR 01/01/12
     HPO4T = CO3T ! tabled values similar to those for carbonate ! SR 01/01/12
     OHT = HCO3T ! tabled values similar to those for bicarbonate ! SR 01/01/12
     H2PO4T = HCO3T ! tabled values similar to those for bicarbonate ! SR 01/01/12
     NH4T = HCO3T ! tabled values similar to those for bicarbonate ! SR 01/01/12
     NH3T = H2CO3T ! neutral species, set coefficient to same as that for carbonic acid ! SR 01/01/12
     H3PO4T = H2CO3T ! neutral species, set coefficient to same as that for carbonic acid ! SR 01/01/12
!**** Temperature adjustment
     KW = 10.0**(-283.971 -0.05069842*T1K +13323.0/T1K +102.24447*LOG10(T1K) -1119669.0/(T1K*T1K))/OHT
     K1 = 10.0**(-356.3094 -0.06091964*T1K +21834.37/T1K +126.8339 *LOG10(T1K) -1684915 /(T1K*T1K))*H2CO3T/HCO3T
     K2 = 10.0**(-107.8871 -0.03252849*T1K + 5151.79/T1K + 38.92561*LOG10(T1K) - 563713.9/(T1K*T1K))*HCO3T/CO3T
     KAMM = 10.0**(-0.09018 -2729.92/T1K)*NH4T/NH3T ! SR 01/01/12
     KP1 = 10.0**(4.5535 -0.013486*T1K -799.31/T1K)*H3PO4T/H2PO4T ! Bates (1951) ! SR 01/21/12
     KP2 = 10.0**(5.3541 -0.019840*T1K -1979.5/T1K)*H2PO4T/HPO4T ! Bates and Acree (1943) ! SR 01/21/12
     KP3 = 10.0**(-12.38) *HPO4T/PO4T ! Dean (1985) ! SR 01/01/12
!**** pH evaluation
     PHT = -PH(K,I)-2.1
     IF (PH(K,I) <= 0.0) PHT = -14.0
     INCR = 10.0
     DO N=1,3
        F = 1.0
        INCR = INCR/10.0
        ITER = 0
        DO WHILE (F > 0.0 .AND. ITER < 12)
          PHT = PHT+INCR
          HION = 10.0**PHT
          F = CART*K1*(HION+2.0*K2)/(HION*HION+K1*HION+K1*K2)+KW/HION-ALKT-HION/HT ! SR 01/01/12
          IF (AMMONIA_BUFFERING) THEN ! SR 01/01/12
            F = F + AMMT*KAMM/(HION+KAMM) ! SR 01/01/12
          END IF ! SR 01/01/12
          IF (PHOSPHATE_BUFFERING) THEN ! SR 01/01/12
            F = F + PHOST*( KP1*KP2*HION + 2*KP1*KP2*KP3 - HION*HION*HION ) &
                /( HION*HION*HION + KP1*HION*HION + KP1*KP2*HION + KP1*KP2*KP3) ! SR 01/01/12
          END IF ! SR 01/01/12
          IF (OM_BUFFERING) THEN ! SR 01/01/12
            DO JA=1,NAG ! SR 01/01/12
              F = F + OMCT*SDEN(JA)*( 1.0/(1.0+HION*(10.0**PK(JA))) - 1.0/(1.0+(10.0**(PK(JA)-4.5))) ) ! SR 01/01/12
            END DO ! SR 01/01/12
          END IF ! SR 01/01/12
          ITER = ITER+1
        END DO
        PHT = PHT-INCR
     END DO
!**** pH, carbon dioxide, bicarbonate, and carbonate concentrations
     HION = 10.0**PHT
     PH(K,I) = -PHT
     CO2(K,I) = TIC(K,I)/(1.0+K1/HION+K1*K2/(HION*HION))
     HCO3(K,I) = TIC(K,I)/(1.0+HION/K1+K2/HION)
     CO3(K,I) = TIC(K,I)/((HION*HION)/(K1*K2)+HION/K2+1.0)
   END DO
 END DO
RETURN


!**********************************************************
!**           SUBROUTINE ZOOPLANKTON                     **
!**********************************************************

ENTRY ZOOPLANKTON
  DO I=IU,ID
    DO K=KT,KB(I)
	  DO JZ = 1, NZP
            ZGZTOT=0.0                                                                                                   ! KV 5/9/2007
	        DO JJZ = 1,NZP
!             ZGZTOT=ZGZTOT+ZGZ(K,I,JZ,JJZ)*ZOO(K,I,JZ)                                                                   ! KV 5/9/2007
            ZGZTOT=ZGZTOT+ZGZ(K,I,JZ,JJZ)                                                                             ! CB 5/26/07
            END DO
        ZOOSS(K,I,JZ)= (ZMU(K,I,JZ)*ZEFF(JZ)-ZRT(K,I,JZ)-ZMT(K,I,JZ))*ZOO(K,I,JZ) - ZGZTOT   ! OMNIVOROUS ZOOPLANKTON    ! KV 5/9/2007
	  END DO
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                              D E R I V E D   C O N S T I T U E N T S                                          **
!***********************************************************************************************************************************

ENTRY DERIVED_CONSTITUENTS
  APR = 0.0; ATOT = 0.0; TOTSS = 0.0; CHLA = 0.0; CBODU=0.0
  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      DO I=CUS(JB),DS(JB)
        DO K=KT,KB(I)
          DO JA=1,NAL
            IF(ALG_CALC(JA))APR(K,I) = APR(K,I)+(AGR(K,I,JA)-ARR(K,I,JA))*ALG(K,I,JA)*H2(K,I)*DAY
          END DO
        END DO
        DO K=KT,KB(I)
          CBODCT = 0.0; CBODNT = 0.0; CBODPT = 0.0; BODTOT = 0.0; ALGP = 0.0; ALGN = 0.0  ! cb 6/6/10
          DO JA=1,NAL
            IF(ALG_CALC(JA))ATOT(K,I) = ATOT(K,I)+ALG(K,I,JA)
          END DO
          DO IBOD=1,NBOD
          IF(BOD_CALC(IBOD))THEN      
            CBODCt  = CBODCt+CBOD(K,I,IBOD)*BODC(IBOD)    ! cb 6/6/10
            CBODNt  = CBODNt+CBODn(K,I,IBOD)              ! cb 6/6/10
            CBODPt  = CBODPt+CBODp(K,I,IBOD)              ! cb 6/6/10
            BODTOT = BODTOT+CBOD(K,I,IBOD)
            IF(CBODS(IBOD)>0.0)TOTSS(K,I) = TOTSS(K,I)+CBOD(K,I,IBOD)/O2OM(JW)               ! SW 9/5/13  Added particulate CBOD to TSS computation
          ENDIF
          END DO
          DOM(K,I) = LDOM(K,I)+RDOM(K,I)
          POM(K,I) = LPOM(K,I)+RPOM(K,I)
          DOC(K,I) = DOM(K,I)*ORGC(JW)+CBODCt             ! cb 6/6/10
          POC(K,I) = POM(K,I)*ORGC(JW)
          DO JA=1,NAL
          IF(ALG_CALC(JA))THEN
            POC(K,I) = POC(K,I)+ALG(K,I,JA)*AC(JA)
            ALGP     = ALGP+ALG(K,I,JA)*AP(JA)
            ALGN     = ALGN+ALG(K,I,JA)*AN(JA)
          ENDIF
          END DO
          IF(ZOOPLANKTON_CALC)THEN
            DO JZ=1,NZP
                POC(K,I)=POC(K,I)+ZC(JZ)*ZOO(K,I,JZ) !MLM BAULK
                ZOOP=ZOO(K,I,JZ)*ZP(JZ) !MLM BAULK
                ZOON=ZOO(K,I,JZ)*ZN(JZ) !MLM BAULK
                CBODU(K,I) = CBODU(K,I) + O2OM(JW)*ZOO(K,I,JZ)
                TOTSS(K,I) = TOTSS(K,I)+ZOO(K,I,JZ)               ! SW 9/5/13  Added zooplankton to TSS computation
	        END DO
	      ENDIF
          TOC(K,I)   = DOC(K,I)+POC(K,I)
          DOP(K,I)   = LDOM(K,I)*ORGPLD(K,I)+RDOM(K,I)*ORGPRD(K,I)+CBODPT    ! CB 6/6/10
          DON(K,I)   = LDOM(K,I)*ORGNLD(K,I)+RDOM(K,I)*ORGNRD(K,I)+CBODNT    ! CB 6/6/10
          POP(K,I)   = LPOM(K,I)*ORGPLP(K,I)+RPOM(K,I)*ORGPRP(K,I)+ALGP+ZOOP
          PON(K,I)   = LPOM(K,I)*ORGNLP(K,I)+RPOM(K,I)*ORGNRP(K,I)+ALGN+ZOOP
          TOP(K,I)   = DOP(K,I)+POP(K,I)
          TON(K,I)   = DON(K,I)+PON(K,I)
          TKN(K,I)   = TON(K,I)+NH4(K,I)
          CBODU(K,I) = CBODU(K,I)+O2OM(JW)*(DOM(K,I)+POM(K,I)+ATOT(K,I))+BODTOT
          TPSS       = 0.0
          DO JS=1,NSS
            TPSS = TPSS+SS(K,I,JS)*PARTP(JW)
          END DO
          TP(K,I)   =  TOP(K,I)+PO4(K,I)+TPSS
          TN(K,I)   =  TON(K,I)+NH4(K,I)+NO3(K,I)
          O2DG(K,I) = (O2(K,I)/SATO(T1(K,I),TDS(K,I),PALT(I),SALT_WATER(JW)))*100.0          
          DO JA=1,NAL
          IF(ALG_CALC(JA))THEN
            CHLA(K,I)  = CHLA(K,I) +ALG(K,I,JA)/ACHLA(JA)
            TOTSS(K,I) = TOTSS(K,I)+ALG(K,I,JA)
          ENDIF
          END DO
          TOTSS(K,I) = TOTSS(K,I)+TISS(K,I)+POM(K,I)
        END DO
      END DO
    END DO
  END DO
RETURN

!********************************************************************************************************************
!**                                             A L K A L I N I T Y                                                **
!********************************************************************************************************************
ENTRY ALKALINITY ! entire subroutine added ! SR 01/01/12
! According to Stumm and Morgan (1996), table 4.5 on page 173:
! Utilization of ammonium during photosynthesis results in an alkalinity decrease: 14 eq. alk per 16 moles ammonium
! Utilization of nitrate during photosynthesis results in an alkalinity increase: 18 eq. alk per 16 moles nitrate
! Production of ammonium during respiration results in an alkalinity increase: 14 eq. alk per 16 moles ammonium
! Nitrification of ammonium results in an alkalinity decrease: 2 eq. alk per 1 mole ammonium
! Denitrification of nitrate (to nitrogen gas) results in an alkalinity increase: 1 eq. alk per 1 mole nitrate
! Alkalinity is represented as mg/L CaCO3 (MW=100.088). CaCO3 has 2 equivalents of alk per mole.
! Nitrogen has an atomic mass of 14.00674. These numbers account for the factor of 50.044/14.00674 used below.

 
 DO I=IU,ID
   DO K=KT,KB(I)
       if(noncon_alkalinity)then
     ALKSS(K,I) = (50.044/14.00674) * ( 14./16.*(NH4AP(K,I)+NH4EP(K,I)+NH4ZR(K,I)+NH4MR(K,I)-NH4MG(K,I)) &
                  + 18./16.*(NO3AG(K,I)+NO3EG(K,I)) &
                  - 2.*NH4D(K,I) + NO3D(K,I) + NO3SED(K,I)*(1-FNO3SED(JW)) )
       else
         alkss(k,i)=0.0   ! NW 2/11/16
       end if
   END DO
 END DO
 
RETURN

ENTRY DEALLOCATE_KINETICS
  DEALLOCATE (OMTRM,  SODTRM, NH4TRM, NO3TRM, DOM, POM, PO4BOD, NH4BOD, TICBOD, ATRM,   ATRMR,  ATRMF, ETRM,   ETRMR,  ETRMF, BIBH2)
  DEALLOCATE (LAM2M)
RETURN
END SUBROUTINE KINETICS

