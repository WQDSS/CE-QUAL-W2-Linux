SUBROUTINE INIT

USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  USE MACROPHYTEC; USE POROSITYC; USE ZOOPLANKTONC
  USE INITIALVELOCITY; USE BIOENERGETICS; USE CEMAVARS, ONLY:CEMARelatedCode, NXTSEDIAG
  IMPLICIT NONE
  EXTERNAL RESTART_OUTPUT
  INTEGER NB

!***********************************************************************************************************************************
!**                                             Task 1.1: Variable Initialization                                                 **
!***********************************************************************************************************************************

!***********************************************************************************************************************************
!**                                                 Task 1.1.1: Zero Variables                                                    **
!***********************************************************************************************************************************
  IceQSS = 0.0d00; WATER_AGE_ACTIVE=.FALSE.   ! SR 7/27/2017
  KB     = 0;   KBR    = 0;   NAC    = 0;   NTAC   = 0;   NACD   = 0;   NACIN  = 0;   NACTR  = 0;   NACDT  = 0;   NACPR  = 0
  NDSP   = 0;   HMAX   = 0;   KBMAX  = 0;   DLXMAX = 0;   KBQIN  = 0;   KTQIN  = 0;   QGT    =0.0;  QSP    =0.0    ! SW 8/26/15 Initialize Qgt and Qsp for screen output on restart
  NAF    = 0;   TISS   = 0.0;  CSHE   = 0.0; CIN    = 0.0; TIN    = 0.0; EV     = 0.0
  DZ     = 0.0D0; ET     = 0.0; CSHE   = 0.0; A      = 0.0D0; F      = 0.0D0; D      = 0.0D0; C      = 0.0D0; ELTMF  = 0.0
  EL     = 0.0; DX     = 0.0D0; ST     = 0.0D0; SB     = 0.0D0; DZQ    = 0.0D0; TSS    = 0.0
  HSEG   = 0.0; QSS    = 0.0D0; HPG    = 0.0D0; HDG    = 0.0D0; VSH    = 0.0D0; QDH1   = 0.0D0; ADMX   = 0.0D0; DECAY  = 0.0D0
  ADMZ   = 0.0D0; UYBR   = 0.0D0; GRAV   = 0.0D0; FETCH  = 0.0D0; FETCHU = 0.0D0; FETCHD = 0.0D0; DLTTVD = 0.0; ICETHU = 0.0; ICETH1 = 0.0
  ICETH2 = 0.0; P      = 0.0D0; CELRTY = 0.0; TAU1   = 0.0D0; TAU2   = 0.0D0; VOLSR  = 0.0; VOLTR  = 0.0; AF     = 0.0; EF     = 0.0
  ELTMS  = 0.0; DM     = 0.0D0; QIN    = 0.0D0; REAER  = 0.0; ST     = 0.0D0; SB     = 0.0D0; ADMX   = 0.0D0
  ADMZ   = 0.0D0; HPG    = 0.0D0; HDG    = 0.0D0; RHO    = 0.0D0; JDAYTS = 0.0; JDAY1  = 0.0; DEPTHB = 0.0D0; DEPTHM = 0.0D0; UXBR   = 0.0D0
  BHRHO  = 0.0D0; DLMR   = 0.0D0; SRON   = 0.0; CSSB   = 0.0; Q      = 0.0D0; BH1    = 0.0D0; BH2    = 0.0D0; BHR1   = 0.0D0; BHR2   = 0.0D0
  AVHR   = 0.0D0; GRAV   = 0.0D0; KBP    = 0  ; DZT    = 0.0D0; AZT    = 0.0D0; KFJW   = 0.0; QC     = 0.0D0; YSS    = 0.0; YSTS   = 0.0
  QWD    = 0.0D0; QDTR   = 0.0D0; TTR    = 0.0; CTR    = 0.0; TDTR   = 0.0; QOLDS  = 0.0D0; DTPS   = 0.0; VSTS   = 0.0; VSS    = 0.0
  EGT2   = 0.0; HAB    = 100.0; sedpinflux=0.0; sedninflux=0.0 
  RS=0.0;RN=0.0;RB=0.0;RE=0.0;RC=0.0;RANLW=0.0              ! SW 10/20/15
  sdfirstadd=.true.   ! cb 9/3/17
  IF (.NOT. RESTART_IN) THEN
    BR_INACTIVE=.FALSE.
    NSPRF  = 0;   IZMIN  = 0;   KTWB   = 2;   KMIN   = 1;   IMIN   = 1
    T1     = 0.0D0; T2     = 0.0D0; C1     = 0.0D0; C2     = 0.0D0; CD     = 0.0; CIN    = 0.0; C1S    = 0.0; KF     = 0.0; CMBRT  = 0.0
    KFS    = 0.0; U      = 0.0D0; W      = 0.0D0; SU     = 0.0D0; SW     = 0.0D0; SAZ    = 0.0D0; AZ     = 0.0D0; ESBR   = 0.0; EPD    = 0.0
    ETBR   = 0.0; EBRI   = 0.0; DLTLIM = 0.0; VOLEV  = 0.0; VOLPR  = 0.0; VOLDT  = 0.0; VOLWD  = 0.0; CURRENT= 0.0; VOLICE =0.0; ICEBANK=0.0
    VOLUH  = 0.0; VOLDH  = 0.0; VOLIN  = 0.0; VOLOUT = 0.0; VOLSBR = 0.0; VOLTRB = 0.0; TSSS   = 0.0; TSSB   = 0.0; EF     = 0.0
    TSSEV  = 0.0; TSSPR  = 0.0; TSSTR  = 0.0; TSSDT  = 0.0; TSSWD  = 0.0; TSSUH  = 0.0; TSSDH  = 0.0; TSSIN  = 0.0; CSSK   = 0.0
    TSSOUT = 0.0; TSSICE = 0.0; TSSUH1 = 0.0; TSSUH2 = 0.0; CSSUH1 = 0.0; CSSUH2 = 0.0; TSSDH1 = 0.0; TSSDH2 = 0.0
    CSSDH1 = 0.0; CSSDH2 = 0.0; QIND   = 0.0; TIND   = 0.0; CIND   = 0.0; SAVH2  = 0.0; SAVHR  = 0.0; VOLUH2 = 0.0
    AVH1   = 0.0; AVH2   = 0.0; VOLDH2 = 0.0; Z      = 0.0D0; QUH1   = 0.0D0; SED    = 0.0; SEDC   = 0.0; SEDN   = 0.0
    VS     = 0.0; YS     = 0.0; YST    = 0.0; VST    = 0.0; DTP    = 0.0; QOLD   = 0.0; QSUM=0.0; VOLTBR=0.0; DLVOL=0.0; EVBR=0.0 ! SW 7/24/2017
    TPOUT=0.0;TPTRIB=0.0;TPDTRIB=0.0;TPWD=0.0;TPPR=0.0;TPIN=0.0;TNOUT=0.0;TNTRIB=0.0;TNDTRIB=0.0;TNWD=0.0;TNPR=0.0;TNIN=0.0;TN_SEDSOD_NH4=0.0;TP_SEDSOD_PO4=0.0   ! SW 2/19/16  TP_SEDBURIAL=0.0;TN_SEDBURIAL=0.0;
    SEDP   = 0.0; ICETH=0.0; PFLUXIN=0.0;NFLUXIN=0.0                                    ! SW 4/19/10
    ZMIN   = -1000.0
    TKE=0.0                                                     ! SG 10/4/07
    SEDP=0.0;SEDC=0.0;SEDN=0.0
    MACMBRT=0.0; MACRC=0.0;SMACRC=0.0;MAC=0.0;SMAC=0.0;MACRM=0.0;EPM=0.0;macss=0.0 ! cb 3/8/16
    KTICOL=.FALSE.
  END IF
  ANLIM = 1.0; APLIM = 1; ASLIM = 1.0; ALLIM = 1.0; ENLIM = 1.0; EPLIM = 1; ESLIM = 1.0; ELLIM = 1.0; KLOC = 1; ILOC = 1
  MNLIM = 1.0; MPLIM = 1; MCLIM = 1.0; MLLIM = 1.0
  ICESW      = 1.0
  HMIN       = 1.0E10
  DLXMIN     = 1.0E10
  LFPR       = BLANK
  CONV       = BLANK
  CONV1      = BLANK1
  CNAME2     = ADJUSTR(CNAME2)
  CDNAME2    = ADJUSTR(CDNAME2)
  KFNAME2    = ADJUSTR(KFNAME2)
  TITLE(11)  = ' '
  TEXT       = ' '
  ICPL       = 0
  IF (.NOT. CONSTITUENTS) THEN
    if(NBOD > 0)DEALLOCATE(NBODC,NBODN,NBODP)
    NAL = 0; NEP = 0; NSS = 0; NBOD = 0;
  END IF
  DO JW=1,NWB
    GAMMA(:,US(BS(JW)):DS(BE(JW))) = EXH2O(JW)
  END DO
  !***********************************************************************************************************************************
!**                                            Task 1.1.2: Miscellaneous Variables                                                **
!***********************************************************************************************************************************

! Logical controls

  NEW_PAGE              = .TRUE.;  VOLUME_WARNING = .TRUE.;  INITIALIZE_GRAPH = .TRUE.;  UPDATE_GRAPH    = .TRUE.
  ICE                   = .FALSE.; FLUX           = .FALSE.; PUMPON           = .FALSE.          
  TDG_GATE              = .FALSE.; TDG_SPILLWAY   = .FALSE.; INTERNAL_WEIR    = .FALSE.; SURFACE_WARNING = .FALSE.
  WARNING_OPEN          = .FALSE.; PRINT_CONST    = .FALSE.; PRINT_DERIVED    = .FALSE.; ERROR_OPEN      = .FALSE.
  LIMITING_FACTOR       = .FALSE.
  HEAD_BOUNDARY         = .FALSE.; PRINT_HYDRO    = .FALSE.; ONE_LAYER        = .FALSE.; ZERO_SLOPE      = .TRUE.
  INTERNAL_FLOW         = .FALSE.; DAM_INFLOW     = .FALSE.; DAM_OUTFLOW      = .FALSE.; HEAD_FLOW       = .FALSE.     !TC 08/03/04
  UPDATE_RATES          = .FALSE.                                                                                       !TC 08/03/04
  WEIR_CALC             =  NIW > 0; GATES    = NGT > 0; PIPES       = NPI > 0
  PUMPS                 =  NPU > 0; SPILLWAY = NSP > 0; TRIBUTARIES = NTR > 0
  WITHDRAWALS           =  NWD > 0
  VOLUME_BALANCE        = VBC         == '      ON'
  PLACE_QIN             = PQC         == '      ON'; EVAPORATION        = EVC    == '      ON'
  ENERGY_BALANCE        = EBC         == '      ON'; RH_EVAP            = RHEVC  == '      ON'
  PRECIPITATION         = PRC         == '      ON'; RESTART_OUT        = RSOC   == '      ON'
  INTERP_TRIBS          = TRIC        == '      ON'; INTERP_DTRIBS      = DTRIC  == '      ON'
  INTERP_HEAD           = HDIC        == '      ON'; INTERP_INFLOW      = QINIC  == '      ON'
  INTERP_OUTFLOW        = STRIC       == '      ON'; INTERP_WITHDRAWAL  = WDIC   == '      ON'
  INTERP_GATE           = GTIC        == '      ON'       ! cb 8/13/2010
  !INTERP_METEOROLOGY    = METIC       == '      ON'; DOWNSTREAM_OUTFLOW = WDOC   == '      ON'
  INTERP_METEOROLOGY    = METIC       == '      ON'
  if(WDOC== '      ON' .or. WDOC   == '     ONH' .or. WDOC        == '     ONS')DOWNSTREAM_OUTFLOW = .TRUE.  ! cb 4/11/18
  CELERITY_LIMIT        = CELC        == '      ON'; VISCOSITY_LIMIT    = VISC   == '      ON'
  HYDRO_PLOT            = HPLTC       == '      ON'; PRINT_HYDRO        = HPRWBC == '      ON'
  LIMITING_DLT          = HPRWBC(1,:) == '      ON'; FETCH_CALC         = FETCHC == '      ON'
  SCREEN_OUTPUT         = SCRC        == '      ON'; SNAPSHOT           = SNPC   == '      ON'
  CONTOUR               = CPLC        == '      ON'; VECTOR             = VPLC   == '      ON'
  PROFILE               = PRFC        == '      ON'                 !; SPREADSHEET        = SPRC   == '      ON'; 
  
  DO JW=1,NWB
  IF(SPRC(JW)   == '      ON')THEN           ! SW 9/28/2018
      SPREADSHEET(JW)=.TRUE.
  ELSEIF(SPRC(JW)   == '     ONV')THEN
      SPREADSHEET(JW)=.TRUE.
  ENDIF
  ENDDO
  
  TIME_SERIES           = TSRC        == '      ON'; READ_RADIATION     = SROC   == '      ON'
  ICE_CALC              = ICEC        == '      ON' .OR.                  ICEC   == '    ONWB'
  INTERP_EXTINCTION     = EXIC        == '      ON'; READ_EXTINCTION    = EXC    == '      ON'
  NO_INFLOW             = QINC        == '     OFF'; NO_OUTFLOW         = QOUTC  == '     OFF'
  NO_HEAT               = HEATC       == '     OFF'; NO_WIND            = WINDC  == '     OFF'
  SPECIFY_QTR           = TRC         == ' SPECIFY'; DIST_TRIBS         = DTRC   == '      ON'
  IMPLICIT_VISC         = AZSLC       == '     IMP'; UPWIND             = SLTRC  == '  UPWIND'
  ULTIMATE              = SLTRC       == 'ULTIMATE'; TERM_BY_TERM       = SLHTC  == '    TERM'
  MANNINGS_N            = FRICC       == '    MANN'; PLACE_QTR          = TRC    == ' DENSITY'
  LATERAL_SPILLWAY      = LATSPC      /= '    DOWN'; LATERAL_PUMP       = LATPUC /= '    DOWN'
  LATERAL_GATE          = LATGTC      /= '    DOWN'; LATERAL_PIPE       = LATPIC /= '    DOWN'
  TRAPEZOIDAL           = GRIDC       == '    TRAP'                                                                    !SW 07/16/04
  EPIPHYTON_CALC        = CONSTITUENTS .AND. EPIC        == '      ON'
  MASS_BALANCE          = CONSTITUENTS .AND. MBC         == '      ON'
  SUSP_SOLIDS           = CONSTITUENTS .AND. CAC(NSSS)   == '      ON'
  OXYGEN_DEMAND         = CONSTITUENTS .AND. CAC(NDO)    == '      ON'
  SEDIMENT_CALC         = CONSTITUENTS .AND. SEDCc        == '      ON'
  zooplankton_CALC      = CONSTITUENTS .AND. cac(nzooS) == '      ON'
  SEDIMENT_RESUSPENSION = CONSTITUENTS .AND. SEDRC       == '      ON'
  DERIVED_PLOT          = CONSTITUENTS .AND. CDPLTC      == '      ON'
  DERIVED_CALC          = CONSTITUENTS .AND. ANY(CDWBC   == '      ON')
  PH_CALC               = CONSTITUENTS .AND. CDWBC(20,:) == '      ON'
  PRINT_EPIPHYTON       = CONSTITUENTS .AND. EPIPRC      == '      ON' .AND. EPIPHYTON_CALC
  PRINT_SEDIMENT        = CONSTITUENTS .AND. SEDPRC      == '      ON' .AND. SEDIMENT_CALC
  SEDIMENT_CALC1         = CONSTITUENTS .AND. SEDCc1        == '      ON'
  SEDIMENT_CALC2         = CONSTITUENTS .AND. SEDCc2        == '      ON'
  PRINT_SEDIMENT1        = CONSTITUENTS .AND. SEDPRC1      == '      ON' .AND. SEDIMENT_CALC1
  PRINT_SEDIMENT2        = CONSTITUENTS .AND. SEDPRC2      == '      ON' .AND. SEDIMENT_CALC2
  FRESH_WATER           = CONSTITUENTS .AND. WTYPEC      == '   FRESH' .AND. CAC(NTDS) == '      ON'
  SALT_WATER            = CONSTITUENTS .AND. WTYPEC      == '    SALT' .AND. CAC(NTDS) == '      ON'
  CONSTITUENT_PLOT      = CONSTITUENTS .AND. CPLTC       == '      ON' .AND. CAC       == '      ON'
  DETAILED_ICE          = ICE_CALC     .AND. SLICEC      == '  DETAIL'
  LEAP_YEAR             = MOD(YEAR,4) == 0
  ICE_COMPUTATION       = ANY(ICE_CALC)
  END_RUN               = JDAY > TMEND
  UPDATE_KINETICS       = CONSTITUENTS
  ! IF (WEIR_CALC) THEN   ! SW MOVED TO AFTER GEOM DEFINITION BECAUSE OF NEW FEATURE 3/18/16
  !  DO JWR=1,NIW
  !    DO K=2,KMX-1
  !      IF ((K >= KTWR(JWR) .AND. K <= KBWR(JWR))) INTERNAL_WEIR(K,IWR(JWR)) = .TRUE.
  !    END DO
  !  END DO
  !END IF
  WHERE (READ_EXTINCTION)
    EXOM = 0.0
    EXSS = 0.0
  ENDWHERE
  IF (CONSTITUENTS) THEN
    JG_AGE=0
    DO JG=1,NGC   !SR 7/27/2017
      IF(CGQ10(JG) == 0.0 .AND. CG0DK(JG)==-1.0 .AND. CG1DK(JG)==0.0 .AND. CGS(JG)==0.0 .AND. CGLDK(JG)==0.0 .AND. CGKLF(JG)==0.0 .AND. CGCS(JG)==0.0)THEN
        WATER_AGE_ACTIVE=.TRUE.
        JG_AGE=JG
        EXIT
      ENDIF
    ENDDO  
  
    SUSP_SOLIDS   = .FALSE.
    FLUX          =  FLXC   == '      ON'
    PRINT_CONST   =  CPRWBC == '      ON'
    PRINT_DERIVED =  CDWBC  == '      ON'
    IF (ANY(CAC(NSSS:NSSE)  == '      ON')) SUSP_SOLIDS  = .TRUE.
    IF (ANY(CAC(NSSS:NCT)   == '      ON')) UPDATE_RATES = .TRUE.
    DO JA=1,NAL
      LIMITING_FACTOR(JA) = CONSTITUENTS .AND. CAC(NAS-1+JA) == '      ON' .AND. LIMC == '      ON'
      ALG_CALC(JA)              =              CAC(NAS-1+JA) == '      ON'
    END DO
    DO NB=1,NBOD
    BOD_CALC(NB)                =             CAC(NBODS-1+NB)== '      ON'
      BOD_CALCP(NB)                =             CAC(NBODS-1+NBOD+NB)== '      ON'   ! cb 5/19/2011
      BOD_CALCN(NB)                =             CAC(NBODS-1+2*NBOD+NB)== '      ON'   ! cb 5/19/2011
    ENDDO
    DSI_CALC                   =              CAC(NDSI)== '      ON'                                   ! cb 10/12/11
    PO4_CALC                   =              CAC(NPO4)== '      ON'                                   ! cb 10/12/11
    N_CALC                     =              CAC(NNH4)== '      ON'  .OR.  CAC(NNO3)== '      ON'     ! cb 10/12/11
    
  END IF
  JBDAM = 0
  CDHS  = DHS
  DO JB=1,NBR
    UP_FLOW(JB)     = UHS(JB) == 0
    DN_FLOW(JB)     = DHS(JB) == 0
    UP_HEAD(JB)     = UHS(JB) /= 0
    UH_INTERNAL(JB) = UHS(JB) >  0
    IF (UP_HEAD(JB)) THEN
      DO JJB=1,NBR
        IF (ABS(UHS(JB)) >= US(JJB) .AND. ABS(UHS(JB)) <= DS(JJB)) THEN
          IF (ABS(UHS(JB)) == DS(JJB)) THEN
            IF (DHS(JJB) == US(JB)) THEN
              UP_FLOW(JB)       = .TRUE.
              HEAD_FLOW(JB)     = .TRUE.
              INTERNAL_FLOW(JB) = .TRUE.
              UP_HEAD(JB)       = .FALSE.
              UH_INTERNAL(JB)   = .FALSE.
            END IF
            IF (UHS(JB) < 0) THEN
              DO JJJB=1,NBR
                IF (ABS(UHS(JB)) == DS(JJJB)) EXIT                                                                     ! CB 1/2/05
              END DO
              UP_FLOW(JB)       = .TRUE.
              DAM_INFLOW(JB)    = .TRUE.                                                                               !TC 08/03/04
              DAM_OUTFLOW(JJJB) = .TRUE.                                                                               !TC 08/03/04
              INTERNAL_FLOW(JB) = .TRUE.
              UP_HEAD(JB)       = .FALSE.
              UHS(JB)           =  ABS(UHS(JB))
              JBDAM(JJJB)       =  JB
            END IF
          END IF
          EXIT
        END IF
      END DO
    END IF
    DH_INTERNAL(JB) = DHS(JB)  >   0; DN_HEAD(JB)     = DHS(JB)  /=  0; UH_EXTERNAL(JB) = UHS(JB)  == -1
    DH_EXTERNAL(JB) = DHS(JB)  == -1; UQ_EXTERNAL(JB) = UHS(JB)  ==  0; DQ_EXTERNAL(JB) = DHS(JB)  ==  0
    DQ_INTERNAL(JB) = DQB(JB)  >   0; UQ_INTERNAL(JB) = UQB(JB)  >   0 .AND. .NOT. DAM_INFLOW(JB)                      !TC 08/03/04
  END DO
  DO JW=1,NWB
    IF (TKELATPRDCONST(JW) > 0.0) TKELATPRD(JW) = .TRUE.
    IF (STRICK(JW)      > 0.0) STRICKON(JW)  = .TRUE.
    DO JB=BS(JW),BE(JW)
      IF (UH_EXTERNAL(JB) .OR. DH_EXTERNAL(JB)) HEAD_BOUNDARY(JW) = .TRUE.
      IF (SLOPE(JB) /= 0.0)                     ZERO_SLOPE(JW)    = .FALSE.
    END DO
  END DO
  WHERE (CAC == '     OFF') CPLTC = '     OFF'

! Kinetic flux variables

  KFNAME(1)  = 'TISS settling in - source, kg/day            '; KFNAME(2)  = 'TISS settling out - sink, kg/day             '
  KFNAME(3)  = 'PO4 algal respiration - source, kg/day       '; KFNAME(4)  = 'PO4 algal growth - sink, kg/day              '
  KFNAME(5)  = 'PO4 algal net- source/sink, kg/day           '; KFNAME(6)  = 'PO4 epiphyton respiration - source, kg/day   '
  KFNAME(7)  = 'PO4 epiphyton growth - sink, kg/day          '; KFNAME(8)  = 'PO4 epiphyton net- source/sink, kg/day       '
  KFNAME(9)  = 'PO4 POM decay - source, kg/day               '; KFNAME(10) = 'PO4 DOM decay - source, kg/day               '
  KF_PO4_SED=12;KF_PO4_SOD=13
  KFNAME(11) = 'PO4 OM decay - source, kg/day                '; KFNAME(KF_PO4_SED) = 'PO4 sediment decay - source, kg/day          '
  KFNAME(KF_PO4_SOD) = 'PO4 SOD release - source, kg/day             '; KFNAME(14) = 'PO4 net settling  - source/sink, kg/day      '
  KFNAME(15) = 'NH4 nitrification - sink, kg/day             '; KFNAME(16) = 'NH4 algal respiration - source, kg/day       '
  KFNAME(17) = 'NH4 algal growth - sink, kg/day              '; KFNAME(18) = 'NH4 algal net - source/sink, kg/day          '
  KFNAME(19) = 'NH4 epiphyton respiration - source, kg/day   '; KFNAME(20) = 'NH4 epiphyton growth - sink, kg/day          '
  KFNAME(21) = 'NH4 epiphyton net - source/sink, kg/day      '; KFNAME(22) = 'NH4 POM decay - source, kg/day               '
  KFNAME(23) = 'NH4 DOM decay  - source, kg/day              '; KFNAME(24) = 'NH4 OM decay - source, kg/day                '
  KF_NH4_SED=25;KF_NH4_SOD=26
  KFNAME(KF_NH4_SED) = 'NH4 sediment decay - source, kg/day          '; KFNAME(KF_NH4_SOD) = 'NH4 SOD release - source, kg/day             '
  KFNAME(27) = 'NO3 denitrification - sink, kg/day           '; KFNAME(28) = 'NO3 algal growth - sink, kg/day              '
  KFNAME(29) = 'NO3 epiphyton growth - sink, kg/day          '; KFNAME(30) = 'NO3 sediment uptake - sink, kg/day           '
  KFNAME(31) = 'DSi algal growth - sink, kg/day              '; KFNAME(32) = 'DSi epiphyton growth - sink, kg/day          '
  KFNAME(33) = 'DSi PBSi decay - source, kg/day              '; KFNAME(34) = 'DSi sediment decay - source, kg/day          '
  KFNAME(35) = 'DSi SOD release  - source, kg/day            '; KFNAME(36) = 'DSi net settling - source/sink, kg/day       '
  KFNAME(37) = 'PBSi algal mortality  - source, kg/day       '; KFNAME(38) = 'PBSi net settling - source/sink, kg/day      '
  KFNAME(39) = 'PBSi decay - sink, kg/day                    '; KFNAME(40) = 'Fe net settling - source/sink, kg/day        '
  KFNAME(41) = 'Fe sediment release - source, kg/day         '; KFNAME(42) = 'LDOM decay - sink, kg/day                    '
  KFNAME(43) = 'LDOM decay to RDOM - sink, kg/day            '; KFNAME(44) = 'RDOM decay - sink, kg/day                    '
  KFNAME(45) = 'LDOM algal mortality - source, kg/day        '; KFNAME(46) = 'LDOM epiphyton mortality - source, kg/day    '
  KFNAME(47) = 'LPOM decay - sink, kg/day                    '; KFNAME(48) = 'LPOM decay to RPOM - sink, kg/day            '
  KFNAME(49) = 'RPOM decay - sink, kg/day                    '; KFNAME(50) = 'LPOM algal production - source, kg/day       '
  KFNAME(51) = 'LPOM epiphyton production - source, kg/day   '; KFNAME(52) = 'LPOM net settling - source/sink, kg/day      '
  KFNAME(53) = 'RPOM net settling - source/sink, kg/day      '; KFNAME(54) = 'CBOD decay - sink, kg/day                    '
  KFNAME(55) = 'DO algal production  - source, kg/day        '; KFNAME(57) = 'DO algal respiration - sink, kg/day          '   ! cb 6/2/2009
  KFNAME(56) = 'DO epiphyton production  - source, kg/day    '; KFNAME(58) = 'DO epiphyton respiration - sink, kg/day      '   ! cb 6/2/2009
  KFNAME(59) = 'DO POM decay - sink, kg/day                  '; KFNAME(60) = 'DO DOM decay - sink, kg/day                  '
  KFNAME(61) = 'DO OM decay - sink, kg/day                   '; KFNAME(62) = 'DO nitrification - sink, kg/day              '
  KFNAME(63) = 'DO CBOD uptake - sink, kg/day                '; KFNAME(64) = 'DO reaeration - source/sink, kg/day          '
  KF_DO_SED=65; KF_DO_SOD=66
  KFNAME(KF_DO_SED) = 'DO sediment uptake - sink, kg/day            '; KFNAME(KF_DO_SOD) = 'DO SOD uptake - sink, kg/day                 '
  KFNAME(67) = 'TIC algal uptake - sink, kg/day              '; KFNAME(68) = 'TIC epiphyton uptake - sink, kg/day          '
  KFNAME(69) = 'Sediment decay - sink, kg/day                '; KFNAME(70) = 'Sediment algal settling - sink, kg/day       '
  KFNAME(71) = 'Sediment LPOM settling - source,kg/day       '; KFNAME(72) = 'Sediment net settling - source/sink, kg/day  '
  KFNAME(73) = 'SOD decay - sink, kg/day                     '

  KFNAME(74) = 'LDOM P algal mortality - source, kg/day      '; KFNAME(75) = 'LDOM P epiphyton mortality - source, kg/day  '
  KFNAME(76) = 'LPOM P algal production- source, kg/day      '; KFNAME(77) = 'LPOM P net settling - source/sink, kg/day    '
  KFNAME(78) = 'RPOM P net settling - source/sink, kg/day    '
  KFNAME(79) = 'LDOM P algal mortality - source, kg/day      '; KFNAME(80) = 'LDOM P epiphyton mortality - source, kg/day  '
  KFNAME(81) = 'LPOM P algal production- source, kg/day      '; KFNAME(82) = 'LPOM P net settling - source/sink, kg/day    '
  KFNAME(83) = 'RPOM P net settling - source/sink, kg/day    '
  KFNAME(84) = 'Sediment P decay - sink, kg/day              '; KFNAME(85) = 'Sediment algal P settling - source, kg/day   '
  KFNAME(86) = 'Sediment P LPOM settling - source,kg/day     '; KFNAME(87) = 'Sediment net P settling - source/sink, kg/day'
  KFNAME(88) = 'Sediment epiphyton P settling - source,kg/day'
  KFNAME(89) = 'Sediment N decay - sink, kg/day              '; KFNAME(90) = 'Sediment algal N settling - source, kg/day   '
  KFNAME(91) = 'Sediment N LPOM settling - source,kg/day     '; KFNAME(92) = 'Sediment net N settling - source/sink, kg/day'
  KFNAME(93) = 'Sediment epiphyton N settling - source,kg/day'
  KFNAME(94) = 'Sediment C decay - sink, kg/day              '; KFNAME(95) = 'Sediment algal C settling - source, kg/day   '
  KFNAME(96) = 'Sediment C LPOM settling - source,kg/day     '; KFNAME(97) = 'Sediment net C settling - source/sink, kg/day'
  KFNAME(98) = 'Sediment epiphyton C settling - source,kg/day'; 
  KFNAME(99) = 'Sediment N denitrification - source, kg/day  ';
  KFNAME(100) = 'PO4 macrophyte resp - source, kg/day         '
  KFNAME(101) = 'PO4 macrophyte growth - sink, kg/day         '
  KFNAME(102) = 'NH4 macrophyte resp - source, kg/day         '
  KFNAME(103) = 'NH4 macrophyte growth - sink, kg/day         '
  KFNAME(104) = 'LDOM macrophyte mort  - source, kg/day       '
  KFNAME(105) = 'LPOM macrophyte mort  - source, kg/day       '
  KFNAME(106) = 'RPOM macrophyte mort  - source, kg/day       '
  KFNAME(107) = 'DO  macrophyte production  - source, kg/day  '
  KFNAME(108) = 'DO  macrophyte respiration - sink, kg/day    '
  KFNAME(109) = 'TIC macrophyte growth/resp  - S/S, kg/day    '
  KFNAME(110) = 'CBOD settling - sink, kg/day                 '
  KFNAME(111) = 'Sediment CBOD settling - source, kg/day      '
  KFNAME(112) = 'Sediment CBOD P settling - source, kg/day    '
  KFNAME(113) = 'Sediment CBOD N settling - source, kg/day    '
  KFNAME(114) = 'Sediment CBOD C settling - source, kg/day    '
  KFNAME(115) = 'Sediment Burial - sink, kg/day               '
  KF_SED_PBURIAL=116;KF_SED_NBURIAL=117
  KFNAME(KF_SED_PBURIAL) = 'Sediment P Burial - sink, kg/day             '
  KFNAME(KF_SED_NBURIAL) = 'Sediment N Burial - sink, kg/day             '
  KFNAME(118) = 'Sediment C Burial - sink, kg/day             '
  KFNAME(119) = 'CBOD P settling - sink, kg/day               '
  KFNAME(120) = 'CBOD N settling - sink, kg/day               '
  KFNAME(121)=  'CO2 gas exchange air/water interface, kg/day '
  KFNAME(122)=  'DO H2S decay - sink, kg/day                  '
  KFNAME(123)=  'DO CH4 decay - sink, kg/day                  '
  KFNAME(124)=  'H2S gas exchange air/water interface, kg/day '
  KFNAME(125)=  'CH4 gas exchange air/water interface, kg/day '
  KFNAME(126)=  'H2S decay - sink, kg/day                     '
  KFNAME(127)=  'CH4 decay - sink, kg/day                     '
  KFNAME(128)=  'C to Sed. Diagenesis module - source, kg/day '
  KFNAME(129)=  'N to Sed. Diagenesis module - source, kg/day '
  KFNAME(130)=  'P to Sed. Diagenesis module - source, kg/day '
  KFNAME(131)=  'DO sediment diagenesis uptake - sink, kg/day '
  KFNAME(132)=  'Fe(II) oxidation water column - sink, kg/day '
  KFNAME(133)=  'DO Fe(II) oxidation water col.- sink, kg/day '
  KFNAME(134)=  'FeOOH settling from water col. - sink, kg/day'
  KFNAME(135)=  'MnO2 settling from water col. - sink, kg/day '
  KFNAME(136)=  'Mn(II) oxidation water column - sink, kg/day '
  KFNAME(137)=  'DO Mn(II) oxidation water col.- sink, kg/day '
  KFNAME(138)=  'Labile standing biomass decay- sink, kg/day  '
  KFNAME(139)=  'Refract. stand. biomass decay- sink, kg/day  '
  
  

! Convert rates from per-day to per-second

  IF (CONSTITUENTS) THEN   
    AE     = AE    /DAY; AM     = AM    /DAY; AR     = AR    /DAY; AG    = AG    /DAY; AS     = AS    /DAY
    EE     = EE    /DAY; EM     = EM    /DAY; ER     = ER    /DAY; EG    = EG    /DAY; EB     = EB    /DAY
    CGS    = CGS   /DAY; CG0DK  = CG0DK /DAY; CG1DK  = CG1DK /DAY; SSS   = SSS   /DAY; FES    = FES   /DAY
    PSIS   = PSIS  /DAY; POMS   = POMS  /DAY; SDK    = SDK   /DAY; NH4DK = NH4DK /DAY; NO3DK  = NO3DK /DAY
    NO3S   = NO3S  /DAY; PSIDK  = PSIDK /DAY; LRDDK  = LRDDK /DAY; LRPDK = LRPDK /DAY; LDOMDK = LDOMDK/DAY
    LPOMDK = LPOMDK/DAY; RDOMDK = RDOMDK/DAY; RPOMDK = RPOMDK/DAY; KBOD  = KBOD  /DAY; seds   = seds/day   !v3.5
    SDK1    = SDK1   /DAY; SDK2    = SDK2   /DAY  ! Amaila
    SEDB=SEDB/DAY   !CB 11/27/06
    CBODS  = CBODS/DAY   !CB 7/23/07
!    SSFLOC = SSFLOC/DAY                                                                                                 !SR 04/21/13
    DO JW=1,NWB
      SOD(US(BS(JW))-1:DS(BE(JW))+1) = (SOD(US(BS(JW))-1:DS(BE(JW))+1)/DAY)*FSOD(JW)
    END DO
    DO J=1,NEP
      EBR(:,:,J) = EB(J)
    END DO

    MG=MG/DAY
    MR=MR/DAY
    MM=MM/DAY
    ZG=ZG/DAY
    ZR=ZR/DAY
    ZM=ZM/DAY


  END IF

! Convert slope to angle alpha in radians

  ALPHA = ATAN(SLOPE)
  SINA  = SIN(ALPHA)
  SINAC = SIN(ATAN(SLOPEC))
  COSA  = COS(ALPHA)

! Time and printout control variables

  IF (.NOT. RESTART_IN) THEN
    JDAY   = TMSTRT
    ELTM   = TMSTRT*DAY
    DLT    = DLTMAX(1)
    DLTS   = DLT
    MINDLT = DLT
    NIT    = 0
    NV     = 0
    DLTDP  = 1; RSODP = 1; TSRDP = 1; SNPDP = 1; VPLDP = 1; PRFDP = 1
    SPRDP  = 1; CPLDP = 1; SCRDP = 1; FLXDP = 1; WDODP = 1
    NXTSEDIAG=TMSTRT
    DO JW=1,NWB
      DO J=1,NOD
        IF (TMSTRT > SNPD(J,JW)) SNPD(J,JW) = TMSTRT; IF (TMSTRT > PRFD(J,JW)) PRFD(J,JW) = TMSTRT
        IF (TMSTRT > SPRD(J,JW)) SPRD(J,JW) = TMSTRT; IF (TMSTRT > CPLD(J,JW)) CPLD(J,JW) = TMSTRT
        IF (TMSTRT > VPLD(J,JW)) VPLD(J,JW) = TMSTRT; IF (TMSTRT > SCRD(J,JW)) SCRD(J,JW) = TMSTRT
        IF (TMSTRT > FLXD(J,JW)) FLXD(J,JW) = TMSTRT
      END DO
      NXTMSN(JW) = SNPD(SNPDP(JW),JW); NXTMPR(JW) = PRFD(PRFDP(JW),JW); NXTMSP(JW) = SPRD(SPRDP(JW),JW)
      NXTMCP(JW) = CPLD(CPLDP(JW),JW); NXTMVP(JW) = VPLD(VPLDP(JW),JW); NXTMSC(JW) = SCRD(SCRDP(JW),JW)
      NXTMFL(JW) = FLXD(FLXDP(JW),JW)
    END DO
    DO J=1,NOD
      IF (TMSTRT > TSRD(J)) TSRD(J) = TMSTRT; IF (TMSTRT > WDOD(J)) WDOD(J) = TMSTRT
      IF (TMSTRT > RSOD(J)) RSOD(J) = TMSTRT; IF (TMSTRT > DLTD(J)) DLTD(J) = TMSTRT
    END DO
    NXTMTS = TSRD(TSRDP);  NXTMRS = RSOD(RSODP); NXTMWD = WDOD(WDODP)
    
    NXTMWD_SEC = WDOD(WDODP)*86400.    ! cb 4/6/18 frequency test seconds 
    
    if(bioexp)then
        BIODP = 1 ! MLM BIOEXP -BIOENERGETICS
        NXBIO = BIOD(1) !BIOD(BIODP) ! MLM BIOEXP  BIOENERGETICS
        NXTBIO = BIOD(1) ! MLM INITIALIZE THE OUTPUT VARIABLES
        BIOD(NBIO+1:NOD) = TMEND+1.0 ! MLM BIOEXP BIOENERGETICS
    endif
  END IF
  TSRD(NTSR+1:NOD) = TMEND+1.0; WDOD(NWDO+1:NOD) = TMEND+1.0; RSOD(NRSO+1:NOD) = TMEND+1.0; DLTD(NDLT+1:NOD) = TMEND+1.0

  DO JW=1,NWB
    SNPD(NSNP(JW)+1:NOD,JW) = TMEND+1.0; PRFD(NPRF(JW)+1:NOD,JW) = TMEND+1.0; SPRD(NSPR(JW)+1:NOD,JW) = TMEND+1.0
    VPLD(NVPL(JW)+1:NOD,JW) = TMEND+1.0; CPLD(NCPL(JW)+1:NOD,JW) = TMEND+1.0; SCRD(NSCR(JW)+1:NOD,JW) = TMEND+1.0
    FLXD(NFLX(JW)+1:NOD,JW) = TMEND+1.0
  END DO
  JDAYG  = JDAY
  JDAYNX = JDAYG+1
  NXTVD  = JDAY
  DO J=1,NOD      ! SW 12/9/2016
      IF(DLTD(J) == DLTD(J+1))THEN
          DLTDP=DLTDP+1
      ELSE
          EXIT
      ENDIF
  ENDDO
  
  DLTMAXX=DLTMAX(DLTDP)
  DLTFF=DLTF(DLTDP)
  CURMAX = DLTMAX(DLTDP)/DLTF(DLTDP)

! Hydraulic structures

  IF (SPILLWAY) THEN
    DO JS=1,NSP
      IF (LATERAL_SPILLWAY(JS)) THEN
        IF (IDSP(JS) /= 0) THEN
          TRIBUTARIES = .TRUE.
          WITHDRAWALS = .TRUE.
        ELSE
          WITHDRAWALS = .TRUE.
        END IF
      END IF
      DO JB=1,NBR
        IF (IUSP(JS) >= US(JB) .AND. IUSP(JS) <= DS(JB)) EXIT
      END DO
      JBUSP(JS) = JB
      IF (IUSP(JS) == DS(JBUSP(JS)) .AND. .NOT. LATERAL_SPILLWAY(JS)) NST = NST+1
      DO JW=1,NWB
        IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
      END DO
      JWUSP(JS) = JW
      IF (IDSP(JS) > 0) THEN
        DO JB=1,NBR
          IF (IDSP(JS) >= US(JB) .AND. IDSP(JS) <= DS(JB)) EXIT
        END DO
        JBDSP(JS) = JB
        DO JW=1,NWB
          IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
        END DO
        JWDSP(JS) = JW
      ELSE                    
      JBDSP(JS)=1
      JWDSP(JS)=1
      END IF
    END DO
  END IF
  IF (PIPES) THEN
    DO JP=1,NPI
      IF (LATERAL_PIPE(JP)) THEN
        IF (IDPI(JP) /= 0) THEN
          TRIBUTARIES = .TRUE.
          WITHDRAWALS = .TRUE.
        ELSE
          WITHDRAWALS = .TRUE.
        END IF
      END IF
      DO JB=1,NBR
        IF (IUPI(JP) >= US(JB) .AND. IUPI(JP) <= DS(JB)) EXIT
      END DO
      JBUPI(JP) = JB
      IF (IUPI(JP) == DS(JBUPI(JP)) .AND. .NOT. LATERAL_PIPE(JP)) NST = NST+1
      DO JW=1,NWB
        IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
      END DO
      JWUPI(JP) = JW
      IF (IDPI(JP) > 0) THEN
        DO JB=1,NBR
          IF (IDPI(JP) >= US(JB) .AND. IDPI(JP) <= DS(JB)) EXIT
        END DO
        JBDPI(JP) = JB
        DO JW=1,NWB
          IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
        END DO
        JWDPI(JP) = JW
        ELSE  
        JBDPI(JP)=1
        JWDPI(JP)=1
      END IF
    END DO
  END IF
  IF (GATES) THEN
    DO JG=1,NGT
      IF (LATERAL_GATE(JG)) THEN
        IF (IDGT(JG) /= 0) THEN
          TRIBUTARIES = .TRUE.
          WITHDRAWALS = .TRUE.
        ELSE
          WITHDRAWALS = .TRUE.
        END IF
      END IF
      DO JB=1,NBR
        IF (IUGT(JG) >= US(JB) .AND. IUGT(JG) <= DS(JB)) EXIT
      END DO
      JBUGT(JG) = JB
      IF (IUGT(JG) == DS(JBUGT(JG)) .AND. .NOT. LATERAL_GATE(JG)) NST = NST+1
      DO JW=1,NWB
        IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
      END DO
      JWUGT(JG) = JW
      IF (IDGT(JG) > 0) THEN
        DO JB=1,NBR
          IF (IDGT(JG) >= US(JB) .AND. IDGT(JG) <= DS(JB)) EXIT
        END DO
        JBDGT(JG) = JB
        DO JW=1,NWB
          IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
        END DO
        JWDGT(JG) = JW
      ELSE  
        JBDGT(JG)=1             ! SW 3/24/10
        JWDGT(JG)=1             ! SW 3/24/10
      END IF
    END DO
  END IF
  IF (PUMPS) THEN
    DO JP=1,NPU
      IF (LATERAL_PUMP(JP)) THEN
        IF (IDPU(JP) /= 0) THEN
          TRIBUTARIES = .TRUE.
          WITHDRAWALS = .TRUE.
        ELSE
          WITHDRAWALS = .TRUE.
        END IF
      END IF
      DO JB=1,NBR
        IF (IUPU(JP) >= US(JB) .AND. IUPU(JP) <= DS(JB)) EXIT
      END DO
      JBUPU(JP) = JB
      IF (IUPU(JP) ==  DS(JBUPU(JP)) .AND. .NOT. LATERAL_PUMP(JP)) NST = NST+1
      DO JW=1,NWB
        IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
      END DO
      JWUPU(JP) = JW
      IF (IDPU(JP) > 0) THEN
        DO JB=1,NBR
          IF (IDPU(JP) >= US(JB) .AND. IDPU(JP) <= DS(JB)) EXIT
        END DO
        JBDPU(JP) = JB
        DO JW=1,NWB
          IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
        END DO
        JWDPU(JP) = JW
      ELSE  
        JBDPU(JP)=1
        JWDPU(JP)=1
      END IF
    END DO
  END IF

  ALLOCATE (ESTR(NST,NBR),WSTR(NST,NBR),QSTR(NST,NBR),KTSW(NST,NBR),KBSW(NST,NBR),SINKC(NST,NBR),POINT_SINK(NST,NBR),QNEW(KMX),tavg(nst,nbr), tavgw(NWD+NSP+NGT+NPI+NPU),CAVG(NST,NBR,NCT),CDAVG(NST,NBR,NDC),CAVGW(NWD+NSP+NGT+NPI+NPU,NCT),CDAVGW(NWD+NSP+NGT+NPI+NPU,NDC))   
  TAVGW=0.0
  TAVG=0.0
  CAVG=0.0
  CDAVG=0.0
  CAVGW=0.0
  CDAVGW=0.0
  
  QSTR = 0.0
  DO JB=1,NBR
    ESTR(1:NSTR(JB),JB)  = ESTRT(1:NSTR(JB),JB)
    KTSW(1:NSTR(JB),JB)  = KTSWT(1:NSTR(JB),JB)
    KBSW(1:NSTR(JB),JB)  = KBSWT(1:NSTR(JB),JB)
    WSTR(1:NSTR(JB),JB)  = WSTRT(1:NSTR(JB),JB)
    SINKC(1:NSTR(JB),JB) = SINKCT(1:NSTR(JB),JB)
    POINT_SINK(1:NSTR(JB),JB) = SINKC(1:NSTR(JB),JB) == '   POINT'   ! SW 9/27/13
  END DO
  DEALLOCATE (ESTRT,KBSWT,KTSWT,WSTRT,SINKCT)  

! Active constituents, derived constituents, and fluxes

  IF (CONSTITUENTS) THEN
    DO JC=1,NCT
      IF (CAC(JC) == '      ON') THEN
        NAC     = NAC+1
        CN(NAC) = JC
      END IF
      DO JB=1,NBR
        IF (CINBRC(JC,JB) == '      ON') THEN
          NACIN(JB)       = NACIN(JB)+1
          INCN(NACIN(JB),JB) = JC
        END IF
        IF (CDTBRC(JC,JB) == '      ON') THEN
          NACDT(JB)       = NACDT(JB)+1
          DTCN(NACDT(JB),JB) = JC
        END IF
        IF (CPRBRC(JC,JB) == '      ON') THEN
          NACPR(JB)       = NACPR(JB)+1
          PRCN(NACPR(JB),JB) = JC
        END IF
      END DO
      DO JT=1,NTR
        IF (CTRTRC(JC,JT) == '      ON') THEN
          NACTR(JT)          = NACTR(JT)+1
          TRCN(NACTR(JT),JT) = JC
        END IF
      END DO
    END DO
    DO JW=1,NWB
      DO JD=1,NDC
        IF (CDWBC(JD,JW) == '      ON') THEN
          NACD(JW)         = NACD(JW)+1
          CDN(NACD(JW),JW) = JD
        END IF
      END DO
      DO JF=1,NFL
        IF (KFWBC(JF,JW) == '      ON') THEN
          NAF(JW)          = NAF(JW)+1
          KFCN(NAF(JW),JW) = JF
        elseif(ph_calc(jw).and. jf==121)then
          NAF(JW)          = NAF(JW)+1
          KFCN(NAF(JW),JW) = JF       
        elseif(jf >= 121 .and. CEMARelatedCode)then   ! CEMA turning on flux output  ! SW 4/14/2017 ONLY TURN ON IF sediment diagenesis IS ON
          NAF(JW)          = NAF(JW)+1
          KFCN(NAF(JW),JW) = JF       
        END IF
      END DO
    END DO
  END IF

! Starting time

  DEG = CHAR(248)//'C'
  ESC = CHAR(027)
  CALL DATE_AND_TIME (CDATE,CCTIME)
  DO JW=1,NWB
    TITLE(11) = 'Model run at '//CCTIME(1:2)//':'//CCTIME(3:4)//':'//CCTIME(5:6)//' on '//CDATE(5:6)//'/'//CDATE(7:8)//'/'//CDATE(3:4)
    IF (RESTART_IN) TITLE(11) = 'Model restarted at '//CCTIME(1:2)//':'//CCTIME(3:4)//':'//CCTIME(5:6)//' on '//CDATE(5:6)//'/'//     &
                                 CDATE(7:8)//'/'//CDATE(3:4)
  END DO

  Call INITGEOM  ! Call initial geometry

! Density related derived constants

  RHOWCP  = RHOW*CP
  RHOICP  = RHOI*CP
  RHOIRL1 = RHOI*RL1
  DLXRHO = 0.0D0
  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      DLXRHO(US(JB):DS(JB)) = 0.5D0/(DLXR(US(JB):DS(JB))*RHOW)
      IF (UP_HEAD(JB)) DLXRHO(US(JB)-1) = 0.5D0/(DLXR(US(JB))*RHOW)
    END DO
  END DO

! Transport interpolation multipliers

  DO JW=1,NWB
    CALL INTERPOLATION_MULTIPLIERS
  END DO

CALL INITCOND    ! CALL ROUTINE TO SET UP IC

  IF (WEIR_CALC) THEN   ! MOVED FROM ABOVE AFTER GEOMETRY SETUP  SW 3/16/18
    DO JWR=1,NIW
        IF (EKTWR(JWR) == 0.0) THEN  
            DO JW=1,NWB
            IF(IWR(JWR) >= US(BS(JW)) .AND. IWR(JWR) <= DS(BE(JW)))THEN
            KTWR(JWR)=KTWB(JW)
            EXIT
            ENDIF
            ENDDO
        ELSE  
          KTWR(JWR) = INT(EKTWR(JWR))  
        END IF 
        IF (EKBWR(JWR) <= 0.0) THEN  
            DO K=KTWR(JWR),KB(IWR(JWR))  
            IF (DEPTHB(K,IWR(JWR)) >= ABS(EKBWR(JWR))) THEN
                KBWR(JWR)=K
                EXIT  
            ENDIF
            END DO   
        ELSE  
          KBWR(JWR) = INT(EKBWR(JWR))  
        END IF  
        
      DO K=2,KMX-1
        IF ((K >= KTWR(JWR) .AND. K <= KBWR(JWR))) INTERNAL_WEIR(K,IWR(JWR)) = .TRUE.
      END DO
    END DO
  END IF


! Saved variables for autostepping

  IF (.NOT. RESTART_IN) THEN
    SZ    = Z
    SU    = U
    SW    = W
    SAZ   = AZ
    SKTI  = KTI
    SBKT  = BKT
    SAVH2 = AVH2
    SAVHR = AVHR
  END IF
  CALL GREGORIAN_DATE
!  CALL TIME_VARYING_DATA
!  if(.not. once_through)
  CALL TIME_VARYING_DATA                 ! w2-ressim
!  CALL READ_INPUT_DATA (NXTVD)
!  if(.not. once_through)
  CALL READ_INPUT_DATA (NXTVD)           ! w2-ressim
  IF (CONSTITUENTS) THEN
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IU = US(JB)
        ID = DS(JB)
        CALL TEMPERATURE_RATES
        CALL KINETIC_RATES
      END DO
    END DO
  END IF

RETURN
END SUBROUTINE INIT
