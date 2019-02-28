
!***********************************************************************************************************************************
!**                                  S U B R O U T I N E   T I M E   V A R Y I N G   D A T A                                      **
!***********************************************************************************************************************************

SUBROUTINE TIME_VARYING_DATA
  USE GLOBAL;  USE SURFHE; USE SCREENC; USE TVDC; USE LOGICC; USE SELWC; USE STRUCTURES; USE NAMESC
  USE KINETIC, ONLY:EXH2O; USE SHADEC; USE MAIN, ONLY: PUMPS, SEGNUM
  IMPLICIT NONE

! Type declaration

  CHARACTER(1)                           :: INFORMAT
  CHARACTER(2)                           :: INFORMAT2
  REAL                                   :: NXQGT2                  
  REAL                                   :: NXQWD1, NXQWD2, NXQGT,  NXTVD, NXQPT
  REAL                                   :: NXWSC
  REAL(R8)                               :: RATIO,QRATIO,TRATIO,CRATIO,HRATIO
  REAL(R8),ALLOCATABLE, DIMENSION(:)     :: QDTRO,  TDTRO,  ELUHO,  ELDHO,  QWDO,   QTRO,   TTRO,   QINO,   TINO
  REAL,    ALLOCATABLE, DIMENSION(:)     :: TAIRNX, TDEWNX, PHINX,  WINDNX, SRONX,  CLOUDNX,BGTNX,XX
  REAL,    ALLOCATABLE, DIMENSION(:)     :: NXEXT1, NXEXT2, EXTNX,  EXTO
  REAL,    ALLOCATABLE, DIMENSION(:)     :: TAIRO,  TDEWO,  PHIO,   WINDO,  SROO,   CLOUDO
  REAL(R8),ALLOCATABLE, DIMENSION(:)     :: QDTRNX, TDTRNX, PRNX,   TPRNX,  ELUHNX, ELDHNX, QWDNX,  QTRNX,  TTRNX,  QINNX,  TINNX
  REAL,    ALLOCATABLE, DIMENSION(:)     :: NXQTR1, NXTTR1, NXCTR1, NXQIN1, NXTIN1, NXCIN1, NXQDT1, NXTDT1, NXCDT1, NXDYNS
  REAL,    ALLOCATABLE, DIMENSION(:)     :: NXPR1,  NXTPR1, NXCPR1, NXEUH1, NXTUH1, NXCUH1, NXEDH1, NXTDH1, NXCDH1, NXQOT1, NXMET1
  REAL,    ALLOCATABLE, DIMENSION(:)     :: NXQTR2, NXTTR2, NXCTR2, NXQIN2, NXTIN2, NXCIN2, NXQDT2, NXTDT2, NXCDT2
  REAL,    ALLOCATABLE, DIMENSION(:)     :: NXPR2,  NXTPR2, NXCPR2, NXEUH2, NXTUH2, NXCUH2, NXEDH2, NXTDH2, NXCDH2, NXQOT2, NXMET2
  REAL,    ALLOCATABLE, DIMENSION(:)     :: WSCNX, BPNX,EPU2,EONPU2,EOFFPU2,QPU2,NXPUMP  
  REAL,    ALLOCATABLE, DIMENSION(:,:)   :: NXESTRT
  REAL(R8),ALLOCATABLE, DIMENSION(:,:)   :: CTRO,   CINO,   QOUTO,  CDTRO,  TUHO,   TDHO,   QSTRO
  REAL(R8),ALLOCATABLE, DIMENSION(:,:)   :: CTRNX,  CINNX,  QOUTNX, CDTRNX, CPRNX,  TUHNX,  TDHNX,  QSTRNX
  REAL,    ALLOCATABLE, DIMENSION(:,:,:) :: CUHO,   CDHO,   CUHNX,  CDHNX
  INTEGER                                :: WDQ,    GTQ,    WSH,    SHD, PIPED, IOPENPIPE,L,NJS                                    
  INTEGER                                :: NPT,J,JT,JAC,JS,K,JG,JWD
  INTEGER, ALLOCATABLE, DIMENSION(:)     :: TRQ,    TRT,    TRC,    INQ,    DTQ,    PRE,    UHE,    DHE,    INFT,   DTT,PUMPD, JJS
  INTEGER, ALLOCATABLE, DIMENSION(:)     :: PRT,    UHT,    DHT,    INC,    DTC,    PRC,    UHC,    DHC,    OTQ,    MET,    EXT, ODYNS,DYNPUMPF
  LOGICAL, ALLOCATABLE, DIMENSION(:)     :: INFLOW_CONST, TRIB_CONST, DTRIB_CONST, PRECIP_CONST,OTQF, TRCF, DTCF, INCF, PRCF, METF,DYNEF, TRQF, TRTF, DTTF, DTQF, INQF, INTF, PRQF, PRTF,EXTF
  LOGICAL                                :: WDQF,WSHF, GATEF   ! SW 9/26/2017
  INTEGER, ALLOCATABLE, DIMENSION(:)     :: EUHF,TUHF,CUHF,EDHF,TDHF,CDHF    ! =0 Old format for head BCs, =1 Time series format no vertical variation, =2 csv format vertical variation             SW 2/28/17
  SAVE
! Allocation declarations

  IF(NAC<1)THEN
      ALLOCATE(XX(1))
  ELSE
      ALLOCATE(XX(NCT))
  ENDIF
  
  ALLOCATE (NXQTR1(NTR), NXTTR1(NTR), NXCTR1(NTR), NXQIN1(NBR), NXTIN1(NBR), NXCIN1(NBR), NXQDT1(NBR), NXTDT1(NBR), NXCDT1(NBR))
  ALLOCATE (NXPR1(NBR),  NXTPR1(NBR), NXCPR1(NBR), NXEUH1(NBR), NXTUH1(NBR), NXCUH1(NBR), NXEDH1(NBR), NXTDH1(NBR), NXCDH1(NBR))
  ALLOCATE (NXQOT1(NBR), NXMET1(NWB), NXQTR2(NTR), NXTTR2(NTR), NXCTR2(NTR), NXQIN2(NBR), NXTIN2(NBR), NXCIN2(NBR), NXQDT2(NBR))
  ALLOCATE (NXTDT2(NBR), NXCDT2(NBR), NXPR2(NBR),  NXTPR2(NBR), NXCPR2(NBR), NXEUH2(NBR), NXTUH2(NBR), NXCUH2(NBR), NXEDH2(NBR))
  ALLOCATE (NXTDH2(NBR), NXCDH2(NBR), NXQOT2(NBR), NXMET2(NWB),DYNPUMPF(NPU))
  ALLOCATE (WSCNX(IMX), PRCF(NBR), METF(NWB), ODYNS(NBR),NXDYNS(NBR),NXESTRT(NST,NBR), PRTF(NBR), PRQF(NBR),EXTF(NWB))
  ALLOCATE (QDTRO(NBR),  TDTRO(NBR),  ELUHO(NBR),  ELDHO(NBR),  QWDO(NWD),   QTRO(NTR),   TTRO(NTR),   QINO(NBR),   TINO(NBR))
  ALLOCATE (QDTRNX(NBR), TDTRNX(NBR), PRNX(NBR),   TPRNX(NBR),  ELUHNX(NBR), ELDHNX(NBR), QWDNX(NWD),  QTRNX(NTR),  TTRNX(NTR))
  ALLOCATE (QINNX(NBR),  TINNX(NBR),  SROO(NWB),   TAIRO(NWB),  TDEWO(NWB),  CLOUDO(NWB), PHIO(NWB),   WINDO(NWB),  TAIRNX(NWB))
  ALLOCATE (TDEWNX(NWB), CLOUDNX(NWB),PHINX(NWB),  WINDNX(NWB), SRONX(NWB),  BGTNX(NGT), BPNX(NPI))
  ALLOCATE (TRQ(NTR),    TRT(NTR),    TRC(NTR),    INQ(NBR),    DTQ(NBR),    PRE(NBR),    UHE(NBR),    DHE(NBR),    INFT(NBR))
  ALLOCATE (DTT(NBR),    PRT(NBR),    UHT(NBR),    DHT(NBR),    INC(NBR),    DTC(NBR),    PRC(NBR),    UHC(NBR),    DHC(NBR))
  ALLOCATE (OTQ(NBR),    MET(NWB),    EXT(NWB))
  ALLOCATE (NXEXT1(NWB), NXEXT2(NWB), EXTNX(NWB),  EXTO(NWB))
  ALLOCATE (CTRO(NCT,NTR),   CINO(NCT,NBR),  QOUTO(KMX,NBR),  CDTRO(NCT,NBR),  TUHO(KMX,NBR),  TDHO(KMX,NBR),  QSTRO(NST,NBR))
  ALLOCATE (CTRNX(NCT,NTR),  CINNX(NCT,NBR), QOUTNX(KMX,NBR), CDTRNX(NCT,NBR), CPRNX(NCT,NBR), TUHNX(KMX,NBR), TDHNX(KMX,NBR))
  ALLOCATE (QSTRNX(NST,NBR), PUMPD(NPU), NXPUMP(NPU),EPU2(NPU),EONPU2(NPU),EOFFPU2(NPU),QPU2(NPU),DYNEF(NBR))
  ALLOCATE (CUHO(KMX,NCT,NBR), CDHO(KMX,NCT,NBR), CUHNX(KMX,NCT,NBR), CDHNX(KMX,NCT,NBR),JJS(NST))
  ALLOCATE (INFLOW_CONST(NBR), TRIB_CONST(NTR),   DTRIB_CONST(NBR),   PRECIP_CONST(NBR), OTQF(NBR), TRCF(NTR), DTCF(NBR), INCF(NBR),TRQF(NTR), TRTF(NTR), DTTF(NBR), DTQF(NBR), INQF(NBR), INTF(NBR))
  ALLOCATE (EUHF(NBR), TUHF(NBR), CUHF(NBR), EDHF(NBR), TDHF(NBR), CDHF(NBR))

  NXPR1  = 0.0; NXQTR1 = 0.0; NXTTR1 = 0.0; NXCTR1 = 0.0; NXQIN1 = 0.0; NXTIN1 = 0.0; NXCIN1 = 0.0; NXQDT1 = 0.0; NXTDT1 = 0.0
  NXCDT1 = 0.0; NXTPR1 = 0.0; NXCPR1 = 0.0; NXEUH1 = 0.0; NXTUH1 = 0.0; NXCUH1 = 0.0; NXEDH1 = 0.0; NXTDH1 = 0.0; NXCDH1 = 0.0
  NXQOT1 = 0.0; NXMET1 = 0.0; QSTRNX = 0.0; CDTRNX = 0.0; CTRNX  = 0.0; CINNX  = 0.0; CPRNX  = 0.0; CUHNX  = 0.0; CDHNX  = 0.0
  QINNX  = 0.0; TINNX  = 0.0; CINNX  = 0.0; NXWSC  = 0.0

! Set logical variables

  INFLOW_CONST = CONSTITUENTS .AND. NACIN > 0; TRIB_CONST   = CONSTITUENTS .AND. NACTR > 0
  DTRIB_CONST  = CONSTITUENTS .AND. NACDT > 0; PRECIP_CONST = CONSTITUENTS .AND. NACPR > 0
  OTQF = .FALSE.
  WSHF = .FALSE.
  WDQF = .FALSE.
  TRCF = .FALSE.
  TRQF = .FALSE.
  TRTF = .FALSE.
  DTCF = .FALSE.
  INCF = .FALSE.
  PRCF = .FALSE.
  PRTF = .FALSE.
  PRQF = .FALSE.
  METF = .FALSE.
  DYNEF= .FALSE.
  DTTF = .FALSE.
  DTQF = .FALSE.
  INQF = .FALSE.
  INTF = .FALSE.
  EXTF = .FALSE.
  DYNPUMPF=.FALSE. ; GATEF = .FALSE.   ! SW 9/26/2017
  EUHF = 0;   TUHF=0;   CUHF=0;   EDHF=0;  TDHF=0;  CDHF=0               ! SW 2/28/17

! Open input files

  NPT = NUNIT
  SHD = NPT; NPT=NPT+1
  OPEN (SHD,FILE=SHDFN,STATUS='OLD')
  READ( SHD,'(A1)')INFORMAT
  IF(INFORMAT=='$')THEN
      READ(SHD,'(/)')
      DO I=1,IMX
          READ (SHD,*)J,SHADEI(I)    ! SW 3/14/2018 ADDED TO BE COMPATIBLE WITH PREPROCESSOR
          IF(SHADEI(I)<0.0)THEN
          BACKSPACE(SHD)
          READ (SHD,*) J,SHADEI(I),TTLB(I),TTRB(I),CLLB(I),CLRB(I),SRLB1(I),SRLB2(I),SRRB1(I),SRRB2(I),(TOPO(I,J),J=1,IANG),SRFJD1(I),SRFJD2(I)
          ENDIF
      ENDDO
  ELSE
  READ (SHD,'(//(8X,29F8.0))') (SHADEI(I),TTLB(I),TTRB(I),CLLB(I),CLRB(I),SRLB1(I),SRLB2(I),SRRB1(I),SRRB2(I),        &
                                      (TOPO(I,J),J=1,IANG),SRFJD1(I),SRFJD2(I),I=1,IMX)
  ENDIF
  SHADE = SHADEI
  WSH = NPT; NPT = NPT+1
  OPEN (WSH,FILE=WSCFN,STATUS='OLD')
  READ( WSH,'(A1)')INFORMAT
  IF(INFORMAT=='$')WSHF=.TRUE.
  IF(WSHF)THEN
  READ (WSH,'(/)')
  READ (WSH,*) NXWSC,(WSCNX(I),I=1,IMX)
  WSC = WSCNX
  READ (WSH,*)    NXWSC,(WSCNX(I),I=1,IMX)  
  ELSE
  READ (WSH,'(//10F8.0:/(8X,9F8.0))') NXWSC,(WSCNX(I),I=1,IMX)
  WSC = WSCNX
  READ (WSH,'(10F8.0:/(8X,9F8.0))')    NXWSC,(WSCNX(I),I=1,IMX)
  ENDIF
  DO JW=1,NWB
    MET(JW) = NPT; NPT = NPT+1
    OPEN (MET(JW),FILE=METFN(JW),STATUS='OLD')
    READ(MET(JW),'(A1)')INFORMAT
    IF(INFORMAT=='$')METF(JW)=.TRUE.
    IF (READ_RADIATION(JW)) THEN
        IF(METF(JW))THEN
        READ (MET(JW),'(/)')
        READ (MET(JW),*) NXMET2(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW),SRONX(JW)
        ELSE
        READ (MET(JW),'(//10F8.0)') NXMET2(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW),SRONX(JW)
        ENDIF
      SRONX(JW) = SRONX(JW)*REFL
      SRON(JW)  = SRONX(JW)
      SROO(JW)  = SRON(JW)
    ELSE
      IF(METF(JW))THEN
      READ (MET(JW),'(/)')
      READ (MET(JW),*) NXMET2(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW)
      ELSE
      READ (MET(JW),'(//10F8.0)') NXMET2(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW)
      ENDIF
    END IF
    TAIR(JW)   = TAIRNX(JW)
    TDEW(JW)   = TDEWNX(JW)
    WIND(JW)   = WINDNX(JW)
    PHI(JW)    = PHINX(JW)
    CLOUD(JW)  = CLOUDNX(JW)
    TAIRO(JW)  = TAIRNX(JW)
    TDEWO(JW)  = TDEWNX(JW)
    WINDO(JW)  = WINDNX(JW)
    PHIO(JW)   = PHINX(JW)
    CLOUDO(JW) = CLOUDNX(JW)
    IF (PHISET > 0) PHI(JW)  = PHISET
    IF (PHISET > 0) PHIO(JW) = PHISET
    IF (READ_RADIATION(JW)) THEN
            IF(METF(JW))THEN
            READ (MET(JW),*) NXMET1(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW),SRONX(JW)
            ELSE
            READ (MET(JW),'(10F8.0)') NXMET1(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW),SRONX(JW)
            ENDIF
      SRONX(JW) = SRONX(JW)*REFL
        ELSE
            IF(METF(JW))THEN
            READ (MET(JW),*) NXMET1(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW)
            ELSE
            READ (MET(JW),'(10F8.0)') NXMET1(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW)
            ENDIF
    END IF
    IF (READ_EXTINCTION(JW)) THEN
      EXT(JW) = NPT; NPT = NPT+1
      OPEN (EXT(JW),FILE=EXTFN(JW),STATUS='OLD')
      READ(EXT(JW),'(A1)')INFORMAT
      IF(INFORMAT=='$')EXTF(JW)=.TRUE.
      IF(EXTF(JW))THEN
      READ (EXT(JW),'(/)')
      READ (EXT(JW),*)NXEXT2(JW), EXTNX(JW)
      ELSE
      READ (EXT(JW),'(///2F8.0)') NXEXT2(JW), EXTNX(JW)
      ENDIF
      
      EXH2O(JW) = EXTNX(JW)
      EXTO(JW)  = EXTNX(JW)
      
      IF(EXTF(JW))THEN
      READ (EXT(JW),*)    NXEXT1(JW), EXTNX(JW)
          ELSE
      READ (EXT(JW),'(2F8.0)')    NXEXT1(JW), EXTNX(JW)
      ENDIF
    END IF
    !DO I=CUS(BS(JW)),DS(BE(JW))   ! SW CODE FIX 5-21-15
    !  WIND2(I) = WIND(JW)*WSC(I)*DLOG(2.0D0/Z0(JW))/DLOG(WINDH(JW)/Z0(JW))
    !END DO
  END DO
  IF (NWD > 0) THEN
    WDQ = NPT; NPT = NPT+1
    OPEN (WDQ,FILE=QWDFN,STATUS='OLD')
    READ(WDQ,'(A1)')INFORMAT
    IF(INFORMAT=='$')WDQF=.TRUE.
    IF(WDQF)THEN
    READ (WDQ,'(/)')
    READ (WDQ,*) NXQWD2,(QWDNX(JW),JW=1,NWD)
    DO JW=1,NWD
      QWD(JW)  = QWDNX(JW)
      QWDO(JW) = QWDNX(JW)
    END DO
    READ (WDQ,*)    NXQWD1,(QWDNX(JW),JW=1,NWD)
    ELSE
    READ (WDQ,'(//10F8.0:/(8X,9F8.0))') NXQWD2,(QWDNX(JW),JW=1,NWD)
    DO JW=1,NWD
      QWD(JW)  = QWDNX(JW)
      QWDO(JW) = QWDNX(JW)
    END DO
    READ (WDQ,'(10F8.0:/(8X,9F8.0))')    NXQWD1,(QWDNX(JW),JW=1,NWD)
    ENDIF
  END IF
  IF (TRIBUTARIES) THEN
    DO JT=1,NTR
      TRQ(JT) = NPT; NPT = NPT+1
      TRT(JT) = NPT; NPT = NPT+1
      OPEN (TRQ(JT),FILE=QTRFN(JT),STATUS='OLD')
      OPEN (TRT(JT),FILE=TTRFN(JT),STATUS='OLD')
      
       READ( TRQ(JT),'(A1)')INFORMAT
       IF(INFORMAT=='$')TRQF(JT)=.TRUE.
            IF(TRQF(JT))THEN
            READ (TRQ(JT),'(/)')
            READ (TRQ(JT),*) NXQTR2(JT),QTRNX(JT)
            ELSE
            READ (TRQ(JT),'(//2F8.0)') NXQTR2(JT),QTRNX(JT)
            ENDIF
      READ( TRT(JT),'(A1)')INFORMAT
       IF(INFORMAT=='$')TRTF(JT)=.TRUE.
            IF(TRTF(JT))THEN
            READ (TRT(JT),'(/)')
            READ (TRT(JT),*) NXTTR2(JT),TTRNX(JT)
            ELSE
            READ (TRT(JT),'(//2F8.0)') NXTTR2(JT),TTRNX(JT)
            ENDIF
      
 !     READ (TRQ(JT),'(///2F8.0)') NXQTR2(JT),QTRNX(JT)
 !     READ (TRT(JT),'(///2F8.0)') NXTTR2(JT),TTRNX(JT)
      IF (TRIB_CONST(JT)) THEN
        TRC(JT) = NPT; NPT = NPT+1
        OPEN (TRC(JT),FILE=CTRFN(JT),STATUS='OLD')
        READ( TRC(JT),'(A1)')INFORMAT
        IF(INFORMAT=='$')TRCF(JT)=.TRUE.
            IF(TRCF(JT))THEN
            READ (TRC(JT),'(/)')
            READ (TRC(JT),*) NXCTR2(JT),(CTRNX(TRCN(JAC,JT),JT),JAC=1,NACTR(JT))
            ELSE
            READ (TRC(JT),'(//1000F8.0)') NXCTR2(JT),(CTRNX(TRCN(JAC,JT),JT),JAC=1,NACTR(JT))
            ENDIF
      END IF
    END DO
    QTR(1:NTR)    = QTRNX(1:NTR)
    QTRO(1:NTR)   = QTRNX(1:NTR)
    TTR(1:NTR)    = TTRNX(1:NTR)
    TTRO(1:NTR)   = TTRNX(1:NTR)
    CTR(:,1:NTR)  = CTRNX(:,1:NTR)
    CTRO(:,1:NTR) = CTRNX(:,1:NTR)
    DO JT=1,NTR
      IF(TRQF(JT))THEN
        READ (TRQ(JT),*) NXQTR1(JT),QTRNX(JT)   ! cb 5/22/14
      ELSE
      READ (TRQ(JT),'(2F8.0)') NXQTR1(JT),QTRNX(JT)
      END IF
      IF(TRTF(JT))THEN
        READ (TRT(JT),*) NXTTR1(JT),TTRNX(JT)  ! cb 5/22/14
      ELSE
      READ (TRT(JT),'(2F8.0)') NXTTR1(JT),TTRNX(JT)
      END IF
      IF (TRIB_CONST(JT)) THEN
        IF(TRCF(JT))THEN
        READ (TRC(JT),*) NXCTR1(JT),(CTRNX(TRCN(JAC,JT),JT),JAC=1,NACTR(JT))
            ELSE
        READ (TRC(JT),'(1000F8.0)') NXCTR1(JT),(CTRNX(TRCN(JAC,JT),JT),JAC=1,NACTR(JT))
        ENDIF
      END IF
    END DO
  END IF
  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      IF (UP_FLOW(JB)) THEN
        IF (.NOT. INTERNAL_FLOW(JB) .AND. .NOT. DAM_INFLOW(JB)) THEN                                                  !TC 08/03/04 RA 1/13/06
          INQ(JB)  = NPT; NPT = NPT+1
          INFT(JB) = NPT; NPT = NPT+1
          OPEN (INQ(JB) ,FILE=QINFN(JB),STATUS='OLD')
          OPEN (INFT(JB),FILE=TINFN(JB),STATUS='OLD')
       READ( INQ(JB),'(A1)')INFORMAT
       IF(INFORMAT=='$')INQF(JB)=.TRUE.
            IF(INQF(JB))THEN
            READ (INQ(JB),'(/)')
            READ (INQ(JB),*) NXQIN2(JB),QINNX(JB)
            ELSE
            READ (INQ(JB), '(//2F8.0)') NXQIN2(JB),QINNX(JB)
            ENDIF
       READ( INFT(JB),'(A1)')INFORMAT
       IF(INFORMAT=='$')INTF(JB)=.TRUE.
            IF(INTF(JB))THEN
            READ (INFT(JB),'(/)')
            READ (INFT(JB),*) NXTIN2(JB),TINNX(JB)
            ELSE
            READ (INFT(JB),'(//2F8.0)') NXTIN2(JB),TINNX(JB)
            ENDIF      
          
     !     READ (INQ(JB), '(///2F8.0)') NXQIN2(JB),QINNX(JB)
     !     READ (INFT(JB),'(///2F8.0)') NXTIN2(JB),TINNX(JB)
          IF (INFLOW_CONST(JB)) THEN
            INC(JB) = NPT; NPT = NPT+1
            OPEN (INC(JB),FILE=CINFN(JB),STATUS='OLD')
            READ( INC(JB),'(A1)')INFORMAT
            IF(INFORMAT=='$')INCF(JB)=.TRUE.
            IF(INCF(JB))THEN
            READ (INC(JB),'(/)')
            READ (INC(JB),*) NXCIN2(JB),(CINNX(INCN(JC,JB),JB),JC=1,NACIN(JB))
            ELSE
            READ (INC(JB),'(//1000F8.0)') NXCIN2(JB),(CINNX(INCN(JC,JB),JB),JC=1,NACIN(JB))
            ENDIF
          END IF
        END IF
        QIN(JB)    = QINNX(JB)
        QIND(JB)   = QINNX(JB)
        QINO(JB)   = QINNX(JB)
        TIN(JB)    = TINNX(JB)
        TIND(JB)   = TINNX(JB)
        TINO(JB)   = TINNX(JB)
        CIN(:,JB)  = CINNX(:,JB)
        CIND(:,JB) = CINNX(:,JB)
        CINO(:,JB) = CINNX(:,JB)
        IF (.NOT. INTERNAL_FLOW(JB) .AND. .NOT. DAM_INFLOW(JB)) THEN                                                  !TC 08/03/04  RA 1/13/06
            IF(INQF(JB))THEN
            READ (INQ(JB),*) NXQIN1(JB),QINNX(JB)
            ELSE
          READ (INQ(JB), '(2F8.0)') NXQIN1(JB),QINNX(JB)
            ENDIF          
            IF(INTF(JB))THEN
            READ (INFT(JB),*) NXTIN1(JB),TINNX(JB)
            ELSE
          READ (INFT(JB),'(2F8.0)') NXTIN1(JB),TINNX(JB)
            ENDIF 

          IF (INFLOW_CONST(JB)) THEN
            IF(INCF(JB))THEN
            READ (INC(JB),*) NXCIN1(JB),(CINNX(INCN(JC,JB),JB),JC=1,NACIN(JB))  
                  ELSE
            READ (INC(JB),'(1000F8.0)') NXCIN1(JB),(CINNX(INCN(JC,JB),JB),JC=1,NACIN(JB))
            ENDIF
          END IF
        END IF
      END IF
      IF (DN_FLOW(JB)) THEN
        IF (NSTR(JB) > 0) THEN
          OTQ(JB) = NPT; NPT = NPT+1
          OPEN (OTQ(JB),FILE=QOTFN(JB),STATUS='OLD')
          READ (OTQ(JB),'(A1)')INFORMAT
          IF(INFORMAT == '$')OTQF(JB)=.TRUE.
           IF(OTQF(JB))THEN
            READ(OTQ(JB),'(/)')
            READ (OTQ(JB),*) NXQOT2(JB),(QSTRNX(JS,JB),JS=1,NSTR(JB))
            QSTR(:,JB)  = QSTRNX(:,JB)
            QSTRO(:,JB) = QSTRNX(:,JB)
            READ (OTQ(JB),*)NXQOT1(JB),(QSTRNX(JS,JB),JS=1,NSTR(JB))
           ELSE
            READ (OTQ(JB),'(//10F8.0:/(8X,9F8.0))') NXQOT2(JB),(QSTRNX(JS,JB),JS=1,NSTR(JB))
            QSTR(:,JB)  = QSTRNX(:,JB)
            QSTRO(:,JB) = QSTRNX(:,JB)
            READ (OTQ(JB),'(10F8.0:/(8X,9F8.0))')    NXQOT1(JB),(QSTRNX(JS,JB),JS=1,NSTR(JB))
           ENDIF
           
           IF(DYNSTRUC(JB) == '      ON')THEN
               ODYNS(JB)=NPT; NPT=NPT+1
               WRITE (SEGNUM,'(I0)') JB     
               SEGNUM = ADJUSTL(SEGNUM)
               L      = LEN_TRIM(SEGNUM)   
               OPEN (ODYNS(JB),FILE='dynselev'//SEGNUM(1:L)//'.npt',STATUS='OLD')
               READ (ODYNS(JB),'(A1)')INFORMAT
               IF(INFORMAT == '$')DYNEF(JB)=.TRUE.
               READ(ODYNS(JB),*)NJS
               DO J=1,NJS
               READ(ODYNS(JB),*)JJS(J)
               ENDDO
               IF(DYNEF(JB))THEN
               READ(ODYNS(JB),*)
               READ (ODYNS(JB),*)NXDYNS(JB),(ESTR(JJS(J),JB), J=1,NJS)
               READ (ODYNS(JB),*)NXDYNS(JB),(NXESTRT(JJS(J),JB), J=1,NJS)
               ELSE
               READ (ODYNS(JB),'(/10F8.0:/(8X,9F8.0))')NXDYNS(JB),(ESTR(JJS(J),JB), J=1,NJS)
               READ (ODYNS(JB),'(10F8.0:/(8X,9F8.0))')NXDYNS(JB),(NXESTRT(JJS(J),JB), J=1,NJS) 
               ENDIF
            ENDIF
        END IF
      END IF
      IF (PRECIPITATION(JW)) THEN
        PRE(JB) = NPT; NPT = NPT+1
        PRT(JB) = NPT; NPT = NPT+1
        OPEN (PRE(JB),FILE=PREFN(JB),STATUS='OLD')
        OPEN (PRT(JB),FILE=TPRFN(JB),STATUS='OLD')
        READ(PRE(JB),'(A1)')INFORMAT
        IF(INFORMAT=='$')PRQF(JB)=.TRUE.
            IF(PRQF(JB))THEN
            READ (PRE(JB),'(/)')
            READ (PRE(JB),*) NXPR2(JB), PRNX(JB)
            ELSE
            READ (PRE(JB), '(//2F8.0)') NXPR2(JB), PRNX(JB)
            ENDIF
        READ(PRT(JB),'(A1)')INFORMAT
        IF(INFORMAT=='$')PRTF(JB)=.TRUE.
            IF(PRTF(JB))THEN
            READ (PRT(JB),'(/)')
            READ (PRT(JB),*) NXTPR2(JB),TPRNX(JB)
            ELSE
            READ (PRT(JB),'(//2F8.0)') NXTPR2(JB),TPRNX(JB)
            ENDIF      
        
        !READ (PRE(JB),'(///2F8.0)') NXPR2(JB), PRNX(JB)
        !READ (PRT(JB),'(///2F8.0)') NXTPR2(JB),TPRNX(JB)
        IF (PRECIP_CONST(JB)) THEN
          PRC(JB) = NPT; NPT = NPT+1
          OPEN (PRC(JB),FILE=CPRFN(JB),STATUS='OLD')
          READ (PRC(JB),'(A1)')INFORMAT
          IF(INFORMAT=='$')PRCF(JB)=.TRUE.
          IF(PRCF(JB))THEN
            READ (PRC(JB),'(/)')
            READ (PRC(JB),*) NXCPR2(JB),(CPRNX(PRCN(JAC,JB),JB),JAC=1,NACPR(JB))
              ELSE
          READ (PRC(JB),'(//1000F8.0)') NXCPR2(JB),(CPRNX(PRCN(JAC,JB),JB),JAC=1,NACPR(JB))
          ENDIF
        END IF
        PR(JB)    = PRNX(JB)
        TPR(JB)   = TPRNX(JB)
        CPR(:,JB) = CPRNX(:,JB)
            IF(PRQF(JB))THEN
            READ (PRE(JB),*) NXPR1(JB), PRNX(JB)
            ELSE
            READ (PRE(JB), '(2F8.0)') NXPR1(JB), PRNX(JB)
            ENDIF
            
            IF(PRTF(JB))THEN
            READ (PRT(JB),*) NXTPR1(JB),TPRNX(JB)
            ELSE
            READ (PRT(JB),'(2F8.0)') NXTPR1(JB),TPRNX(JB)
            ENDIF  
        !READ (PRE(JB),'(2F8.0)') NXPR1(JB), PRNX(JB)
        !READ (PRT(JB),'(2F8.0)') NXTPR1(JB),TPRNX(JB)
        IF (PRECIP_CONST(JB)) THEN
         IF(PRCF(JB))THEN
          READ (PRC(JB),*) NXCPR1(JB),(CPRNX(PRCN(JAC,JB),JB),JAC=1,NACPR(JB))
             ELSE
          READ (PRC(JB),'(1000F8.0)') NXCPR1(JB),(CPRNX(PRCN(JAC,JB),JB),JAC=1,NACPR(JB))
          ENDIF
        END IF
      END IF
      IF (DIST_TRIBS(JB)) THEN
        DTQ(JB) = NPT; NPT = NPT+1
        DTT(JB) = NPT; NPT = NPT+1
        OPEN (DTQ(JB),FILE=QDTFN(JB),STATUS='OLD')
        OPEN (DTT(JB),FILE=TDTFN(JB),STATUS='OLD')
        
        READ( DTQ(JB),'(A1)')INFORMAT
       IF(INFORMAT=='$')DTQF(JB)=.TRUE.
            IF(DTQF(JB))THEN
            READ (DTQ(JB),'(/)')
            READ (DTQ(JB),*) NXQDT2(JB),QDTRNX(JB)
            ELSE
            READ (DTQ(JB),'(//2F8.0)') NXQDT2(JB),QDTRNX(JB)
            ENDIF
      READ( DTT(JB),'(A1)')INFORMAT
       IF(INFORMAT=='$')DTTF(JB)=.TRUE.
            IF(DTTF(JB))THEN
            READ (DTT(JB),'(/)')
            READ (DTT(JB),*) NXTDT2(JB),TDTRNX(JB)
            ELSE
            READ (DTT(JB),'(//2F8.0)') NXTDT2(JB),TDTRNX(JB)
            ENDIF
        
        
  !      READ (DTQ(JB),'(///2F8.0)') NXQDT2(JB),QDTRNX(JB)
  !      READ (DTT(JB),'(///2F8.0)') NXTDT2(JB),TDTRNX(JB)
        IF (DTRIB_CONST(JB)) THEN
          DTC(JB) = NPT; NPT = NPT+1
          OPEN (DTC(JB),FILE=CDTFN(JB),STATUS='OLD')
          READ( DTC(JB),'(A1)')INFORMAT
          IF(INFORMAT=='$')DTCF(JB)=.TRUE.
          IF(DTCF(JB))THEN
          READ(DTC(JB),'(/)')
          READ (DTC(JB),*) NXCDT2(JB),(CDTRNX(DTCN(JAC,JB),JB),JAC=1,NACDT(JB))
          ELSE
          READ (DTC(JB),'(//1000F8.0)') NXCDT2(JB),(CDTRNX(DTCN(JAC,JB),JB),JAC=1,NACDT(JB))
          ENDIF
        END IF
        QDTR(JB)    = QDTRNX(JB)
        QDTRO(JB)   = QDTRNX(JB)
        TDTR(JB)    = TDTRNX(JB)
        TDTRO(JB)   = TDTRNX(JB)
        CDTR(:,JB)  = CDTRNX(:,JB)
        CDTRO(:,JB) = CDTRNX(:,JB)
            IF(DTQF(JB))THEN
            READ (DTQ(JB),*) NXQDT1(JB),QDTRNX(JB)
            ELSE
        READ (DTQ(JB),'(2F8.0)') NXQDT1(JB),QDTRNX(JB)
            ENDIF
               
            IF(DTTF(JB))THEN
            READ (DTT(JB),*) NXTDT1(JB),TDTRNX(JB)
            ELSE
        READ (DTT(JB),'(2F8.0)') NXTDT1(JB),TDTRNX(JB)
            ENDIF

        IF (DTRIB_CONST(JB)) THEN
        IF(DTCF(JB))THEN
        READ (DTC(JB),*) NXCDT1(JB),(CDTRNX(DTCN(JAC,JB),JB),JAC=1,NACDT(JB))
        ELSE
        READ (DTC(JB),'(1000F8.0)') NXCDT1(JB),(CDTRNX(DTCN(JAC,JB),JB),JAC=1,NACDT(JB))
        ENDIF
        END IF
      END IF
      IF (UH_EXTERNAL(JB)) THEN
        UHE(JB) = NPT; NPT = NPT+1
        UHT(JB) = NPT; NPT = NPT+1
        OPEN (UHE(JB),FILE=EUHFN(JB),STATUS='OLD')
        OPEN (UHT(JB),FILE=TUHFN(JB),STATUS='OLD')
        
          READ( UHE(JB),'(A1)')INFORMAT
          IF(INFORMAT=='$')THEN
              EUHF(JB)=1
          ENDIF
          
          READ( UHT(JB),'(A2)')INFORMAT2
          IF(INFORMAT2=='$T')THEN
              TUHF(JB)=1
          ELSEIF(INFORMAT2(1:1)=='$')THEN
              TUHF(JB)=2
          ENDIF
          
          IF(EUHF(JB)>0)THEN
          READ(UHE(JB),'(/)')
          READ (UHE(JB),*)  NXEUH2(JB), ELUHNX(JB)
          ELSE
          READ (UHE(JB),'(//2F8.0)')NXEUH2(JB), ELUHNX(JB)
          ENDIF
          
          IF(TUHF(JB)==1)THEN
          READ(UHT(JB),'(/)')
          READ (UHT(JB),*)  NXTUH2(JB),XX(1)
          TUHNX(2:KB(US(JB)),JB)=XX(1)
          ELSEIF(TUHF(JB)==2)THEN
          READ(UHT(JB),'(/)')
          READ (UHT(JB),*)  NXTUH2(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
          ELSE
          READ (UHT(JB),'(//10F8.0:/(8X,9F8.0))')NXTUH2(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
          ENDIF
        
       ! READ (UHE(JB),'(///2F8.0)')              NXEUH2(JB), ELUHNX(JB)
       ! READ (UHT(JB),'(///10F8.0:/(8X,9F8.0))') NXTUH2(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
        IF (CONSTITUENTS) THEN
          UHC(JB) = NPT; NPT = NPT+1
          OPEN (UHC(JB),FILE=CUHFN(JB),STATUS='OLD')
          
          READ( UHC(JB),'(A2)')INFORMAT2
          IF(INFORMAT2=='$T')THEN
              CUHF(JB)=1
          ELSEIF(INFORMAT2(1:1)=='$')THEN
              CUHF(JB)=2
          ENDIF
          
          READ (UHC(JB),'(/)')
          
        !  READ (UHC(JB),'(//)')
          IF(CUHF(JB)==1)THEN
              READ (UHC(JB),*)  NXCUH2(JB),(XX(CN(JAC)),JAC=1,NAC)
              DO JAC=1,NAC
              CUHNX(2:KB(US(JB)),CN(JAC),JB)=XX(CN(JAC))
              ENDDO
              
          ELSE
              
          
          DO JAC=1,NAC
            !IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (UHC(JB),'(10F8.0:/(8X,9F8.0))') NXCUH2(JB),(CUHNX(K,CN(JAC),JB),     &
            !                                                  K=2,KB(US(JB)))
            
            IF(CUHF(JB)==2)THEN
            READ (UHC(JB),*)  NXCUH2(JB),(CUHNX(K,CN(JAC),JB),K=2,KB(US(JB)))
            ELSE
            READ (UHC(JB),'(10F8.0:/(8X,9F8.0))') NXCUH2(JB),(CUHNX(K,CN(JAC),JB),K=2,KB(US(JB)))
            ENDIF
            
            
          END DO
          END IF
          
          ENDIF
        ELUH(JB)     = ELUHNX(JB)
        ELUHO(JB)    = ELUHNX(JB)
        TUH(:,JB)    = TUHNX(:,JB)
        TUHO(:,JB)   = TUHNX(:,JB)
        CUH(:,:,JB)  = CUHNX(:,:,JB)
        CUHO(:,:,JB) = CUHNX(:,:,JB)
        
          IF(EUHF(JB)>0)THEN
          READ (UHE(JB),*)  NXEUH1(JB), ELUHNX(JB)
          ELSE
          READ (UHE(JB),'(2F8.0)')NXEUH1(JB), ELUHNX(JB)
          ENDIF
        
        !READ (UHE(JB),'(2F8.0)')              NXEUH1(JB), ELUHNX(JB)
          
          IF(TUHF(JB)==1)THEN
          READ (UHT(JB),*)  NXTUH1(JB),XX(1)
          TUHNX(2:KB(US(JB)),JB)=XX(1)
          ELSEIF(TUHF(JB)==2)THEN
          READ (UHT(JB),*)  NXTUH1(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
          ELSE
          READ (UHT(JB),'(10F8.0:/(8X,9F8.0))')NXTUH1(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
          ENDIF
          
       ! READ (UHT(JB),'(10F8.0:/(8X,9F8.0))') NXTUH1(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
        IF (CONSTITUENTS) THEN
            
           IF(CUHF(JB)==1)THEN
              READ (UHC(JB),*)  NXCUH1(JB),(XX(CN(JAC)),JAC=1,NAC)
              DO JAC=1,NAC
              CUHNX(2:KB(US(JB)),CN(JAC),JB)=XX(CN(JAC))
              ENDDO
              
          ELSE
              
          
          DO JAC=1,NAC
 
            IF(CUHF(JB)==2)THEN
            READ (UHC(JB),*)  NXCUH1(JB),(CUHNX(K,CN(JAC),JB),K=2,KB(US(JB)))
            ELSE
            READ (UHC(JB),'(10F8.0:/(8X,9F8.0))') NXCUH1(JB),(CUHNX(K,CN(JAC),JB),K=2,KB(US(JB)))
            ENDIF
            
          END DO
          END IF
            
            
          !DO JAC=1,NAC
          !  IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (UHC(JB),'(10F8.0:/(8X,9F8.0))') NXCUH1(JB),(CUHNX(K,CN(JAC),JB),     &
          !                                                    K=2,KB(US(JB)))
          !END DO
        END IF
      END IF
      IF (DH_EXTERNAL(JB)) THEN
        DHE(JB) = NPT; NPT = NPT+1
        DHT(JB) = NPT; NPT = NPT+1
        OPEN (DHE(JB),FILE=EDHFN(JB),STATUS='OLD')
        OPEN (DHT(JB),FILE=TDHFN(JB),STATUS='OLD')
        
         READ( DHE(JB),'(A1)')INFORMAT
          IF(INFORMAT=='$')THEN
              EDHF(JB)=1
          ENDIF
          
          READ( DHT(JB),'(A2)')INFORMAT2
          IF(INFORMAT2=='$T')THEN
              TDHF(JB)=1
          ELSEIF(INFORMAT2(1:1)=='$')THEN
              TDHF(JB)=2
          ENDIF
          
          IF(EDHF(JB)>0)THEN
          READ(DHE(JB),'(/)')
          READ (DHE(JB),*)  NXEDH2(JB), ELDHNX(JB)
          ELSE
          READ (DHE(JB),'(//2F8.0)')NXEDH2(JB), ELDHNX(JB)
          ENDIF
          
          IF(TDHF(JB)==1)THEN
          READ(DHT(JB),'(/)')
          READ (DHT(JB),*)  NXTDH2(JB),XX(1)
          TDHNX(2:KB(DS(JB)),JB)=XX(1)
          ELSEIF(TDHF(JB)==2)THEN
          READ(DHT(JB),'(/)')
          READ (DHT(JB),*)  NXTDH2(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
          ELSE
          READ (DHT(JB),'(//10F8.0:/(8X,9F8.0))')NXTDH2(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
          ENDIF     
        
        !READ (DHE(JB),'(///10F8.0)')             NXEDH2(JB),ELDHNX(JB)
        !READ (DHT(JB),'(///10F8.0:/(8X,9F8.0))') NXTDH2(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
        IF (CONSTITUENTS) THEN
          DHC(JB) = NPT; NPT = NPT+1
          OPEN (DHC(JB),FILE=CDHFN(JB),STATUS='OLD')
          
           READ( DHC(JB),'(A2)')INFORMAT2
          IF(INFORMAT2=='$T')THEN
              CDHF(JB)=1
          ELSEIF(INFORMAT2(1:1)=='$')THEN
              CDHF(JB)=2
          ENDIF
          
          READ (DHC(JB),'(/)')
          
          IF(CDHF(JB)==1)THEN
              READ (DHC(JB),*)  NXCDH2(JB),(XX(CN(JAC)),JAC=1,NAC)
              DO JAC=1,NAC
              CDHNX(2:KB(DS(JB)),CN(JAC),JB)=XX(CN(JAC)) 
              ENDDO
          ELSE
          DO JAC=1,NAC          
            IF(CDHF(JB)==2)THEN
            READ (DHC(JB),*)  NXCDH2(JB),(CDHNX(K,CN(JAC),JB),K=2,KB(DS(JB)))
            ELSE
            READ (DHC(JB),'(10F8.0:/(8X,9F8.0))') NXCDH2(JB),(CDHNX(K,CN(JAC),JB),K=2,KB(DS(JB)))
            ENDIF
          END DO
          END IF

          !READ (DHC(JB),'(//)')
          !DO JAC=1,NAC
          !  IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (DHC(JB),'(10F8.0:/(8X,9F8.0))') NXCDH2(JB),(CDHNX(K,CN(JAC),JB),     &
          !                                                    K=2,KB(DS(JB)))
          !END DO
        END IF
        ELDH(JB)     = ELDHNX(JB)
        ELDHO(JB)    = ELDHNX(JB)
        TDH(:,JB)    = TDHNX(:,JB)
        TDHO(:,JB)   = TDHNX(:,JB)
        CDH(:,:,JB)  = CDHNX(:,:,JB)
        CDHO(:,:,JB) = CDHNX(:,:,JB)
        
          IF(EDHF(JB)>0)THEN
          READ (DHE(JB),*)  NXEDH1(JB), ELDHNX(JB)
          ELSE
          READ (DHE(JB),'(2F8.0)')NXEDH1(JB), ELDHNX(JB)
          ENDIF
          
          IF(TDHF(JB)==1)THEN
          READ (DHT(JB),*)  NXTDH1(JB),XX(1)
          TDHNX(2:KB(DS(JB)),JB)=XX(1)
          ELSEIF(TDHF(JB)==2)THEN
          READ (DHT(JB),*)  NXTDH1(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
          ELSE
          READ (DHT(JB),'(10F8.0:/(8X,9F8.0))')NXTDH1(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
          ENDIF     
        
        !READ (DHE(JB),'(10F8.0)')             NXEDH1(JB),ELDHNX(JB)
        !READ (DHT(JB),'(10F8.0:/(8X,9F8.0))') NXTDH1(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
        IF (CONSTITUENTS) THEN
            
          IF(CDHF(JB)==1)THEN
              READ (DHC(JB),*)  NXCDH1(JB),(XX(CN(JAC)),JAC=1,NAC)
              DO JAC=1,NAC
              CDHNX(2:KB(DS(JB)),CN(JAC),JB)=XX(CN(JAC))      
              ENDDO
          ELSE
          DO JAC=1,NAC          
            IF(CDHF(JB)==2)THEN
            READ (DHC(JB),*)  NXCDH1(JB),(CDHNX(K,CN(JAC),JB),K=2,KB(DS(JB)))
            ELSE
            READ (DHC(JB),'(10F8.0:/(8X,9F8.0))') NXCDH1(JB),(CDHNX(K,CN(JAC),JB),K=2,KB(DS(JB)))
            ENDIF
          END DO
          END IF
            
          !DO JAC=1,NAC
          !  IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (DHC(JB),'(10F8.0:/(8X,9F8.0))') NXCDH1(JB),(CDHNX(K,CN(JAC),JB),     &
          !                                                    K=2,KB(DS(JB)))
          !END DO
        END IF
      END IF
    END DO
  END DO
  IF (GATES) THEN
    GTQ = NPT; NPT = NPT+1
    OPEN (GTQ,FILE=QGTFN,STATUS='OLD')
    READ( GTQ,'(A1)')INFORMAT
     IF(INFORMAT=='$')GATEF=.TRUE.
    ! Added code is to allow for computing flow based on BGT elevation - this can be a target water level in the reservoir - but allowing that flow to be removed from a different elevation
    !READ(GTQ,*)                           ! SW 2/25/11
     IF(GATEF)THEN
             READ(GTQ,'(A8)')GT2CHAR
    
    IF(GT2CHAR == 'EGT2ELEV')THEN
    REWIND(GTQ)
    READ(GTQ,*) 
    READ(GTQ,*)GT2CHAR,(EGT2(JG),JG=1,NGT)
    ENDIF
    
    READ(GTQ,*)
    READ (GTQ,*) NXQGT2,(BGTNX(JG),JG=1,NGT)                                         
    WHERE (DYNGTC == '     ZGT')
      EGT  = BGTNX
	  egto=bgtnx
      BGT  = 1.0
      G1GT = 1.0
      G2GT = 1.0
    ELSEWHERE
      BGT = BGTNX
	  bgto = bgtnx
    END WHERE
    READ (GTQ,*)    NXQGT,(BGTNX(JG),JG=1,NGT)

     ELSE
         
    READ(GTQ,'(A8)')GT2CHAR
    
    IF(GT2CHAR == 'EGT2ELEV')THEN
    REWIND(GTQ)
    READ(GTQ,*) 
    READ(GTQ,'(8X,1000F8.0)')(EGT2(JG),JG=1,NGT)
    ENDIF
    
    READ(GTQ,*)
    READ (GTQ,'(1000F8.0)') NXQGT2,(BGTNX(JG),JG=1,NGT)    
!    READ (GTQ,'(///1000F8.0)') NXQGT2,(BGTNX(JG),JG=1,NGT)                                       
    WHERE (DYNGTC == '     ZGT')
      EGT  = BGTNX
	  egto=bgtnx
      BGT  = 1.0
      G1GT = 1.0
      G2GT = 1.0
    ELSEWHERE
      BGT = BGTNX
	  bgto = bgtnx
    END WHERE
    READ (GTQ,'(1000F8.0)')    NXQGT,(BGTNX(JG),JG=1,NGT)
 ENDIF
     
  END IF
   IF (PIPES)THEN                                         ! SW 5/5/10
    IOPENPIPE=0
    do j=1,npi
     if(DYNPIPE(j) == '      ON')then
     iopenpipe=1
     exit
     endif
    enddo
     if(iopenpipe == 1)then
      PIPED = NPT; NPT = NPT+1
      OPEN (PIPED,FILE='dynpipe.npt',STATUS='OLD')
      READ (PIPED,'(///1000F8.0)') NXQPT,(BPNX(J),J=1,NPI)
      BP = BPNX
      READ (PIPED,'(1000F8.0)')    NXQPT,(BPNX(J),J=1,NPI)
     END IF
  ENDIF
  IF (PUMPS)THEN                                         
    DO J=1,NPU
     if(DYNPUMP(j) == '      ON')then
     WRITE (SEGNUM,'(I0)') J      
     SEGNUM = ADJUSTL(SEGNUM)
     L      = LEN_TRIM(SEGNUM)   
     PUMPD(J) = NPT; NPT = NPT+1  
     OPEN (PUMPD(J),FILE='dynpump'//SEGNUM(1:L)//'.npt',STATUS='OLD') 
     READ( PUMPD(J),'(A1)')INFORMAT
     IF(INFORMAT=='$')DYNPUMPF(J)=.TRUE.
     If(DYNPUMPF(J))then
        READ (PUMPD(J),'(/)')
        READ (PUMPD(J),*) NXPUMP(J),EPU2(J),EONPU2(J),EOFFPU2(J),QPU2(J)
        EPU(J)=EPU2(J)
        EONPU(J)=EONPU2(J)
        EOFFPU(J)=EOFFPU2(J)
        QPU(J)=QPU2(J)
        READ (PUMPD(J),*) NXPUMP(J),EPU2(J),EONPU2(J),EOFFPU2(J),QPU2(J)
     ELSE
        READ (PUMPD(J),'(//1000F8.0)') NXPUMP(J),EPU2(J),EONPU2(J),EOFFPU2(J),QPU2(J)
        EPU(J)=EPU2(J)
        EONPU(J)=EONPU2(J)
        EOFFPU(J)=EOFFPU2(J)
        QPU(J)=QPU2(J)
        READ (PUMPD(J),'(1000F8.0)') NXPUMP(J),EPU2(J),EONPU2(J),EOFFPU2(J),QPU2(J)
    ENDIF
      !READ (PUMPD(J),'(///1000F8.0)') NXPUMP(J),EPU2(J),EONPU2(J),EOFFPU2(J),QPU2(J)
      !  EPU(J)=EPU2(J)
      !  EONPU(J)=EONPU2(J)
      !  EOFFPU(J)=EOFFPU2(J)
      !  QPU(J)=QPU2(J)
      !READ (PUMPD(J),'(1000F8.0)') NXPUMP(J),EPU2(J),EONPU2(J),EOFFPU2(J),QPU2(J)
     END IF
    ENDDO
  ENDIF
  
  NOPEN         = NPT-1
  DYNAMIC_SHADE = SHADEI < 0
  NUNIT=NPT
RETURN

!***********************************************************************************************************************************
!**                                                  R E A D  I N P U T  D A T A                                                  **
!***********************************************************************************************************************************

ENTRY READ_INPUT_DATA (NXTVD)
  NXTVD = 1.0E10

! Meteorological data

  DO WHILE (JDAY >= NXWSC)
    WSC = WSCNX
    IF(WSHF)THEN
    READ (WSH,*) NXWSC,(WSCNX(I),I=1,IMX)
    ELSE
    READ (WSH,'(10F8.0:/(8X,9F8.0))') NXWSC,(WSCNX(I),I=1,IMX)
    ENDIF
  END DO
  DO JW=1,NWB
    DO WHILE (JDAY >= NXMET1(JW))
      TDEW(JW)   = TDEWNX(JW)
      TDEWO(JW)  = TDEWNX(JW)
      WIND(JW)   = WINDNX(JW)
      WINDO(JW)  = WINDNX(JW)
      PHI(JW)    = PHINX(JW)
      PHIO(JW)   = PHINX(JW)
      IF (PHISET > 0) PHI(JW)  = PHISET
      IF (PHISET > 0) PHIO(JW) = PHISET
      TAIR(JW)   = TAIRNX(JW)
      TAIRO(JW)  = TAIRNX(JW)
      CLOUD(JW)  = CLOUDNX(JW)
      CLOUDO(JW) = CLOUDNX(JW)
      NXMET2(JW) = NXMET1(JW)
      IF (READ_RADIATION(JW)) THEN
        SRON(JW)  = SRONX(JW)
        SROO(JW)  = SRON(JW)
        IF(METF(JW))THEN
        READ (MET(JW),*) NXMET1(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW),SRONX(JW)    
        ELSE
        READ (MET(JW),'(7F8.0)') NXMET1(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW),SRONX(JW)
        ENDIF
        SRONX(JW) = SRONX(JW)*REFL
      ELSE
        IF(METF(JW))THEN  
        READ (MET(JW),*) NXMET1(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW)  
        ELSE
        READ (MET(JW),'(6F8.0)') NXMET1(JW),TAIRNX(JW),TDEWNX(JW),WINDNX(JW),PHINX(JW),CLOUDNX(JW)
        ENDIF
      END IF
    END DO
    NXTVD = MIN(NXTVD,NXMET1(JW))
    IF (READ_EXTINCTION(JW)) THEN
      DO WHILE (JDAY >= NXEXT1(JW))
        EXH2O(JW)  = EXTNX(JW)
        EXTO(JW)   = EXTNX(JW)
        NXEXT2(JW) = NXEXT1(JW)
          IF(EXTF(JW))THEN
          READ (EXT(JW),*)    NXEXT1(JW), EXTNX(JW)
          ELSE
          READ (EXT(JW),'(2F8.0)') NXEXT1(JW),EXTNX(JW)
          ENDIF
      END DO
    END IF
    !DO I=CUS(BS(JW)),DS(BE(JW))   ! SW CODE FIX 5-21-15
    !  WIND2(I) = WIND(JW)*WSC(I)*DLOG(2.0D0/Z0(JW))/DLOG(WINDH(JW)/Z0(JW))    ! old value  z0 == 0.003
    !END DO
 END DO

! Withdrawals

  IF (NWD > 0) THEN
    DO WHILE (JDAY >= NXQWD1)
      NXQWD2 = NXQWD1
      DO JWD=1,NWD
        QWD(JWD)  = QWDNX(JWD)
        QWDO(JWD) = QWDNX(JWD)
      END DO
      IF(WDQF)THEN
      READ (WDQ,*) NXQWD1,(QWDNX(JWD),JWD=1,NWD)  
      ELSE
      READ (WDQ,'(10F8.0:/(8X,9F8.0))') NXQWD1,(QWDNX(JWD),JWD=1,NWD)
      ENDIF
    END DO
    NXTVD = MIN(NXTVD,NXQWD1)
  END IF

! Tributaries

  IF (TRIBUTARIES) THEN
    DO JT=1,NTR

!**** Inflow

      DO WHILE (JDAY >= NXQTR1(JT))
        QTR(JT)    = QTRNX(JT)
        QTRO(JT)   = QTRNX(JT)
        NXQTR2(JT) = NXQTR1(JT)
        
          IF(TRQF(JT))THEN
          READ (TRQ(JT),*) NXQTR1(JT),QTRNX(JT)
          ELSE
        READ (TRQ(JT),'(2F8.0)') NXQTR1(JT),QTRNX(JT)
          ENDIF
        
        
  !      READ (TRQ(JT),'(2F8.0)') NXQTR1(JT),QTRNX(JT)
      END DO
      NXTVD = MIN(NXTVD,NXQTR1(JT))

!**** Inflow temperatures

      IF (JDAY >= NXTTR1(JT)) THEN
        DO WHILE (JDAY >= NXTTR1(JT))
          TTR(JT)    = TTRNX(JT)
          TTRO(JT)   = TTRNX(JT)
          NXTTR2(JT) = NXTTR1(JT)
          
          IF(TRTF(JT))THEN
          READ (TRT(JT),*) NXTTR1(JT),TTRNX(JT)
          ELSE
          READ (TRT(JT),'(2F8.0)') NXTTR1(JT),TTRNX(JT)
          ENDIF
          
          
   !       READ (TRT(JT),'(2F8.0)') NXTTR1(JT),TTRNX(JT)
        END DO
      END IF
      NXTVD = MIN(NXTVD,NXTTR1(JT))

!**** Inflow constituent concentrations

      IF (TRIB_CONST(JT)) THEN
        DO WHILE (JDAY >= NXCTR1(JT))
          CTR(TRCN(1:NACTR(JT),JT),JT)  = CTRNX(TRCN(1:NACTR(JT),JT),JT)
          CTRO(TRCN(1:NACTR(JT),JT),JT) = CTRNX(TRCN(1:NACTR(JT),JT),JT)
          NXCTR2(JT)                    = NXCTR1(JT)
          IF(TRCF(JT))THEN
          READ (TRC(JT),*) NXCTR1(JT),(CTRNX(TRCN(JAC,JT),JT),JAC=1,NACTR(JT))
          ELSE
          READ (TRC(JT),'(1000F8.0)') NXCTR1(JT),(CTRNX(TRCN(JAC,JT),JT),JAC=1,NACTR(JT))
          ENDIF
        END DO
        NXTVD = MIN(NXTVD,NXCTR1(JT))
      END IF
    END DO
  END IF

! Branch related inputs

  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)

!**** Inflow

      IF (UP_FLOW(JB)) THEN
        IF (.NOT. INTERNAL_FLOW(JB) .AND. .NOT. DAM_INFLOW(JB)) THEN                                                  !TC 08/03/04 RA 1/13/06
          DO WHILE (JDAY >= NXQIN1(JB))
            QIND(JB)   = QINNX(JB)
            QINO(JB)   = QINNX(JB)
            NXQIN2(JB) = NXQIN1(JB)
            
          IF(INQF(JB))THEN
          READ (INQ(JB),*) NXQIN1(JB),QINNX(JB)
          ELSE
            READ (INQ(JB),'(2F8.0)') NXQIN1(JB),QINNX(JB)
          ENDIF 
 
        !    READ (INQ(JB),'(2F8.0)') NXQIN1(JB),QINNX(JB)
          END DO
          NXTVD = MIN(NXTVD,NXQIN1(JB))

!******** Inflow temperature

          DO WHILE (JDAY >= NXTIN1(JB))
            TIND(JB)   = TINNX(JB)
            TINO(JB)   = TINNX(JB)
            NXTIN2(JB) = NXTIN1(JB)
            
          IF(INTF(JB))THEN
          READ (INFT(JB),*) NXTIN1(JB),TINNX(JB)
          ELSE
            READ (INFT(JB),'(2F8.0)') NXTIN1(JB),TINNX(JB)
          ENDIF 
            
     !       READ (INFT(JB),'(2F8.0)') NXTIN1(JB),TINNX(JB)
          END DO
          NXTVD = MIN(NXTVD,NXTIN1(JB))

!******** Inflow constituent concentrations

          IF (INFLOW_CONST(JB)) THEN
            DO WHILE (JDAY >= NXCIN1(JB))
              CIND(INCN(1:NACIN(JB),JB),JB) = CINNX(INCN(1:NACIN(JB),JB),JB)
              CINO(INCN(1:NACIN(JB),JB),JB) = CINNX(INCN(1:NACIN(JB),JB),JB)
              NXCIN2(JB)                    = NXCIN1(JB)
              IF(INCF(JB))THEN
              READ (INC(JB),*) NXCIN1(JB),(CINNX(INCN(JAC,JB),JB),JAC=1,NACIN(JB))    
              ELSE
              READ (INC(JB),'(1000F8.0)') NXCIN1(JB),(CINNX(INCN(JAC,JB),JB),JAC=1,NACIN(JB))
              ENDIF
            END DO
            NXTVD = MIN(NXTVD,NXCIN1(JB))
          END IF
        END IF
      END IF

!**** Outflow

      IF (DN_FLOW(JB) .AND. NSTR(JB) > 0) THEN
        DO WHILE (JDAY >= NXQOT1(JB))
          QSTR(1:NSTR(JB),JB)  = QSTRNX(1:NSTR(JB),JB)
          QSTRO(1:NSTR(JB),JB) = QSTRNX(1:NSTR(JB),JB)
          NXQOT2(JB)           = NXQOT1(JB)
          IF(OTQF(JB))THEN
          READ (OTQ(JB),*) NXQOT1(JB),(QSTRNX(JS,JB),JS=1,NSTR(JB))
          ELSE
          READ (OTQ(JB),'(10F8.0:/(8X,9F8.0))') NXQOT1(JB),(QSTRNX(JS,JB),JS=1,NSTR(JB))
          ENDIF
        END DO
        IF(DYNSTRUC(JB) == '      ON')THEN
        DO WHILE (JDAY >= NXDYNS(JB))
            DO J=1,NJS
            ESTR(JJS(J),JB)=NXESTRT(JJS(J),JB)
            ENDDO
            IF(DYNEF(JB))THEN
            READ (ODYNS(JB),*)NXDYNS(JB),(NXESTRT(JJS(J),JB), J=1,NJS)
            ELSE
            READ (ODYNS(JB),'(10F8.0:/(8X,9F8.0))')NXDYNS(JB),(NXESTRT(JJS(J),JB), J=1,NJS)
            ENDIF
        ENDDO
        NXTVD = MIN(NXTVD,NXDYNS(JB))
        ENDIF
        NXTVD = MIN(NXTVD,NXQOT1(JB))
      END IF

!**** Distributed tributaries

      IF (DIST_TRIBS(JB)) THEN

!****** Inflow

        DO WHILE (JDAY >= NXQDT1(JB))
          QDTR(JB)   = QDTRNX(JB)
          QDTRO(JB)  = QDTRNX(JB)
          NXQDT2(JB) = NXQDT1(JB)
          
          IF(DTQF(JB))THEN
          READ (DTQ(JB),*) NXQDT1(JB),QDTRNX(JB)
          ELSE
          READ (DTQ(JB),'(2F8.0)') NXQDT1(JB),QDTRNX(JB)
          ENDIF        
          
  !        READ (DTQ(JB),'(2F8.0)') NXQDT1(JB),QDTRNX(JB)
        END DO
        NXTVD = MIN(NXTVD,NXQDT1(JB))

!****** Temperature

        DO WHILE (JDAY >= NXTDT1(JB))
          TDTR(JB)   = TDTRNX(JB)
          TDTRO(JB)  = TDTRNX(JB)
          NXTDT2(JB) = NXTDT1(JB)
          
          IF(DTTF(JB))THEN
          READ (DTT(JB),*) NXTDT1(JB),TDTRNX(JB)
          ELSE
          READ (DTT(JB),'(2F8.0)') NXTDT1(JB),TDTRNX(JB)
          ENDIF 
          
 !         READ (DTT(JB),'(2F8.0)') NXTDT1(JB),TDTRNX(JB)
        END DO
        NXTVD = MIN(NXTVD,NXTDT1(JB))

!****** Constituent concentrations

        IF (DTRIB_CONST(JB)) THEN
          DO WHILE (JDAY >= NXCDT1(JB))
            CDTR(DTCN(1:NACDT(JB),JB),JB)  = CDTRNX(DTCN(1:NACDT(JB),JB),JB)
            CDTRO(DTCN(1:NACDT(JB),JB),JB) = CDTRNX(DTCN(1:NACDT(JB),JB),JB)
            NXCDT2(JB)                     = NXCDT1(JB)
            IF(DTCF(JB))THEN
            READ (DTC(JB),*) NXCDT1(JB),(CDTRNX(DTCN(JAC,JB),JB),JAC=1,NACDT(JB))   
                ELSE
            READ (DTC(JB),'(1000F8.0)') NXCDT1(JB),(CDTRNX(DTCN(JAC,JB),JB),JAC=1,NACDT(JB))
            ENDIF
          END DO
          NXTVD = MIN(NXTVD,NXCDT1(JB))
        END IF
      END IF

!**** Precipitation

      IF (PRECIPITATION(JW)) THEN
        DO WHILE (JDAY >= NXPR1(JB))
          PR(JB)    = PRNX(JB)
          NXPR2(JB) = NXPR1(JB)
          
          IF(PRQF(JB))THEN
          READ (PRE(JB),*) NXPR1(JB),PRNX(JB)
          ELSE
          READ (PRE(JB),'(2F8.0)') NXPR1(JB),PRNX(JB)
          ENDIF        
          
          !READ (PRE(JB),'(2F8.0)') NXPR1(JB),PRNX(JB)
        END DO
        NXTVD = MIN(NXTVD,NXPR1(JB))

!****** Temperature

        DO WHILE (JDAY >= NXTPR1(JB))
          TPR(JB)    = TPRNX(JB)
          NXTPR2(JB) = NXTPR1(JB)
          
          IF(PRTF(JB))THEN
          READ (PRT(JB),*) NXTPR1(JB),TPRNX(JB)
          ELSE
          READ (PRT(JB),'(2F8.0)') NXTPR1(JB),TPRNX(JB)
          ENDIF        
          
          !READ (PRT(JB),'(2F8.0)') NXTPR1(JB),TPRNX(JB)
        END DO
        NXTVD = MIN(NXTVD,NXTPR1(JB))

!****** Constituent concentrations

        IF (PRECIP_CONST(JB)) THEN
          DO WHILE (JDAY >= NXCPR1(JB))
            CPR(PRCN(1:NACPR(JB),JB),JB) = CPRNX(PRCN(1:NACPR(JB),JB),JB)
            NXCPR2(JB)                   = NXCPR1(JB)
            IF(PRCF(JB))THEN
            READ (PRC(JB),*) NXCPR1(JB),(CPRNX(PRCN(JAC,JB),JB),JAC=1,NACPR(JB))   
                ELSE
            READ (PRC(JB),'(1000F8.0)') NXCPR1(JB),(CPRNX(PRCN(JAC,JB),JB),JAC=1,NACPR(JB))
            ENDIF
          END DO
          NXTVD = MIN(NXTVD,NXCPR1(JB))
        END IF
      END IF

!**** Upstream head conditions

      IF (UH_EXTERNAL(JB)) THEN

!****** Elevations

        DO WHILE (JDAY >= NXEUH1(JB))
          ELUH(JB)   = ELUHNX(JB)
          ELUHO(JB)  = ELUHNX(JB)
          NXEUH2(JB) = NXEUH1(JB)
          
          IF(EUHF(JB)>0)THEN
          READ (UHE(JB),*)  NXEUH1(JB), ELUHNX(JB)
          ELSE
          READ (UHE(JB),'(2F8.0)')NXEUH1(JB), ELUHNX(JB)
          ENDIF
          !READ (UHE(JB),'(2F8.0)') NXEUH1(JB),ELUHNX(JB)
        END DO
        NXTVD = MIN(NXTVD,NXEUH1(JB))

!****** Temperatures

        DO WHILE (JDAY >= NXTUH1(JB))
          DO K=2,KMX-1
            TUH(K,JB)  = TUHNX(K,JB)
            TUHO(K,JB) = TUHNX(K,JB)
          END DO
          NXTUH2(JB) = NXTUH1(JB)
          
          IF(TUHF(JB)==1)THEN
          READ (UHT(JB),*)  NXTUH1(JB),XX(1)
          TUHNX(2:KB(US(JB)),JB)=XX(1)
          ELSEIF(TUHF(JB)==2)THEN
          READ (UHT(JB),*)  NXTUH1(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
          ELSE
          READ (UHT(JB),'(10F8.0:/(8X,9F8.0))')NXTUH1(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
          ENDIF
          
          
          !READ (UHT(JB),'(10F8.0:/(8X,9F8.0))') NXTUH1(JB),(TUHNX(K,JB),K=2,KB(US(JB)))
        END DO
        NXTVD = MIN(NXTVD,NXTUH1(JB))

!****** Constituent concentrations

        IF (CONSTITUENTS) THEN
          DO WHILE (JDAY >= NXCUH1(JB))
            DO K=2,KMX-1
              CUH(K,CN(1:NAC),JB)  = CUHNX(K,CN(1:NAC),JB)
              CUHO(K,CN(1:NAC),JB) = CUHNX(K,CN(1:NAC),JB)
            END DO
            NXCUH2(JB) = NXCUH1(JB)
            
            IF(CUHF(JB)==1)THEN
              READ (UHC(JB),*)  NXCUH1(JB),(XX(CN(JAC)),JAC=1,NAC)
              DO JAC=1,NAC
              CUHNX(2:KB(US(JB)),CN(JAC),JB)=XX(CN(JAC))
              ENDDO
            ELSE
              
            DO JAC=1,NAC
            IF(CUHF(JB)==2)THEN
            READ (UHC(JB),*)  NXCUH1(JB),(CUHNX(K,CN(JAC),JB),K=2,KB(US(JB)))
            ELSE
            READ (UHC(JB),'(10F8.0:/(8X,9F8.0))') NXCUH1(JB),(CUHNX(K,CN(JAC),JB),K=2,KB(US(JB)))
            ENDIF
            END DO
           END IF
            !DO JAC=1,NAC
            !  IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (UHC(JB),'(10F8.0:/(8X,9F8.0))') NXCUH1(JB),(CUHNX(K,CN(JAC),JB),   &
            !                                                    K=2,KB(US(JB)))
            !END DO
          END DO
          NXTVD = MIN(NXTVD,NXCUH1(JB))
        END IF
      END IF

!**** Downstream head

      IF (DH_EXTERNAL(JB)) THEN

!****** Elevation

        DO WHILE (JDAY >= NXEDH1(JB))
          ELDH(JB)   = ELDHNX(JB)
          ELDHO(JB)  = ELDHNX(JB)
          NXEDH2(JB) = NXEDH1(JB)
          
          IF(EDHF(JB)>0)THEN
          READ (DHE(JB),*)  NXEDH1(JB), ELDHNX(JB)
          ELSE
          READ (DHE(JB),'(2F8.0)')NXEDH1(JB), ELDHNX(JB)
          ENDIF
          !READ (DHE(JB),'(2F8.0)') NXEDH1(JB),ELDHNX(JB)
        END DO
        NXTVD = MIN(NXTVD,NXEDH1(JB))

!****** Temperature

        DO WHILE (JDAY >= NXTDH1(JB))
          DO K=2,KMX-1
            TDH(K,JB)  = TDHNX(K,JB)
            TDHO(K,JB) = TDHNX(K,JB)
          END DO
          NXTDH2(JB) = NXTDH1(JB)
          IF(TDHF(JB)==1)THEN
          READ (DHT(JB),*)  NXTDH1(JB),XX(1)
          TDHNX(2:KB(DS(JB)),JB)=XX(1)
          ELSEIF(TDHF(JB)==2)THEN
          READ (DHT(JB),*)  NXTDH1(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
          ELSE
          READ (DHT(JB),'(10F8.0:/(8X,9F8.0))')NXTDH1(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
          ENDIF     
          !READ (DHT(JB),'(10F8.0:/(8X,9F8.0))') NXTDH1(JB),(TDHNX(K,JB),K=2,KB(DS(JB)))
        END DO
        NXTVD = MIN(NXTVD,NXTDH1(JB))

!****** Constituents

        IF (CONSTITUENTS) THEN
          DO WHILE (JDAY >= NXCDH1(JB))
            DO K=2,KMX-1
              CDH(K,CN(1:NAC),JB)  = CDHNX(K,CN(1:NAC),JB)
              CDHO(K,CN(1:NAC),JB) = CDHNX(K,CN(1:NAC),JB)
            END DO
            NXCDH2(JB) = NXCDH1(JB)
            IF(CDHF(JB)==1)THEN
              READ (DHC(JB),*)  NXCDH1(JB),(XX(CN(JAC)),JAC=1,NAC)
              DO JAC=1,NAC
              CDHNX(2:KB(DS(JB)),CN(JAC),JB)=XX(CN(JAC))    
              ENDDO
          ELSE
          DO JAC=1,NAC          
            IF(CDHF(JB)==2)THEN
            READ (DHC(JB),*)  NXCDH1(JB),(CDHNX(K,CN(JAC),JB),K=2,KB(DS(JB)))
            ELSE
            READ (DHC(JB),'(10F8.0:/(8X,9F8.0))') NXCDH1(JB),(CDHNX(K,CN(JAC),JB),K=2,KB(DS(JB)))
            ENDIF
          END DO
          END IF
            !DO JAC=1,NAC
            !  IF (ADJUSTL(CNAME2(CN(JAC))) /= 'AGE     ') READ (DHC(JB),'(10F8.0:/(8X,9F8.0))') NXCDH1(JB),(CDHNX(K,CN(JAC),JB),   &
            !                                                    K=2,KB(DS(JB)))
            !END DO
          END DO
          NXTVD = MIN(NXTVD,NXCDH1(JB))
        END IF
      END IF
    END DO
  END DO

! Gate height opening

  IF (GATES) THEN
    DO WHILE (JDAY >= NXQGT)
	  nxqgt2=nxqgt
      WHERE (DYNGTC == '     ZGT')
        EGT  = BGTNX
		egto = bgtnx
        BGT  = 1.0
        G1GT = 1.0
        G2GT = 1.0
      ELSEWHERE
        BGT = BGTNX
		bgto = bgtnx
      ENDWHERE
      IF(GATEF)THEN
                READ (GTQ,*) NXQGT,(BGTNX(JG),JG=1,NGT)  
      ELSE
                READ (GTQ,'(1000F8.0)') NXQGT,(BGTNX(JG),JG=1,NGT)
      ENDIF
    END DO
    NXTVD = MIN(NXTVD,NXQGT)
  END IF

! Pipe reduction factor

  IF (PIPES .and. iopenpipe == 1) THEN
    DO WHILE (JDAY >= NXQPT)
        BP = BPNX
      READ (PIPED,'(1000F8.0)') NXQPT,(BPNX(J),J=1,NPI)
    END DO
    NXTVD = MIN(NXTVD,NXQPT)
  END IF

! DYNAMIC PUMPS

 IF (PUMPS) THEN

 DO J=1,NPU
  IF(DYNPUMP(J) == '      ON')THEN
     DO WHILE (JDAY >= NXPUMP(J))
        EPU(J)=EPU2(J)
        EONPU(J)=EONPU2(J)
        EOFFPU(J)=EOFFPU2(J)
        QPU(J)=QPU2(J)
         IF(DYNPUMPF(J))then
        READ (PUMPD(J),*) NXPUMP(J),EPU2(J),EONPU2(J),EOFFPU2(J),QPU2(J)
        ELSE
        READ (PUMPD(J),'(1000F8.0)') NXPUMP(J),EPU2(J),EONPU2(J),EOFFPU2(J),QPU2(J)
        ENDIF
 !     READ (PUMPD(J),'(1000F8.0)') NXPUMP(J),EPU2(J),EONPU2(J),EOFFPU2(J),QPU2(J)
    END DO
    NXTVD = MIN(NXTVD,NXPUMP(J))
   ENDIF
  ENDDO
 ENDIF

! Dead sea case

  DO JW=1,NWB
    IF (NO_INFLOW(JW)) THEN
      QIN(BS(JW):BE(JW))    = 0.0
      QINO(BS(JW):BE(JW))   = 0.0
      QIND(BS(JW):BE(JW))   = 0.0
      QINNX(BS(JW):BE(JW))  = 0.0
      QDTR(BS(JW):BE(JW))   = 0.0
      QDTRO(BS(JW):BE(JW))  = 0.0
      QDTRNX(BS(JW):BE(JW)) = 0.0
      PR(BS(JW):BE(JW))     = 0.0
      PRNX(BS(JW):BE(JW))   = 0.0
    END IF
    IF (NO_OUTFLOW(JW)) THEN
      QSTR(:,BS(JW):BE(JW))   = 0.0
      QSTRO(:,BS(JW):BE(JW))  = 0.0
      QSTRNX(:,BS(JW):BE(JW)) = 0.0
    END IF
  END DO
  WHERE (NO_WIND)
    WIND   = 0.0
    WINDO  = 0.0
    WINDNX = 0.0
  ENDWHERE
  WHERE (READ_RADIATION .AND. NO_HEAT)
    SRON  = 0.0
    SROO  = 0.0
    SRONX = 0.0
  ENDWHERE
  IF (ANY(NO_INFLOW)) THEN
    QTR   = 0.0
    QTRO  = 0.0
    QTRNX = 0.0
    QWD   = 0.0
    QWDO  = 0.0
    QWDNX = 0.0
  END IF
RETURN

!***********************************************************************************************************************************
!**                                              I N T E R P O L A T E  I N P U T S                                               **
!***********************************************************************************************************************************

ENTRY INTERPOLATE_INPUTS

! Meteorological/light extinction data

  DO JW=1,NWB
    IF (INTERP_METEOROLOGY(JW)) THEN
      RATIO     = (NXMET1(JW)-JDAY)/(NXMET1(JW)-NXMET2(JW))
      TDEW(JW)  = (1.0-RATIO)*TDEWNX(JW)+RATIO*TDEWO(JW)
      WIND(JW)  = (1.0-RATIO)*WINDNX(JW)+RATIO*WINDO(JW)
      ! CONVERT PHIO AND PHINX TO LESS THAN 2*PI     SW 2/13/15
      DO WHILE(PHIO(JW)>2.*PI)
          PHIO(JW)=PHIO(JW)-2.*PI
      ENDDO
      DO WHILE(PHINX(JW)>2.*PI)
          PHINX(JW)=PHINX(JW)-2.*PI
      ENDDO     
      IF (PHIO(JW)-PHINX(JW) > PI) THEN                      
        PHI(JW) = (1.0-RATIO)*(PHINX(JW)+2.0*PI)+RATIO*PHIO(JW)
      ELSEIF (PHIO(JW)-PHINX(JW) < -PI) THEN                       ! WX 2/13/15
        PHI(JW) = (1.0-RATIO)*PHINX(JW)+RATIO*(PHIO(JW) +2.0*PI)   ! WX 2/13/15
      ELSE
        PHI(JW) = (1.0-RATIO)*PHINX(JW)+RATIO*PHIO(JW)
      END IF

      !IF (ABS(PHIO(JW)-PHINX(JW)) > PI) THEN
      !  PHI(JW) = (1.0-RATIO)*(PHINX(JW)+2.0*PI)+RATIO*PHIO(JW)
      !ELSE
      !  PHI(JW) = (1.0-RATIO)*PHINX(JW)+RATIO*PHIO(JW)
      !END IF
      TAIR(JW)  = (1.0-RATIO)*TAIRNX(JW) +RATIO*TAIRO(JW)
      CLOUD(JW) = (1.0-RATIO)*CLOUDNX(JW)+RATIO*CLOUDO(JW)
      IF (READ_RADIATION(JW)) SRON(JW) = (1.0-RATIO)*SRONX(JW)+RATIO*SROO(JW)
    END IF
    IF (READ_EXTINCTION(JW).AND.INTERP_EXTINCTION(JW)) THEN    ! 6/30/15 SW
      RATIO     = (NXEXT1(JW)-JDAY)/(NXEXT1(JW)-NXEXT2(JW))
      EXH2O(JW) = (1.0-RATIO)*EXTNX(JW)+RATIO*EXTO(JW)
    END IF
  END DO

! Withdrawals

  IF (NWD > 0) THEN
    QRATIO = (NXQWD1-JDAY)/(NXQWD1-NXQWD2)
    DO JWD=1,NWD
      IF (INTERP_WITHDRAWAL(JWD)) QWD(JWD) = (1.0-QRATIO)*QWDNX(JWD)+QRATIO*QWDO(JWD)
    END DO
  END IF

! Gates  adding interpolation cb 8/13/2010  
  IF (gates) THEN
    QRATIO = (NXQgt-JDAY)/(NXQgt-NXQgt2)
    DO Jg=1,ngt
      IF (INTERP_gate(Jg)) then
	    if (DYNGTC(Jg) .eq. '     ZGT')then
	      egt(jg) = (1.0-QRATIO)*bgtNX(jg)+QRATIO*egtO(jg)
		else
          bgt(jg) = (1.0-QRATIO)*bgtNX(Jg)+QRATIO*bgtO(Jg)
		end if   
	  end if
    END DO
  END IF

! Tributaries

  IF (NTR > 0) THEN
    DO JT=1,NTR
      IF (INTERP_TRIBS(JT)) THEN
        QRATIO = (NXQTR1(JT)-JDAY)/(NXQTR1(JT)-NXQTR2(JT))
        TRATIO = (NXTTR1(JT)-JDAY)/(NXTTR1(JT)-NXTTR2(JT))
        IF (TRIB_CONST(JT)) CRATIO = (NXCTR1(JT)-JDAY)/(NXCTR1(JT)-NXCTR2(JT))
        QTR(JT)                      = (1.0-QRATIO)*QTRNX(JT)                     +QRATIO*QTRO(JT)
        TTR(JT)                      = (1.0-TRATIO)*TTRNX(JT)                     +TRATIO*TTRO(JT)
        CTR(TRCN(1:NACTR(JT),JT),JT) = (1.0-CRATIO)*CTRNX(TRCN(1:NACTR(JT),JT),JT)+CRATIO*CTRO(TRCN(1:NACTR(JT),JT),JT)
      END IF
    END DO
  END IF

! Branch related inputs

  DO JB=1,NBR

!** Inflow

    IF (UP_FLOW(JB)) THEN
      IF (.NOT. INTERNAL_FLOW(JB) .AND. .NOT. DAM_INFLOW(JB)) THEN                                                    !TC 08/03/04 RA 1/13/06
        IF (INTERP_INFLOW(JB)) THEN
          QRATIO = (NXQIN1(JB)-JDAY)/(NXQIN1(JB)-NXQIN2(JB))
          TRATIO = (NXTIN1(JB)-JDAY)/(NXTIN1(JB)-NXTIN2(JB))
          IF (INFLOW_CONST(JB))  CRATIO = (NXCIN1(JB)-JDAY)/(NXCIN1(JB)-NXCIN2(JB))
          QIND(JB)                      = (1.0-QRATIO)*QINNX(JB)                  +QRATIO*QINO(JB)
          TIND(JB)                      = (1.0-TRATIO)*TINNX(JB)                  +TRATIO*TINO(JB)
          CIND(INCN(1:NACIN(JB),JB),JB) = (1.0-CRATIO)*CINNX(INCN(1:NACIN(JB),JB),JB)+CRATIO*CINO(INCN(1:NACIN(JB),JB),JB)
        END IF
      END IF
    END IF

!** Outflow

    IF (DN_FLOW(JB) .AND. NSTR(JB) > 0) THEN
      QRATIO = (NXQOT1(JB)-JDAY)/(NXQOT1(JB)-NXQOT2(JB))
      DO JS=1,NSTR(JB)
        IF (INTERP_OUTFLOW(JS,JB)) QSTR(JS,JB) = (1.0-QRATIO)*QSTRNX(JS,JB)+QRATIO*QSTRO(JS,JB)
      END DO
    END IF

!** Distributed tributaries

    IF (DIST_TRIBS(JB)) THEN
      IF (INTERP_DTRIBS(JB)) THEN
        QRATIO = (NXQDT1(JB)-JDAY)/(NXQDT1(JB)-NXQDT2(JB))
        TRATIO = (NXTDT1(JB)-JDAY)/(NXTDT1(JB)-NXTDT2(JB))
        IF (DTRIB_CONST(JB))   CRATIO = (NXCDT1(JB)-JDAY)/(NXCDT1(JB)-NXCDT2(JB))
        QDTR(JB)                      = (1.0-QRATIO)*QDTRNX(JB)                  +QRATIO*QDTRO(JB)
        TDTR(JB)                      = (1.0-TRATIO)*TDTRNX(JB)                  +TRATIO*TDTRO(JB)
        CDTR(DTCN(1:NACDT(JB),JB),JB) = (1.0-CRATIO)*CDTRNX(DTCN(1:NACDT(JB),JB),JB)+CRATIO*CDTRO(DTCN(1:NACDT(JB),JB),JB)
      END IF
    END IF

!** Upstream head

    IF (UH_EXTERNAL(JB)) THEN
      IF (INTERP_HEAD(JB)) THEN
        HRATIO   = (NXEUH1(JB)-JDAY)/(NXEUH1(JB)-NXEUH2(JB))
        TRATIO   = (NXTUH1(JB)-JDAY)/(NXTUH1(JB)-NXTUH2(JB))
        IF (CONSTITUENTS) CRATIO = (NXCUH1(JB)-JDAY)/(NXCUH1(JB)-NXCUH2(JB))
        ELUH(JB) = (1.0-HRATIO)*ELUHNX(JB)+HRATIO*ELUHO(JB)
        DO K=2,KMX-1
          TUH(K,JB)           = (1.0-TRATIO)*TUHNX(K,JB)          +TRATIO*TUHO(K,JB)
          CUH(K,CN(1:NAC),JB) = (1.0-CRATIO)*CUHNX(K,CN(1:NAC),JB)+CRATIO*CUHO(K,CN(1:NAC),JB)
        END DO
      END IF
    END IF

!** Downstream head

    IF (DH_EXTERNAL(JB)) THEN
      IF (INTERP_HEAD(JB)) THEN
        HRATIO = (NXEDH1(JB)-JDAY)/(NXEDH1(JB)-NXEDH2(JB))
        TRATIO = (NXTDH1(JB)-JDAY)/(NXTDH1(JB)-NXTDH2(JB))
        IF (CONSTITUENTS) CRATIO = (NXCDH1(JB)-JDAY)/(NXCDH1(JB)-NXCDH2(JB))
        ELDH(JB) = (1.0-HRATIO)*ELDHNX(JB)+HRATIO*ELDHO(JB)
        DO K=2,KMX-1
          TDH(K,JB)           = (1.0-TRATIO)*TDHNX(K,JB)          +TRATIO*TDHO(K,JB)
          CDH(K,CN(1:NAC),JB) = (1.0-CRATIO)*CDHNX(K,CN(1:NAC),JB)+CRATIO*CDHO(K,CN(1:NAC),JB)
        END DO
      END IF
    END IF
  END DO
RETURN
ENTRY DEALLOCATE_TIME_VARYING_DATA
  DEALLOCATE (NXQTR1, NXTTR1, NXCTR1, NXQIN1, NXTIN1, NXCIN1, NXQDT1, NXTDT1, NXCDT1, NXPR1,  NXTPR1, NXCPR1, NXEUH1, NXTUH1)
  DEALLOCATE (NXCUH1, NXEDH1, NXTDH1, NXCDH1, NXQOT1, NXMET1, NXQTR2, NXTTR2, NXCTR2, NXQIN2, NXTIN2, NXCIN2, NXQDT2, NXTDT2)
  DEALLOCATE (NXCDT2, NXPR2,  NXTPR2, NXCPR2, NXEUH2, NXTUH2, NXCUH2, NXEDH2, NXTDH2, NXCDH2, NXQOT2, NXMET2, WSCNX,DYNPUMPF)
  DEALLOCATE (QDTRO,  TDTRO,  ELUHO,  ELDHO,  QWDO,   QTRO,   TTRO,   QINO,   TINO,   QDTRNX, TDTRNX, PRNX,   TPRNX,  ELUHNX)
  DEALLOCATE (ELDHNX, QWDNX,  QTRNX,  TTRNX,  QINNX,  TINNX,  SROO,   TAIRO,  TDEWO,  CLOUDO, PHIO,   WINDO,  TAIRNX, BGTNX,BPNX)
  DEALLOCATE (TDEWNX, CLOUDNX,PHINX,  WINDNX, SRONX,  TRQ,    TRT,    TRC,    INQ,    DTQ,    PRE,    UHE,    DHE,    INFT)
  DEALLOCATE (DTT,    PRT,    UHT,    DHT,    INC,    DTC,    PRC,    UHC,    DHC,    OTQ,    MET,    EXT,    EXTNX,  EXTO, EXTF)
  DEALLOCATE (NXEXT1, NXEXT2, CTRO,   CINO,   QOUTO,  CDTRO,  TUHO,   TDHO,   QSTRO,  CTRNX,  CINNX,  QOUTNX, CDTRNX, CPRNX)
  DEALLOCATE (TUHNX, TDHNX,   QSTRNX, CUHO,   CDHO,   CUHNX,  CDHNX, INCF, DTCF, PRCF, METF, ODYNS, NXDYNS,NXESTRT,DYNEF,JJS, PRTF, PRQF)
  DEALLOCATE (EUHF,TUHF,CUHF,EDHF,TDHF,CDHF,XX)
  DEALLOCATE (INFLOW_CONST,   TRIB_CONST,     DTRIB_CONST,    PRECIP_CONST, PUMPD, NXPUMP,EPU2,EONPU2,EOFFPU2,QPU2, OTQF, TRCF,TRQF, TRTF,DTTF, DTQF, INQF, INTF)
RETURN
END SUBROUTINE TIME_VARYING_DATA
