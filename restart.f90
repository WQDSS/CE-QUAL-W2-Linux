!***********************************************************************************************************************************
!*                                       S U B R O U T I N E    R E S T A R T   O U T P U T                                       **
!***********************************************************************************************************************************

SUBROUTINE RESTART_OUTPUT (RSOFN)
  USE GLOBAL; USE SCREENC; USE RSTART; USE GDAYC; USE GEOMC; USE KINETIC, ONLY:EPM,EPD,SEDC,SEDN,SEDP,PH, SDKV; USE TVDC, ONLY:QSUM     
  USE KINETIC, ONLY:SED, PFLUXIN,NFLUXIN ; USE ZOOPLANKTONC, ONLY: ZOO; USE EDDY, ONLY: TKE; USE MAIN, ONLY:ENVIRPC;USE ENVIRPMOD; USE LOGICC;USE STRUCTURES
  USE MACROPHYTEC
  IMPLICIT NONE
  
  CHARACTER(*) :: RSOFN
  OPEN  (RSO,FILE=RSOFN,FORM='UNFORMATTED',STATUS='UNKNOWN')
  WRITE (RSO) NIT,    NV,     KMIN,   IMIN,   NSPRF,  CMBRT,  ZMIN,   IZMIN,  START,  CURRENT
  WRITE (RSO) DLTDP,  SNPDP,  TSRDP,  VPLDP,  PRFDP,  CPLDP,  SPRDP,  RSODP,  SCRDP,  FLXDP,  WDODP
  WRITE (RSO) JDAY,    ELTM,   ELTMF,  DLT,    DLTAV,  DLTS,   MINDLT, JDMIN,  CURMAX
  WRITE (RSO) NXTMSN, NXTMTS, NXTMPR, NXTMCP, NXTMVP, NXTMRS, NXTMSC, NXTMSP, NXTMFL, NXTMWD
  WRITE (RSO) VOLIN,  VOLOUT, VOLUH,  VOLDH,  VOLPR,  VOLTRB, VOLDT,  VOLWD,  VOLEV,  VOLSBR, VOLTR, VOLSR,VOLICE,ICEBANK
  WRITE (RSO) TSSEV,  TSSPR,  TSSTR,  TSSDT,  TSSWD,  TSSIN,  TSSOUT, TSSS,   TSSB,   TSSICE
  WRITE (RSO) TSSUH,  TSSDH,  TSSUH2, TSSDH2, CSSUH2, CSSDH2, VOLUH2, VOLDH2, QUH1
  WRITE (RSO) ESBR,   ETBR,   EBRI
  WRITE (RSO) Z,      SZ,     ELWS,   SAVH2,  SAVHR,  H2
  WRITE (RSO) KTWB,   KTI,    SKTI,   SBKT
  WRITE (RSO) ICE,    ICETH,  CUF,    QSUM
  WRITE (RSO) U,      W,      SU,     SW,     AZ,     SAZ,    DLTLIM
  WRITE (RSO) T1,     T2,     C1,     C2,     C1S,    SED,    KFS,    CSSK
  WRITE (RSO) EPD,   EPM
  WRITE (RSO) MACMBRT,MACRC,  SMACRC, MAC,    SMAC,   MACRM,  MACSS
  WRITE (RSO) SEDC, SEDN, SEDP, ZOO, CD  ! mlm 10/06
  WRITE(RSO) SDKV                        ! MLM 6/10/07
  WRITE(RSO) TKE                         ! SW 10/4/07
  WRITE(RSO) BR_INACTIVE                 ! SW 8/1/2018
  IF(ENVIRPC == '      ON')WRITE(RSO)T_CLASS,V_CLASS,C_CLASS,CD_CLASS,T_TOT,T_CNT,SUMVOLT,V_CNT,V_TOT,C_TOT,C_CNT,CD_TOT,CD_CNT
  IF(PIPES)WRITE(RSO)YS,VS,VST,YST,DTP,QOLD
  WRITE(RSO)TPOUT,TPTRIB,TPDTRIB,TPWD,TPPR,TPIN,TP_SEDSOD_PO4,PFLUXIN,TNOUT,TNTRIB,TNDTRIB,TNWD,TNPR,TNIN,TN_SEDSOD_NH4,NFLUXIN    !TP_SEDBURIAL,TN_SEDBURIAL,
  CLOSE (RSO)
END SUBROUTINE RESTART_OUTPUT