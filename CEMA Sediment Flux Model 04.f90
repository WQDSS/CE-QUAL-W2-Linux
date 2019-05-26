    RECURSIVE Subroutine CEMASedimentDiagenesis
	Use MAIN 
    Use GLOBAL
    Use GEOMC
    Use SCREENC
    Use RSTART
    Use PREC
    Use EDDY
    Use LOGICC
    Use TVDC
    Use KINETIC
    Use SURFHE
    Use CEMAVars
    Use TRANS
	!
    
    ! Type declarations
    IMPLICIT NONE
	!Local Variables
    Byte NFCle, NFLog
	Integer(2)	iTemp, RegnNum, iter, InitConRegn
	Real(8) SD_POC_L_Fr, SD_POC_R_Fr, SD_POC_I_Fr, SD_PON_L_Fr         
    Real(8) SD_PON_R_Fr, SD_PON_I_Fr, SD_PW_DiffCoeff, SD_Ox_Threshol      
    Real(8) SD_POP_L_Fr, SD_POP_R_Fr, SD_POP_I_Fr
    Real(8) SD_Ae_NH3_NO3_L, SD_Ae_NH3_NO3_H, SD_Ae_NO3_N2_L, SD_Ae_NO3_N2_H      
    Real(8) SD_An_NO3_N2, SD_Ae_CH4_CO2, SD_Ae_HS_NH4_Nit, SD_Ae_HS_O2_Nit     
    Real(8) SD_Theta_PW, SD_Theta_NH3_NO3, SD_Theta_NO3_N2, SD_Theta_CH4_CO2    
    Real(8) SD_Sulfate_CH4_H2S, SD_Ae_H2S_SO4, SD_Theta_H2S_SO4, SD_NormConst_H2S_SO4
    Real(8) SD_MinRate_PON_Lab, SD_MinRate_PON_Ref, SD_MinRate_PON_Ine, SD_MinRate_POC_Lab   
    Real(8) SD_MinRate_POC_Ref, SD_MinRate_POC_Ine, SD_Theta_PON_Lab, SD_Theta_PON_Ref    
    Real(8) SD_Theta_PON_Ine, SD_Theta_POC_Lab, SD_Theta_POC_Ref, SD_Theta_POC_Ine
    Real(8) SD_MinRate_POP_Lab, SD_MinRate_POP_Ref, SD_MinRate_POP_Ine
    Real(8) SD_Theta_POP_Lab, SD_Theta_POP_Ref, SD_Theta_POP_Ine
    Real(8) SD_tc, CellThickness, FastFluxCarbon, SlowFluxCarbon      !CellArea,   SW 10/3/2017
    Real(8) SD_Jcin, SD_Jpin, SD_Jnin, SD_O20, SD_JFeOOHin, SD_JMnO2in
    Real(8) SD_Depth, SD_Tw, SD_NH30, SD_NO30, SD_TIC0, SD_ALK0    
    !Real(8) SD_PO40, SD_CH40, SD_SO4, SD_SOD, SD_SO40
    Real(8) SD_PO40, SD_CH40, SD_SOD, SD_SO40  ! cb 7/26/18
    Real(8) SD_Fe20, SD_FeOOH0, SD_Mn20, SD_MnO20 
    Real(8) SD_JNH4, SD_JNO3, SD_JPO4, SD_JTIC, SD_JALK, SD_JFe2, SD_JFeOOH
    Real(8) SD_JMn2, SD_JT, SD_Ksw, SD_rhowcp, SD_tsed
    Real(8) SD_Ae_NH3_NO3, SD_Ae_NO3_N2
    Real(8) SD_BEN_STRp, SD_BEN_STRp2
    Real(8) SD_POCT2, SD_PONT2, SD_POPT2, SDPOCT1
    Real(8) SD_H1, SD_H2, SD_JC, SD_JN, SD_JP
    Real(8) maxit, SD_es, SD_KL12, SD_SODold, SD_JSOD
    Real(8) SD_CH4SAT, SD_s, SD_NH3toNO3, SD_fd1, SD_fp1, SD_fd2, SD_fp2
    Real(8) SD_Mn2toMnO2, SD_MnO2toMn2
    Real(8) SD_m1, SD_m2, SD_KdNH3, SD_a11, SD_a12, SD_b1, SD_a21, SD_a22, SD_b2
    Real(8) SD_w12, SD_w2, SD_nh3conv, SD_NSOD, SD_JDenitT, SD_JO2NO3T, SD_JC_O2equiv
    Real(8) SD_CSODmax, SD_SECH_ARG, SD_CSOD, SD_CH4toCO2, SD_KdH2S1, SD_KdH2S2
    Real(8) SD_xappd1, SD_xappp1, SD_k1h1d, SD_k1h1p, SD_k2h2d, SD_k2h2p, SD_F12
    Real(8) SD_F21, SD_xk1, SD_xk2, SD_KappaH2Sp1, SD_ea
    Real(8) SD_Rho, SD_Porosity
    Real(8) SD_SJCH4, SD_JCH4g, SD_JCH4, SD_JHS, SD_BEN_STR
    Real(8) SD_Ammonia, SD_Ammonium, SD_SulfiMinus, SD_Sulfide
    Real(8) MW_Constituent
    Real(8) Dissolved_CH4_Src, Dissolved_NH3_Src, Dissolved_H2S_Src
    Real(8) Dissolved_NO3_Src, Dissolved_SO4_Src, Dissolved_CO2_Src
    Real(8) Dissolved_ALK_Src, Dissolved_PO4_Src, Dissolved_Fe2_Src
    Real(8) Dissolved_Mn2_Src, Sediment_Heat_Src
    Real(8) Dissolved_O2_Snk, CO2ProducedSrc1L1, SO4ProducedSrc1L1
    Real(8) LPOM_Resuspension, RPOM_Resuspension, LPOMP_Resuspension, RPOMP_Resuspension
    Real(8) LPOMN_Resuspension, RPOMN_Resuspension
    Real(8) SO4ConsumedSnk1L2, CO2ProducedSrc2L2, SD1_Ammonia, SD1_Ammonium
	Real(8) SD2_Ammonia, SD2_Ammonium, SedTemp, VolWater, SedTemp1, SedTemp2
	Real(8) AmmoniaG_SD1, AmmoniaD_SD1, AmmoniaG_SD2, AmmoniaD_SD2
	Real(8) SD1_SulfiMinus, SD1_Sulfide, SD2_SulfiMinus, SD2_Sulfide
	Real(8) SulfideG_SD1, SulfideD_SD1, SulfideG_SD2, SulfideD_SD2
	Real(8) CO2ProducedCon1L1, CO2ProducedCon2L2, SD_JSO4
    Integer SD_CH4CompMethod, SD_POMResuspMethod
	Integer GenBODNum
	Real(8) Dissolved_BOD_Src(NumGenBODConstituents,2),SD_BODD(NumGenBODConstituents,2),SD_BOD
	Real(8) SedGenBODConc0(NumGenBODConstituents), PW_RelRate, PW_RelRate1, PW_RelRate2
    Real(8) SD_POM, SD_E
    REAL(8) FETCHW, U2, COEF1,COEF2,COEF3,COEF4,HS,TS,COEF,UORB,TAU, epsilon
    REAL(8) LW, LW0, LW1, Vscour, shields, molvisc_h2o
    REAL(8) SD_taubot, reyn_resusp, c_bottom, c_bottom2,CO2ProducedSrc1L2, SD_Fe2toFeOOH
    REAL(8) SD_FeOOHtoFe2, SD_Jctest, SD_JinFeOOH, SD_POCT1, SedimentHeat_Src
    
    REAL(8), ALLOCATABLE, DIMENSION(:) :: CellArea
    
    ! CEMA testing variables start
    Real(8) aO2n, o2ss,xnh4,z1ss,z2ss,z3ss,z4ss,z5ss,z5ass,z5bss,z6ss,xjn,acn,jnin,jpin
    ! CEMA testing end
    
    ! sediment pH start
    REAL(8) CART,ALKT,T1K,S2,SQRS2,DH1,DH2,H2CO3T,CO3T,PHT,F1,HION,HCO3T,AMMT,DH3,DHH
    REAL(8) KW,     INCR1,   OH,     K1,     K2, bicart, H2PO4T, H3PO4T, HPO4T, HT
    REAL(8) KAMM, KP1, KP2, KP3, NH3T, NH4T, OHT, OMCT, PHOST, PO4T
    Real(8) t1sed,ticsed,alksed,nh4sed,po4sed,pocsed,tdssed, phsed, CO2sed, HCO3sed, CO3sed
    INTEGER N,ITER1
    ! sediment pH end
                      
	SAVE
    ! 
    !  INPUTS 
    !  SD_Jcin = flux to sediments from settling organic carbon 
    !         from phytoplankton and detritus in oxygen equivalent units (gO2/m2/d) 
    !         (NOTE: gO2/m2/d = gC/m2/d * 2.67 gO2/gC)                            [SDINC] in gC/m2/d from WATER QUALITY SEDIMENTC Subroutine
    !  SD_Jnin = nitrogen flux in settling phytoplankton and detritus (gN/m2/d)   [SDINN]
    !  SD_Jpin = phosphorus flux in settling phytoplankton and detritus (gP/m2/d) [SDINP]
    !  SD_O20 = dissolved oxygen in water overlying the sediment (mgO2/L) 
    !  SD_depth = total water depth overlying the sediment (m) (used to calculate methane saturation concentration at in situ pressure)

    !  SD_Tw = temperature in water overlying the sediment (deg C) 
    !  SD_NH30 = ammonia N in water overlying the sediment (mgN/L) 
    !  SD_NO30 = nitrate N in water overlying the sediment (mgN/L) 
    !  SD_PO40 = soluble reactive P in water overlying the sediment (mgP/L) 
    !  SD_CH40 = fast reacting dissolved organic carbon and CBODu in the water overlying the sediment 
    !         in oxygen equivalent units (mgO2/L) 
    !         (NOTE: mgO2/L = mC/L * 2.67 mgO2/mgC) 
    !  SD_SALw = salinity in the water overlying the sediment (ppt) 
    ! 
    !  OUTPUTS 
    !  SOD = sediment oxygen demand flux of dissolved oxygen between the water and sediment (gO2/m2/d) 
    !        (positive is loss of O2 from water column) 
    !  Jnh4 = flux of ammonia N between the water and sediment (gN/m2/d) 
    !        (positive is source of NH4-N to water column) 
    !  Jno3 = flux of nitrate N between the water and sediment (gN/m2/d) 
    !        (positive is source of NO3-N to water column) 
    !  Jch4 = flux of dissolved methane, fast reacting C, and CBODu between water and sediment in O2 equivalent units (gO2/m2/d)
    !
    !        (positive is source of CBOD to water column) 
    !        (NOTE: gO2/m2/d = gC/m2/d * 2.67 gO2/gC) 
    !        (methane is not produced in salt water) 
    !  Jch4g = flux of methane gas bubbles between the water and sediment in O2 equivalent units (gO2/m2/d) 
    !        (positive is source of CH4 bubbles to water column) 
    !        (NOTE: gO2/m2/d = gC/m2/d * 2.67 gO2/gC) 
    !        (methane is not produced in salt water) 
    !  Jhs = flux of dissolved hydrogen sulfide (COD) between water and sediment in O2 equivalent units (gO2/m2/d) 
    !        (positive is source of COD to water column) 
    !        (hydrogen sulfide is not produced in freshwater) 
    !  Jpo4 = flux of soluble reactive P between the water and sedmiment (gP/m2/d) 
    !        (positive is source of PO4-P to water column) 
    !  NH3(1) and NH3(2) = ammonia N in the sediment layers 1 and 2 (mgN/L) 
    !  NO3(1) and NO3(2) = nitrate N in the sediment layers 1 and 2 (mgN/L) 
    !  CH4(1) = dissolved methane in the aerobic sediment layer 1 (O2 equivalent units mgO2/L) 
    !  HS(1) and HS(2) = dissolved sulfide in the sediment layers 1 and 2 (O2 equivalent units mgO2/L) 
    !  PO4(1) and PO4(2) = soluble reactive P in the sediment layers 1 and 2 (mgP/L) 
	! 
	!loop start for each control volume 
	!Set local default waterquality rates and constants
	!
	!Set initial time for averaging

	If(FirstTimeinCEMAMFTSedDiag)Then
	 
        ALLOCATE(CellArea(IMX))
        GasReleaseCH4=0.0   ! SW 10/10/2017
        !GasReleaseCO2=0.0   ! SW 10/19/2017
        ! determining width of sediment cells  ! cb 6/11/16 
        !Do JW=1, NWB            SW 9/18/2017
        !  Do JB=BS(JW),BE(JW)
        !    !IU = CUS(JB)
        !    IU = us(jb)
        !    ID = DS(JB)
        !    do SegNumI=iu,id
        !        sedcellwidth(SegNumI)=B(2,SegNumI)                      
        !    end do
        !  end do
        !end do
        CALL CEMAOutputRoutines   ! FIRST TIME CALL TO INITIALIZE VARIABLES  SW 3/9/16
	    !Get Rates and constants Region Numbers for each cell
	    CEMAMFT_RandC_RegN = 1
	    Do RegnNum = 1, NumRegnsSedimentDiagenesis
	        Do SegNumI = SedBedDiaRCRegSegSt(RegnNum), SedBedDiaRCRegSegEn(RegnNum)
	            CEMAMFT_RandC_RegN(SegNumI) = RegnNum
	        End Do !SegNumI
	    End Do !RegnNum
	    
	    !Get initial condition Region Numbers for each cell
	    CEMAMFT_InCond_RegN = 1
	    Do RegnNum = 1, NumRegnsSedimentBedComposition
	        Do SegNumI = SedBedInitRegSegSt(RegnNum), SedBedInitRegSegEn(RegnNum)
	            CEMAMFT_InCond_RegN(SegNumI) = RegnNum
	        End Do !SegNumI
	    End Do !RegnNum
	    	    
	    !Initialize variables
	    MFTSedFlxVars = 0.d00
	    Do JW=1, NWB
            KT = KTWB(JW)
            Do JB=BS(JW),BE(JW)
                !IU = CUS(JB)
                IU = us(jb)
                ID = DS(JB)
                Do SegNumI = IU, ID
                    CellArea(SegNumI) = B(KB(SegNumI),SegNumI)*DLX(SegNumI)     ! SW 10/3/2017
		            RegnNum = CEMAMFT_RandC_RegN(SegNumI)
		            InitConRegn = CEMAMFT_InCond_RegN(SegNumI)
    		        
		            Call CEMAMFTRatesandConstants
    		        
                    MFTSedFlxVars(SegNumI,16)   =   SD_POC_L_Fr*SDRegnPOC_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,17)   =   SD_POC_R_Fr*SDRegnPOC_T(InitConRegn)
	                MFTSedFlxVars(SegNumI,18)   =   SD_POC_I_Fr*SDRegnPOC_T(InitConRegn)
	                MFTSedFlxVars(SegNumI,19)   =   SD_PON_L_Fr*SDRegnPON_T(InitConRegn)
	                MFTSedFlxVars(SegNumI,20)   =   SD_PON_R_Fr*SDRegnPON_T(InitConRegn)
	                MFTSedFlxVars(SegNumI,21)   =   SD_PON_I_Fr*SDRegnPON_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,22)   =   SD_POP_L_Fr*SDRegnPOP_T(InitConRegn)
	                MFTSedFlxVars(SegNumI,23)   =   SD_POP_R_Fr*SDRegnPOP_T(InitConRegn)
	                MFTSedFlxVars(SegNumI,24)   =   SD_POP_I_Fr*SDRegnPOP_T(InitConRegn)
	                MFTSedFlxVars(SegNumI,25)   =   SDRegnSul_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,48)   =   SDRegnSul_T(InitConRegn)
	                !SD_SO4Conc(SegNumI)         =   SDRegnSul_T(InitConRegn)
	                !SD_pHValue(SegNumI)         =   CEMAPWpH
	                MFTSedFlxVars(SegNumI,7)    =   SDRegnCH4_T(InitConRegn)
	                !MFTSedFlxVars(SegNumI,5)    =   SDRegnH2S_T(InitConRegn)
	                !MFTSedFlxVars(SegNumI,6)    =   SDRegnH2S_T(InitConRegn)
	                !MFTSedFlxVars(SegNumI,14)   =   SDRegnNH3_T(InitConRegn)
	                !MFTSedFlxVars(SegNumI,15)   =   SDRegnNH3_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,5)   =   SDRegnNH3_T(InitConRegn)        ! cb 2/18/13
	                MFTSedFlxVars(SegNumI,6)   =   SDRegnNH3_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,14)    =   SDRegnH2S_T(InitConRegn)
	                MFTSedFlxVars(SegNumI,15)    =   SDRegnH2S_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,31)    =   SDRegnTIC_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,32)    =   SDRegnTIC_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,33)    =   SDRegnALK_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,34)    =   SDRegnALK_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,10)    =   SDRegnPO4_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,11)    =   SDRegnPO4_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,35)	=   CEMAPWpH
                    MFTSedFlxVars(SegNumI,36)   =   CEMAPWpH
                    
                    MFTSedFlxVars(SegNumI,37)    =   SDRegnFe2_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,38)    =   SDRegnFe2_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,39)    =   SDRegnFeOOH_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,40)    =   SDRegnFeOOH_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,41)    =   SDRegnMn2_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,42)    =   SDRegnMn2_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,43)    =   SDRegnMnO2_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,44)    =   SDRegnMnO2_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,45)    =   SDRegnT_T(InitConRegn)
                    MFTSedFlxVars(SegNumI,46)    =   SDRegnT_T(InitConRegn)
	                
                    Call UpdateCEMASedimentFluxVariables    ! cb 2/13/15
    	            
	            End Do !SegNumI
	        End Do !JB
        End Do !JW
        
        SD_w12=PartMixVel   ! Particle mixing velocity
        SD_w2=BurialVel     ! Burial Velocity
	    
	    !Generic BOD code
	    If(IncludeCEMAGenBODConstituents)Then
	        !Apply BOD constituent initial concentrations
	        Do i = 1, NumGenBODConstituents
	            Do RegnNum = 1, NumGenBODInitRegns
	                Do SegNumI = SedGenBODRegSegSt(RegnNum), SedGenBODRegSegEn(RegnNum)
	                    SedGenBODConc(i,SegNumI,:) = SedGenBODInit(i,RegnNum)
	                End Do !SegNum    
	            End Do !RegnNum
	        End Do !i  
	        
	        !Apply BOD constituent consumption rate
	        Do i = 1, NumGenBODConstituents
	            Do RegnNum = 1, NumGenBODConsumptionRegions
	                Do SegNumI = SedGenBODConsRegSegSt(RegnNum), SedGenBODConsRegSegEn(RegnNum)
	                    SedGenBODConsRate(i,SegNumI) = SedGenBODRegnRate(i,RegnNum)
	                    SedGenBODConsTCoeff(i,SegNumI) = SedGenBODRegnTCoeff(i,RegnNum)
	                End Do !SegNum    
	            End Do !RegnNum
	        End Do !i    
	    End If
	    !End generic BOD code        
	    
	    !Open Output Snapshot File
        If(WriteCEMAMFTSedFlx)Then
            Open(CEMASedFlxFilN1, File = "CEMADiagenesisConstituent.opt")
            Open(CEMASedFlxFilN2, File = "CEMADiagenesisSOD.opt")
            Open(CEMASedFlxFilN3, File = "CEMADiagenesisAerobicLayer.opt")
            Open(CEMASedFlxFilN4, File = "SOD.opt")
            Write(CEMASedFlxFilN4,'(" SOD (gO2/m2/d)")')
            Write(CEMASedFlxFilN4,'("    JDAY",<IMX>i10)')(SegNumI, SegNumI = 1, IMX)
            Open(CEMASedFlxFilN5, File = "POC_Sediments.opt")
            Write(CEMASedFlxFilN5,'(" POC (mg/l)")')
            Write(CEMASedFlxFilN5,'("    JDAY",<IMX>i10)')(SegNumI, SegNumI = 1, IMX)
            Open(CEMASedFlxFilN6, File = "PON_Sediments.opt")
            Write(CEMASedFlxFilN6,'(" PON (mg/l)")')
            Write(CEMASedFlxFilN6,'("    JDAY",<IMX>i10)')(SegNumI, SegNumI = 1, IMX)
            Open(CEMASedFlxFilN7, File = "POP_Sediments.opt")
            Write(CEMASedFlxFilN7,'(" POP (mg/l)")')
            Write(CEMASedFlxFilN7,'("    JDAY",<IMX>i10)')(SegNumI, SegNumI = 1, IMX)
            Open(CEMASedFlxFilN8, File = "CH4_Sediments.opt")
            Write(CEMASedFlxFilN8,'(" Dissolved CH4 (mg/l)")')
            Write(CEMASedFlxFilN8,'("    JDAY",<IMX>i10)')(SegNumI, SegNumI = 1, IMX)
            Open(CEMASedFlxFilN9, File = "H2S_Sediments.opt")
            Write(CEMASedFlxFilN9,'(" Dissolved H2S (mg/l)")')
            Write(CEMASedFlxFilN9,'("    JDAY",<IMX>i10)')(SegNumI, SegNumI = 1, IMX)
            ! testing start
            !open(4478,file='CEMA_testing_output.opt',status='unknown')
            !write(4478,5598)
!5598        format("    JDAY","      SOD2      SOD3      SOD4     CSOD2     CSOD3     CSOD4     NSOD2     NSOD3     NSOD4      POC2      POC3      POC4      PON2      PON3      PON4      SO42      SO43      SO44   TIC2_L1   TIC2_L2   TIC3_L1   TIC3_L2   TIC4_L1   TIC4_L2   ALK2_L1   ALK2_L2   ALK3_L1   ALK3_L2   ALK4_L1   ALK4_L2  PO4T2_L1  PO4T2_L2  PO4T3_L1  PO4T3_L2  PO4T4_L1  PO4T4_L2   PO42_L1   PO42_L2   PO43_L1   PO43_L2   PO44_L1   PO44_L2    PH2_L1    PH2_L2    PH3_L1    PH3_L2    PH4_L1    PH4_L2")            
!            open(4479,file='CEMA_testing_output_metals.opt',status='unknown')
!            write(4479,5599)
!5599        format("    JDAY","  FeII2_L1  FeII2_L2  FeII3_L1  FeII3_L2  FeII4_L1  FeII4_L2 FeOOH2_L1 FeOOH2_L2 FeOOH3_L1 FeOOH3_L2 FeOOH4_L1 FeOOH4_L2  MnII2_L1  MnII2_L2  MnII3_L1  MnII3_L2  MnII4_L1  MnII4_L2  MnO22_L1  MnO22_L2  MnO23_L1  MnO23_L2  MnO24_L1  MnO24_L2")            
!            open(4480,file='CEMA_testing_output_temperature.opt',status='unknown')
!            write(4480,5601)
!5601        format("    JDAY","     T2_L1     T2_L2     T3_L1     T3_L2     T4_L1     T4_L2")            
!            open(5387, file='resuspension_testing_wind.opt',status='unknown')
!            write(5387,5991)
!5991        format("        JDAY        NIT         SD_E      FetchW        Wind          U2          HS          TS          LW        UORB     EPSILON         TAU       TAUCR      DEPTHB       COEF1       COEF2       COEF3       COEF4")
!            open(5388, file='resuspension_testing_scour.opt',status='unknown')
!            write(5388,5992)
!5992        format("        JDAY        NIT         SD_E    c_bottom   c_bottom2      vscour     shields   crshields molvisc_h20 reyn_resusp")
!            ! testing end
            
        End If
        If(IncludeCEMAGenBODConstituents)Then
            Open(CEMAOutFilN5, File = "CEMAGenBOD.opt")
        End If
        
! moved from SetupCEMASedimentModel
!Open Output Snapshot File
        If(WriteBESnp .or. WriteBESnp)Then
            Open(CEMASNPOutFilN, File = "CEMAOutput.opt")
            !Open(CEMATSR1OutFilN, File = "CEMATSR1Output.opt")
            Write(CEMASNPOutFilN,'(a)')"$The file contains following output:"
            !Write(CEMATSR1OutFilN,'(a)')"$JDAY, KTWB(1), KB(15), z(15), BedElevation(15), BedElevationLayer(15), PorewaterRelRate(15), BedPorosity(15), BedConsolidRate(15)"
            If(WriteBESnp)Write(CEMASNPOutFilN,'(a)')"$Bed elevation (m)"
            If(WritePWSnp)Write(CEMASNPOutFilN,'(a)')"$Bed porosity (%)"
        End If
        
        !Bed Elevation
        If(WriteBESnp)Then
            Write(CEMASNPOutFilN,'("Bed elevation(m)")')
            Write(CEMASNPOutFilN,'("JDAY = ",f8.2)')JDAY
            Write(CEMASNPOutFilN,'(<IMX>i10)')(SegNumI, SegNumI = 1, IMX)
            Write(CEMASNPOutFilN,'(<IMX>f10.5)')(BedElevation(SegNumI), SegNumI = 1, IMX)
            
            Write(CEMASNPOutFilN,'("Bed elevation layer(m)")')
            Write(CEMASNPOutFilN,'("JDAY = ",f8.2)')JDAY
            Write(CEMASNPOutFilN,'(<IMX>i10)')(SegNumI, SegNumI = 1, IMX)
            Write(CEMASNPOutFilN,'(<IMX>f10.5)')(BedElevationLayer(SegNumI), SegNumI = 1, IMX)
        End If
            
        !Bed porosity
        If(WritePWSnp)Then
            Write(CEMASNPOutFilN,'("Bed porosity (%)")')
            Write(CEMASNPOutFilN,'("JDAY = ",f8.2)')JDAY
            Write(CEMASNPOutFilN,'(<IMX>i10)')(SegNumI, SegNumI = 1, IMX)
            Write(CEMASNPOutFilN,'(<IMX>f10.5)')(BedPorosity(SegNumI)*100, SegNumI = 1, IMX)
            Write(CEMASNPOutFilN,'("Total Volume of Sediments = ",e14.6," m3")')TotalSedimentsInBed
            Write(CEMASNPOutFilN,'("Total Porewater Volume = ",e14.6," m3")')TotalPoreWatVolume
            Write(CEMASNPOutFilN,'("Total Porewater Removed = ",e14.6," m3")')TotalPoreWatRemoved
        End If
        
        !Bottom Layer
        Open(CEMABtmLayFilN, File = "CEMABottomLayer.opt")
        Write(CEMABtmLayFilN,'(<IMX>i5)')(SegNumI, SegNumI = 1, IMX)
        TempCntr1 = 0
        EndBedConsolidation = .FALSE.
        CEMALayerAdded = .FALSE.
        CEMASSApplied = .FALSE.
        
        !Other files
        Open(CEMAOutFilN1, File = "CEMABubbles.opt")
        Open(CEMAOutFilN2, File = "CEMASedimentGas.opt")
        Open(CEMAOutFilN3, File = "CEMABubblesAtmosphereRelease.opt")
        Open(CEMAOutFilN6, File = "CEMASedDissGas.opt")
        Open(CEMAOutFilBub,File = "CEMABubbleReleaseSummary.csv")
        Write(CEMAOutFilN1,'(a)')"CEMA Sediment Bubbles Output"
        Write(CEMAOutFilN2,'(a)')"CEMA Sediment Bubbles Gas Output"
        Write(CEMAOutFilN3,'(a)')"CEMA Bubbles Release to Atmosphere Output"
        Write(CEMAOutFilN6,'(a)')"CEMA Sediment Dissolved Gas Output"
        Write(CEMAOutFilBub,*)'JDAY,WB,H2S(kg),CH4(kg),NH3(kg),CO2(kg),CH4(kg)-Method2'
        !Open(CEMAOutFilN4, File = "CEMABottomLayer.opt")
        


	     
	End If
	
	Do JW=1, NWB
        KT = KTWB(JW)
        Do JB=BS(JW),BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            Do SegNumI = IU, ID
            
                If(EndBedConsolidation(SegNumI))Then
                    Call EndSedimentFluxVariables
                End If
                
            End Do !SegNumI
        End Do !JB
    End Do !JW	
	
	SD_tc = dlt/86400.0     !Time step in days
	
	Do JW=1, NWB
        KT = KTWB(JW)
        Do JB=BS(JW),BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            Do SegNumI = IU, ID
                
                If(EndBedConsolidation(SegNumI))Cycle
		    
		        !CellArea = B(KB(SegNumI),SegNumI)*DLX(SegNumI)    ! SW 10/3/2017 MOVED TO INITIAL LOOP
                !CellArea =sedcellwidth(SegNumI)*DLX(SegNumI)
                !CellArea(SegNumI) = B(KB(SegNumI),SegNumI)*DLX(SegNumI)
                RegnNum = CEMAMFT_RandC_RegN(SegNumI)
                SD_H2 = BedElevation(SegNumI)
                CellThickness = H1(KB(SegNumI),SegNumI)
                
                Call CEMAMFTRatesandConstants
		
        		Call CEMAMFTSedimentFluxModel
                !if(SegNumI == 3)then  ! writing out resuspension variables for testing
                  !if(nit == 1 .or. nit/10*10 == nit)then
                 !     if(SD_POMResuspMethod == 0)then            
                 !       write(5387,'(f12.4,i12,30e12.4)')JDAY,nit,SD_E,FetchW,WIND(JW)*WSC(SegNumI),U2,HS,TS,LW,UORB,EPSILON,TAU,TAUCRPOM,DEPTHB(Kb(SegNumI),SegNumI),coef1,coef2,coef3,coef4
                 !     else
                 !       write(5388,'(f12.4,i12,30e12.4)')JDAY,nit,SD_E,c_bottom,c_bottom2,vscour,shields,crshield,smolvisc_h20,reyn_resusp
                 !     end if
                 ! end if
                !end if
        		
        	End Do !SegNumI
        End Do !JB
    End Do !JW
    
    If(FirstTimeinCEMAMFTSedDiag)FirstTimeinCEMAMFTSedDiag = .FALSE.
    If(FirstTimeInBubbles)FirstTimeInBubbles = .FALSE. 

	Return
    
    Entry CEMAMFTRatesandConstants
    
        SD_POC_L_Fr             =     SDRegnPOC_L_Fr(RegnNum)
        SD_POC_R_Fr             =     SDRegnPOC_R_Fr(RegnNum)
        SD_POC_I_Fr             =     1 - SD_POC_L_Fr - SD_POC_R_Fr
        SD_PON_L_Fr             =     SDRegnPON_L_Fr(RegnNum)
        SD_PON_R_Fr             =     SDRegnPON_R_Fr(RegnNum)
        SD_PON_I_Fr             =     1 - SD_PON_L_Fr - SD_PON_R_Fr
        SD_POP_L_Fr             =     SDRegnPOP_L_Fr(RegnNum)
        SD_POP_R_Fr             =     SDRegnPOP_R_Fr(RegnNum)
        SD_POP_I_Fr             =     1 - SD_POP_L_Fr - SD_POP_R_Fr
        SD_PW_DiffCoeff         =     SDRegnPW_DiffCoeff(RegnNum)    ! m^2/d
        SD_Ox_Threshol          =     SDRegnOx_Threshold(RegnNum)
        SD_Ae_NH3_NO3_L         =     SDRegnAe_NH3_NO3_L(RegnNum)
        SD_Ae_NH3_NO3_H         =     SDRegnAe_NH3_NO3_H(RegnNum)
        SD_Ae_NO3_N2_L          =     SDRegnAe_NO3_N2_L(RegnNum)
        SD_Ae_NO3_N2_H          =     SDRegnAe_NO3_N2_H(RegnNum)
        SD_An_NO3_N2            =     SDRegnAn_NO3_N2(RegnNum)
        SD_Ae_CH4_CO2           =     SDRegnAe_CH4_CO2(RegnNum)
        SD_Ae_HS_NH4_Nit        =     SDRegnAe_HS_NH4_Nit(RegnNum)
        SD_Ae_HS_O2_Nit         =     SDRegnAe_HS_O2_Nit(RegnNum)
        SD_Theta_PW             =     SDRegn_Theta_PW(RegnNum)   
        SD_Theta_NH3_NO3        =     SDRegn_Theta_NH3_NO3(RegnNum)   
        SD_Theta_NO3_N2         =     SDRegn_Theta_NO3_N2(RegnNum)  
        SD_Theta_CH4_CO2        =     SDRegn_Theta_CH4_CO2(RegnNum)
        SD_Sulfate_CH4_H2S      =     SDRegn_Sulfate_CH4_H2S(RegnNum)
        SD_Ae_H2S_SO4           =     SDRegnAe_H2S_SO4(RegnNum)
        SD_Theta_H2S_SO4        =     SDRegn_Theta_H2S_SO4(RegnNum)
        SD_NormConst_H2S_SO4    =     SDRegn_NormConst_H2S_SO4(RegnNum)
        SD_MinRate_PON_Lab      =     SDRegn_MinRate_PON_Lab(RegnNum)
        SD_MinRate_PON_Ref      =     SDRegn_MinRate_PON_Ref(RegnNum)
        SD_MinRate_PON_Ine      =     SDRegn_MinRate_PON_Ine(RegnNum)
        SD_MinRate_POC_Lab      =     SDRegn_MinRate_POC_Lab(RegnNum)
        SD_MinRate_POC_Ref      =     SDRegn_MinRate_POC_Ref(RegnNum)
        SD_MinRate_POC_Ine      =     SDRegn_MinRate_POC_Ine(RegnNum)
        SD_MinRate_POP_Lab      =     SDRegn_MinRate_POP_Lab(RegnNum)
        SD_MinRate_POP_Ref      =     SDRegn_MinRate_POP_Ref(RegnNum)
        SD_MinRate_POP_Ine      =     SDRegn_MinRate_POP_Ine(RegnNum)
        SD_Theta_PON_Lab        =     SDRegn_Theta_PON_Lab(RegnNum)
        SD_Theta_PON_Ref        =     SDRegn_Theta_PON_Ref(RegnNum)
        SD_Theta_PON_Ine        =     SDRegn_Theta_PON_Ine(RegnNum)
        SD_Theta_POC_Lab        =     SDRegn_Theta_POC_Lab(RegnNum)
        SD_Theta_POC_Ref        =     SDRegn_Theta_POC_Ref(RegnNum)
        SD_Theta_POC_Ine        =     SDRegn_Theta_POC_Ine(RegnNum)
        SD_Theta_POP_Lab        =     SDRegn_Theta_POP_Lab(RegnNum)
        SD_Theta_POP_Ref        =     SDRegn_Theta_POP_Ref(RegnNum)
        SD_Theta_POP_Ine        =     SDRegn_Theta_POP_Ine(RegnNum)
        SD_CH4CompMethod        =     SDRegn_CH4CompMethod(RegnNum)
        SD_POMResuspMethod      =     SDRegn_POMResuspMethod(RegnNum)
        
    Return

	Entry CEMAMFTSedimentFluxModel
	
	    Call UpdateCEMASedimentFluxVariables
		
		!All fluxes from water column are zero
		
		FastFluxCarbon	= 0.d0
	    SlowFluxCarbon	= 0.d0
        
        
 ! CEMA testing start
 
     !--Compute jn !(C) 
 !    xknh4= kappnh4 *thtanh4**(tempzz-20.) 

    !--making the conversion factor for jn to sod a function of temperature 
 !   ao2n = 106./16.*12./14.*32./12.*1.068**(tempzz-20.)
 !   o2ss = o2zz 

    !--convert from mg/m2-d to g/m2-d 
 !   xjnh4 = jnh4zz/1000.      ! jnh4zz is the ammonia flux to the sediments - this is the only term used to set all the other diagenesis terms
 !   xnh4= nh4zz/1000.         ! nh4zz is overlying water ammonia conc in mg/m3, this converts it to g/m3

    !--compute intermediate variables 
 !   z1ss =xjnh4**2                              ! xjnh4: rate of ammonia flux to overlying water
 !   z2ss =xknh4**2 
 !   z3ss =ao2n**2 
 !   z4ss = 1./z3ss 
 !   z5ss =o2ss**2 
 !   z5bss = (-4.*ao2n*o2ss*xknh4**4*xnh4**3 & 
 !      -z3ss*z1ss*z2ss*xnh4**2+18.*ao2n*o2ss*z1ss*z2ss*xnh4+ 27.*z5ss*z1ss*z2ss+4.*z3ss*xjnh4**4) !

    !--check for negative root. If this occurs, assume that the olw nh4=0 
  !  if(z5bss .le. 0.0) xnh4=0.0 
  !    z5ass= (z4ss*o2ss*xknh4*sqrt(-4.0*ao2n*o2ss*xknh4**4*xnh4**3-z3ss*z1ss*z2ss*xnh4**2+18.*ao2n*o2ss*z1ss*z2ss*xnh4+ &
  !    27.*z5ss*z1ss*z2ss+4.0*z3ss*xjnh4**4)/sqrt(3.)/6.0+z4ss*(9.*ao2n*o2ss*xjnh4*z2ss*xnh4+27.*z5ss*xjnh4*z2ss  &
  !    + 2.*z3ss*xjnh4**3)/54.0)                                ! Eq 16.5b
  !   z6ss = z5ass/abs(z5ass)*(abs(z5ass))**(1.0/3.0)                      ! Eq 16.5a
  !   xjn = z6ss+ (3.*o2ss*z2ss*xnh4+ao2n*z1ss)/(ao2n*z6ss)/9.0+xjnh4/3.0    ! xjn: rate at which organic matter is mineralized and ammonia is released to pore water 
  !   acn= 106./16.*12./14.                                       
  !   jcinzz= acn*xjn*1000.             !(D)   Convert ammonia flux to carbon diagenesis flux
     
  !   jnin= jcinzz/(frpon1+frpon2)/ratiocn
  !   jpin= jcinzz/(frpop1+frpop2)/ratiocp
   
 ! CEMA testing end
        
		!
		!Flux of caron
		!SD_Jcin         = (FastFluxCarbon + SlowFluxCarbon)*2.67 !(gC/m2/sec to gO2/m2/sec)
        !SD_Jcin         =  sdinc(kb(segnumi),segnumi)*2.67*h(kb(segnumi),jw)          !  cb 2/28/13 (gC/m3/sec to gO2/m2/sec)
        SD_Jcin=0.0
        !do k=ktwb(jw),kb(segnumi)
          SD_Jcin         =  sdinc(kb(segnumi),segnumi)*2.67*h(kb(segnumi),jw)     !+SD_Jcin     !  cb 2/28/13 (gC/m3/sec to gO2/m2/sec)
        !end do
        !SD_Jcin=jcinzz  ! CEMA testing
		!
		!Flux of nitrogen
		SD_Jnin			= 0.d0     
        !do k=ktwb(jw),kb(segnumi)
            SD_Jnin			= sdinn(kb(segnumi),segnumi)*h(kb(segnumi),jw) !+SD_Jnin
        !end do
	!SD_Jnin			= sdinn(kb(segnumi),segnumi)*h(kb(segnumi),jw)
    !    SD_Jnin=jnin   ! CEMA testing
		!
		!Flux of phosphorous
		SD_Jpin			= 0.d0
        !do k=ktwb(jw),kb(segnumi)  ! SW 8/31/2017
            SD_Jpin			= sdinp(kb(segnumi),segnumi)*h(kb(segnumi),jw) !+ SD_Jpin
        !end do
                !SD_Jpin			= sdinp(kb(segnumi),segnumi)*h(kb(segnumi),jw)
!                SD_Jpin=jpin
         SD_JFeOOHin=sdinFeOOH(kb(segnumi),segnumi)*h(kb(segnumi),jw)
         SD_JMnO2in=sdinMnO2(kb(segnumi),segnumi)*h(kb(segnumi),jw)
		
		SD_Jcin		=	SD_Jcin*86400.0d+00   ! conversions not needed for testing...
		SD_Jnin		=	SD_Jnin*86400.0d+00
		SD_Jpin		=	SD_Jpin*86400.0d+00
        SD_JinFeOOH	=	SD_JinFeOOH*86400.0d+00

		!Obtain properties from the water column
		SD_O20      =   MAX(O2(KB(SegNumI),SegNumI),0.0)         !c(m1(ii,jj),k,I_DO)
        !SD_O20=o2zz  ! testing
		SD_Depth    =   DEPTHB(KB(SegNumI),SegNumI)
		SD_Tw       =   MAX(T1(KB(SegNumI),SegNumI),0.0)         !c(m1(ii,jj),k,I_Temp)
!        SD_Tw       = tempzz
		SD_NH30     =   MAX(NH4(KB(SegNumI),SegNumI),0.0)        !c(m1(ii,jj),k,I_NH3)
!        SD_NH30     =nh4zz/1000.0   ! testing
		SD_NO30     =   MAX(NO3(KB(SegNumI),SegNumI),0.0)        !c(m1(ii,jj),k,I_NO3)
!        SD_NO30     =no3zz/1000.0
		SD_PO40     =   MAX(PO4(KB(SegNumI),SegNumI),0.0)        !c(m1(ii,jj),k,I_PO4) 
!        SD_PO40     =po4zz/1000.0
		
        SD_CH40     =   MAX(C1(KB(SegNumI),SegNumI,nch4),0.0)
        SD_SO40     =   MAX(C1(KB(SegNumI),SegNumI,nso4),0.0)
        sd_tic0     =   MAX(C1(KB(SegNumI),SegNumI,ntic),0.0)
        sd_alk0     =   MAX(C1(KB(SegNumI),SegNumI,nalk),0.0)
        sd_Fe20     =   MAX(C1(KB(SegNumI),SegNumI,nfe2),0.0)
        sd_FeOOH0   =   MAX(C1(KB(SegNumI),SegNumI,nfeOOH),0.0)
        sd_Mn20     =   MAX(C1(KB(SegNumI),SegNumI,nfe2),0.0)
        sd_MnO20    =   MAX(C1(KB(SegNumI),SegNumI,nMnO2),0.0)
        SD_Tsed     =   tsed(jw)
        SD_Ksw      =   cbhe(jw)
        SD_taubot   =   sb(KB(SegNumI), SegNumI)/b(kb(SegNumI),SegNumI)
		
        if(IncludeCEMAGenBODConstituents)then
		  Do GenBODNum = 1, NumGenBODConstituents
		    SedGenBODConc0(GenBODNum) = MAX(CBOD(KB(SegNumI),SegNumI,GenBODNum),0.0)
          End Do !GenBODNum
        end if
		
		!SD_SO4      =   SD_SO4Conc(SegNumI)             !GET SULFATE VALUE IN SEDIMENT
		
		If(SD_O20 <= 1.0e-10) Then
			SD_O20 = 1.0e-10
        End If
		
        SD_Rho=CEMASedimentDensity
        SD_Porosity=BedPorosity(SegnumI)
        
        Do iTemp = 1, 3			
			SD_POCT2 = SD_POCT2 + SD_POC22(iTemp)
			SD_PONT2 = SD_PONT2 + SD_PON22(iTemp)
			SD_POPT2 = SD_POPT2 + SD_POP22(iTemp)
		End Do
        
		!Use steady state approach to initialize the values
		Call DoCEMAMFTSedimentDiagen
		               
		!
		!Update source/sink terms
		DOSS(KB(SegNumI),SegNumI)       =     DOSS(KB(SegNumI),SegNumI)     - Dissolved_O2_Snk      !DO
        dosedia(KB(SegNumI),SegNumI) = - Dissolved_O2_Snk ! cb 3/19/13 storing for flux output
		NH4SS(KB(SegNumI),SegNumI)      =     NH4SS(KB(SegNumI),SegNumI)    + Dissolved_NH3_Src     !NH4
		NO3SS(KB(SegNumI),SegNumI)      =     NO3SS(KB(SegNumI),SegNumI)    + Dissolved_NO3_Src     !NO3
		CGSS(KB(SegNumI),SegNumI,ngh2s)     =     CGSS(KB(SegNumI),SegNumI,ngh2s)   + Dissolved_H2S_Src     !H2S
		CGSS(KB(SegNumI),SegNumI,ngch4)     =     CGSS(KB(SegNumI),SegNumI,ngch4)   + Dissolved_CH4_Src     !CH4
        CGSS(KB(SegNumI),SegNumI,ngso4)     =     CGSS(KB(SegNumI),SegNumI,ngso4)   + Dissolved_SO4_Src     !SO4
        CGSS(KB(SegNumI),SegNumI,ngFe2)     =     CGSS(KB(SegNumI),SegNumI,ngFe2)   + Dissolved_Fe2_Src     !Fe(II)
        CGSS(KB(SegNumI),SegNumI,ngMn2)     =     CGSS(KB(SegNumI),SegNumI,ngMn2)   + Dissolved_Mn2_Src     !Fe(II)
        TICSS(KB(SegNumI),SegNumI)     =     TICSS(KB(SegNumI),SegNumI)   + Dissolved_CO2_Src     !CO2
        ALKSS(KB(SegNumI),SegNumI)     =     ALKSS(KB(SegNumI),SegNumI)   + Dissolved_ALK_Src     !ALkalinity
        PO4SS(KB(SegNumI),SegNumI)     =     PO4SS(KB(SegNumI),SegNumI)   + Dissolved_PO4_Src     !PO4     
        SDPFLUX(JW)                        =     Dissolved_PO4_Src*DLT*VOL(KB(SEGNUMI),SEGNUMI)/1000.             + SDPFLUX(JW)      ! SW 8/31/2017 kg
        SDNH4FLUX(JW)                      =     Dissolved_NH3_Src*DLT*VOL(KB(SEGNUMI),SEGNUMI)/1000.             + SDNH4FLUX(JW)
        SDNO3FLUX(JW)                      =     Dissolved_NO3_Src*DLT*VOL(KB(SEGNUMI),SEGNUMI)/1000.             + SDNO3FLUX(JW)
        if(IncludeCEMAGenBODConstituents)then
          Do GenBODNum = 1, NumGenBODConstituents
		    CBODSS(KB(SegNumI),SegNumI,GenBODNum) = CBODSS(KB(SegNumI),SegNumI,GenBODNum) + SedGenBODSS(GenBODNum)  !BOD
          End Do !GenBODNum
        end if
        LPOMSS(KB(SegNumI),SegNumI)  =     LPOMSS(KB(SegNumI),SegNumI) + LPOM_Resuspension
        RPOMSS(KB(SegNumI),SegNumI)  =     RPOMSS(KB(SegNumI),SegNumI) + RPOM_Resuspension
        LPOMPSS(KB(SegNumI),SegNumI)  =     LPOMPSS(KB(SegNumI),SegNumI) + LPOMP_Resuspension
        RPOMPSS(KB(SegNumI),SegNumI)  =     RPOMPSS(KB(SegNumI),SegNumI) + RPOMP_Resuspension
        LPOMNSS(KB(SegNumI),SegNumI)  =     LPOMNSS(KB(SegNumI),SegNumI) + LPOMP_Resuspension
        RPOMNSS(KB(SegNumI),SegNumI)  =     RPOMNSS(KB(SegNumI),SegNumI) + RPOMP_Resuspension
        
        TSS(KB(SegNumI),SegNumI)     =     TSS(KB(SegNumI),SegNumI)   + SedimentHeat_Src     !Heat
        
        Call WriteCEMASedimentFluxVariables
		
	Return

	Entry UpdateCEMASedimentFluxVariables

		SD_NO3p2(1)		=	MFTSedFlxVars(SegNumI,1)			
		SD_NO3p2(2)		=	MFTSedFlxVars(SegNumI,2)	
		SD_NH3p2(1)		=	MFTSedFlxVars(SegNumI,3)	
		SD_NH3p2(2)		=	MFTSedFlxVars(SegNumI,4)	
		SD_NH3Tp2(1)	=	MFTSedFlxVars(SegNumI,5)	
		SD_NH3Tp2(2)	=	MFTSedFlxVars(SegNumI,6)	
		SD_CH4p2(1)		=	MFTSedFlxVars(SegNumI,7)	
		!SD_PO4p2(1)		=	MFTSedFlxVars(SegNumI,8)	
		!SD_PO4p2(2)		=	MFTSedFlxVars(SegNumI,9)
        SD_PO4(1)		=	MFTSedFlxVars(SegNumI,8)	
		SD_PO4(2)		=	MFTSedFlxVars(SegNumI,9)
		!SD_PO4Tp2(1)	=	MFTSedFlxVars(SegNumI,10)	
		!SD_PO4Tp2(2)	=	MFTSedFlxVars(SegNumI,11)	
        SD_PO4T(1)  	=	MFTSedFlxVars(SegNumI,10)	
		SD_PO4T(2)	   =	MFTSedFlxVars(SegNumI,11)	
		SD_HSp2(1)		=	MFTSedFlxVars(SegNumI,12)	
		SD_HSp2(2)		=	MFTSedFlxVars(SegNumI,13)	
		SD_HSTp2(1)		=	MFTSedFlxVars(SegNumI,14)	
		SD_HSTp2(2)		=	MFTSedFlxVars(SegNumI,15)	


		SD_POC22(1)		=	MFTSedFlxVars(SegNumI,16)	
		SD_POC22(2)		=	MFTSedFlxVars(SegNumI,17)
		SD_POC22(3)		=	MFTSedFlxVars(SegNumI,18)
		SD_PON22(1)		=	MFTSedFlxVars(SegNumI,19)
		SD_PON22(2)		=	MFTSedFlxVars(SegNumI,20)
		SD_PON22(3)		=	MFTSedFlxVars(SegNumI,21)
		SD_POP22(1)		=	MFTSedFlxVars(SegNumI,22)
		SD_POP22(2)		=	MFTSedFlxVars(SegNumI,23)
		SD_POP22(3)		=	MFTSedFlxVars(SegNumI,24)
		!SD_SO4          =   MFTSedFlxVars(SegNumI,25)
        SD_SO4(1)       =   MFTSedFlxVars(SegNumI,25)       ! cb 7/26/18
        SD_SO4(2)      =   MFTSedFlxVars(SegNumI,48)
		
		SD_JSOD         =   MFTSedFlxVars(SegNumI,26)
		SD_CSOD         =   MFTSedFlxVars(SegNumI,27)
		SD_NSOD         =   MFTSedFlxVars(SegNumI,28)
        
        SD_TIC(1)  		=   MFTSedFlxVars(SegNumI,31)	
        SD_TIC(2)  		=   MFTSedFlxVars(SegNumI,32)
        SD_ALK(1)  		=   MFTSedFlxVars(SegNumI,33)	
        SD_ALK(2)  		=   MFTSedFlxVars(SegNumI,34)
        SD_PH(1)  		=   MFTSedFlxVars(SegNumI,35)	
        SD_PH(2)  		=   MFTSedFlxVars(SegNumI,36)
        SD_Fe2T(1)      =   MFTSedFlxVars(SegNumI,37)	
        SD_Fe2T(2) 		=   MFTSedFlxVars(SegNumI,38)
        SD_FeOOH(1)		=   MFTSedFlxVars(SegNumI,39)	
        SD_FeOOH(2)		=   MFTSedFlxVars(SegNumI,40)
        SD_Mn2T(1)      =   MFTSedFlxVars(SegNumI,41)	
        SD_Mn2T(2) 		=   MFTSedFlxVars(SegNumI,42)
        SD_MnO2(1)		=   MFTSedFlxVars(SegNumI,43)	
        SD_MnO2(2)		=   MFTSedFlxVars(SegNumI,44)
        SD_T(1)	    	=   MFTSedFlxVars(SegNumI,45)	
        SD_T(2)		    =   MFTSedFlxVars(SegNumI,46)

	Return
	
	Entry WriteCEMASedimentFluxVariables
		!
		!  NH3(1) and NH3(2) = ammonia N in the sediment layers 1 and 2 (mgN/L) 
		!  NO3(1) and NO3(2) = nitrate N in the sediment layers 1 and 2 (mgN/L) 
		!  CH4(1) = dissolved methane in the aerobic sediment layer 1 (O2 equivalent units mgO2/L) 
		!  HS(1) and HS(2) = dissolved sulfide in the sediment layers 1 and 2 (O2 equivalent units mgO2/L) 
		!  PO4(1) and PO4(2) = soluble reactive P in the sediment layers 1 and 2 (mgP/L) 
		MFTSedFlxVars(SegNumI,1)		    =   SD_NO3p2(1)	
		MFTSedFlxVars(SegNumI,2)		    =   SD_NO3p2(2)	
		MFTSedFlxVars(SegNumI,3)		    =   SD_NH3p2(1)	
		MFTSedFlxVars(SegNumI,4)		    =   SD_NH3p2(2)	
		MFTSedFlxVars(SegNumI,5)		    =   SD_NH3Tp2(1)	
		MFTSedFlxVars(SegNumI,6)		    =   SD_NH3Tp2(2)	
		MFTSedFlxVars(SegNumI,7)		    =   SD_CH4p2(1)		
		MFTSedFlxVars(SegNumI,8)		    =   SD_PO4(1)	
		MFTSedFlxVars(SegNumI,9)		    =   SD_PO4(2)	
		MFTSedFlxVars(SegNumI,10)	        =   SD_PO4T(1)	
		MFTSedFlxVars(SegNumI,11)	        =   SD_PO4T(2)	
		MFTSedFlxVars(SegNumI,12)	        =   SD_HSp2(1)	
		MFTSedFlxVars(SegNumI,13)	        =   SD_HSp2(2)	
		MFTSedFlxVars(SegNumI,14)	        =   SD_HSTp2(1)	
		MFTSedFlxVars(SegNumI,15)	        =   SD_HSTp2(2)	

		MFTSedFlxVars(SegNumI,16)           =   SD_POC22(1)	
		MFTSedFlxVars(SegNumI,17)           =   SD_POC22(2)	
		MFTSedFlxVars(SegNumI,18)           =   SD_POC22(3)	
		MFTSedFlxVars(SegNumI,19)           =   SD_PON22(1)	
		MFTSedFlxVars(SegNumI,20)           =   SD_PON22(2)	
		MFTSedFlxVars(SegNumI,21)           =   SD_PON22(3)	
		MFTSedFlxVars(SegNumI,22)           =   SD_POP22(1)	
		MFTSedFlxVars(SegNumI,23)           =   SD_POP22(2)	
		MFTSedFlxVars(SegNumI,24)           =   SD_POP22(3)
		!MFTSedFlxVars(SegNumI,25)           =   SD_SO4	
        MFTSedFlxVars(SegNumI,25)           =   SD_SO4(1)	! cb 7/28/18
        MFTSedFlxVars(SegNumI,48)           =   SD_SO4(2)
		
		MFTSedFlxVars(SegNumI,26)           =   SD_JSOD
		MFTSedFlxVars(SegNumI,27)           =   SD_CSOD
		MFTSedFlxVars(SegNumI,28)           =   SD_NSOD
		
		MFTSedFlxVars(SegNumI,29)           =   SD_POC22(1)+ SD_POC22(2)+ SD_POC22(3)
		MFTSedFlxVars(SegNumI,30)           =   SD_PON22(1)+ SD_PON22(2)+ SD_PON22(3)
        MFTSedFlxVars(SegNumI,31)           =   SD_TIC(1)
        MFTSedFlxVars(SegNumI,32)           =   SD_TIC(2)
        MFTSedFlxVars(SegNumI,33)           =   SD_ALK(1)
        MFTSedFlxVars(SegNumI,34)           =   SD_ALK(2)
        MFTSedFlxVars(SegNumI,35)           =   SD_PH(1)
        MFTSedFlxVars(SegNumI,36)           =   SD_PH(2)
        MFTSedFlxVars(SegNumI,37)           =   SD_Fe2T(1)
        MFTSedFlxVars(SegNumI,38)           =   SD_Fe2T(2)
        MFTSedFlxVars(SegNumI,39)           =   SD_FeOOH(1)
        MFTSedFlxVars(SegNumI,40)           =   SD_FeOOH(2)
        MFTSedFlxVars(SegNumI,41)           =   SD_Mn2T(1)
        MFTSedFlxVars(SegNumI,42)           =   SD_Mn2T(2)
        MFTSedFlxVars(SegNumI,43)           =   SD_MnO2(1)
        MFTSedFlxVars(SegNumI,44)           =   SD_MnO2(2)
        MFTSedFlxVars(SegNumI,45)           =   SD_T(1)
        MFTSedFlxVars(SegNumI,46)           =   SD_T(2)
        MFTSedFlxVars(SegNumI,47)           =   SD_POP22(1)+ SD_POP22(2)+ SD_POP22(3)
		
		!For output
		
		CEMA_SD_Vars(SegNumI,1)             =   Dissolved_NH3_Src*CellThickness*86400.0
        CEMA_SD_Vars(SegNumI,2)             =   Dissolved_NO3_Src*CellThickness*86400.0
        CEMA_SD_Vars(SegNumI,3)             =   Dissolved_CH4_Src*CellThickness*86400.0
        CEMA_SD_Vars(SegNumI,4)             =   Dissolved_SO4_Src*CellThickness*86400.0
		CEMA_SD_Vars(SegNumI,5)             =   Dissolved_CO2_Src*CellThickness*86400.0
		CEMA_SD_Vars(SegNumI,6)             =   Dissolved_H2S_Src*CellThickness*86400.0
		CEMA_SD_Vars(SegNumI,7)             =   Dissolved_O2_Snk*CellThickness*86400.0
		
		CEMA_SD_Vars(SegNumI,8)             =   SD_POCT2
		CEMA_SD_Vars(SegNumI,9)             =   SD_PONT2   
		!CEMA_SD_Vars(SegNumI,10)             =   SD_SO4   
        CEMA_SD_Vars(SegNumI,10)             =   SD_SO4(1)       ! cb 7/26/18
        CEMA_SD_Vars(SegNumI,22)             =   SD_SO4(2)
	
	Return
	
	ENTRY DoCEMAMFTSedimentDiagen
	!Adapted from original work by Steve Chapra in November 2003 version of QUAL2K. 
    !Equation references are to DiToro 2001. Sediment Flux Modeling. Wiley-Interscience. 
    !Greg Pelletier re-organized the specification of the rate constants and 
    !added salinity-dependent nitrification/denitrification, inorganic P partitioning, 
    !and sulfide production/flux from carbon diagenesis based on the work of DiToro 2001. 
    !Significant modifications by Pelletier are marked by 'gp 
    ! 
    !  OUTPUTS 
    !  SOD = sediment oxygen demand flux of dissolved oxygen between the water and sediment (gO2/m2/d) 
    !        (positive is loss of O2 from water column) 
    !  Jnh4 = flux of ammonia N between the water and sediment (gN/m2/d) 
    !        (positive is source of NH4-N to water column) 
    !  Jno3 = flux of nitrate N between the water and sediment (gN/m2/d) 
    !        (positive is source of NO3-N to water column) 
    !  Jch4 = flux of dissolved methane, fast reacting C, and CBODu between water and sediment in O2 equivalent units (gO2/m2/d)
    !
    !        (positive is source of CBOD to water column) 
    !        (NOTE: gO2/m2/d = gC/m2/d * 2.67 gO2/gC) 
    !        (methane is not produced in salt water) 
    !  Jch4g = flux of methane gas bubbles between the water and sediment in O2 equivalent units (gO2/m2/d) 
    !        (positive is source of CH4 bubbles to water column) 
    !        (NOTE: gO2/m2/d = gC/m2/d * 2.67 gO2/gC) 
    !        (methane is not produced in salt water) 
    !  Jhs = flux of dissolved hydrogen sulfide (COD) between water and sediment in O2 equivalent units (gO2/m2/d) 
    !        (positive is source of COD to water column) 
    !        (hydrogen sulfide is not produced in freshwater) 
    !  Jpo4 = flux of soluble reactive P between the water and sedmiment (gP/m2/d) 
    !        (positive is source of PO4-P to water column) 
    !  NH3(1) and NH3(2) = ammonia N in the sediment layers 1 and 2 (mgN/L) 
    !  NO3(1) and NO3(2) = nitrate N in the sediment layers 1 and 2 (mgN/L) 
    !  CH4(1) = dissolved methane in the aerobic sediment layer 1 (O2 equivalent units mgO2/L) 
    !  HS(1) and HS(2) = dissolved sulfide in the sediment layers 1 and 2 (O2 equivalent units mgO2/L) 
    !  PO4(1) and PO4(2) = soluble reactive P in the sediment layers 1 and 2 (mgP/L) 
    ! 
    !  REFERENCES: 
    ! 
    !  Cerco, C.F. and T. Cole. 1995. User's guide to the CE-QUAL-ICM three-dimensional eutrophication model. 
    !  Release version 1.0. U.S. Army Corps of Engineers. Waterways Experiment Station. Vicksburg, MS. 
    ! 
    !  Chapra, S.C. and Pelletier, G.J. 2003. QUAL2K: A Modeling Framework for Simulating River and Stream 
    !  Water Quality (Beta Version): Documentation and Users Manual. 
    !  Civil and Environmental Engineering Dept., Tufts University, Medford, MA. 
    ! 
    !  DiToro, D.M. 2001. Sediment Flux Modeling. Wiley-Interscience. New York, NY. 
    !***** Assign rate constants (NOTE: these could be available for adjusting during model calibration) 
    !reaction velocities (m/d) 
		!
		!gp use salinity threshold to assign constants for the 
		!  nitrification/denitrification velocities for aerobic layer 1 
		!  and factor to increase aerobic layer partition for inorganic P 
		
		If(SD_O20 > SD_Ox_Threshol)Then
            SD_Ae_NH3_NO3   =   SD_Ae_NH3_NO3_H
            SD_Ae_NO3_N2    =   SD_Ae_NO3_N2_H
        Else
            SD_Ae_NH3_NO3   =   SD_Ae_NH3_NO3_L
            SD_Ae_NO3_N2    =   SD_Ae_NO3_N2_L
        End If
        
		!
		!gp        
		Do iTemp = 1, 3
			SD_poc2(iTemp) = SD_poc22(iTemp)
			SD_pon2(iTemp) = SD_pon22(iTemp)
			SD_pop2(iTemp) = SD_pop22(iTemp)           
		End Do
		Do iTemp = 1, 2
			SD_NH3Tp(iTemp) = SD_NH3Tp2(iTemp)
			SD_NO3p(iTemp)  = SD_NO3p2(iTemp)
		!	SD_PO4Tp(iTemp) = SD_PO4Tp2(iTemp)
			SD_HSTp(iTemp)  = SD_HSTp2(iTemp)
		End Do
		SD_BEN_STRp = SD_BEN_STRp2
		
		!gp assign constants for G class 1 and 2 PON, POC, and POP and calculate G class 3 as 1-fpox1-fpox2 
		SD_fpon(1) = SD_PON_L_Fr
		SD_fpon(2) = SD_PON_R_Fr
		SD_fpon(3) = SD_PON_I_Fr 
		SD_fpoc(1) = SD_POC_L_Fr
		SD_fpoc(2) = SD_POC_R_Fr
		SD_fpoc(3) = SD_POC_I_Fr
        SD_fpop(1) = SD_POP_L_Fr
		SD_fpop(2) = SD_POP_R_Fr
		SD_fpop(3) = SD_POP_I_Fr
		
		!
		!gp assign constants for G class 1, 2, and 3 mineralization of PON, POC, POP 
		SD_kdiaPON(1)   = SD_MinRate_PON_Lab
		SD_ThtaPON(1)   = SD_Theta_PON_Lab
		SD_kdiaPON(2)   = SD_MinRate_PON_Ref
		SD_ThtaPON(2)   = SD_Theta_PON_Ref
		SD_kdiaPON(3)   = SD_MinRate_PON_Ine
		SD_ThtaPON(3)   = SD_Theta_PON_Ine
		SD_kdiaPOC(1)   = SD_MinRate_POC_Lab
		SD_ThtaPOC(1)   = SD_Theta_POC_Lab
		SD_kdiaPOC(2)   = SD_MinRate_POC_Ref
		SD_ThtaPOC(2)   = SD_Theta_POC_Ref
		SD_kdiaPOC(3)   = SD_MinRate_POC_Ine
		SD_ThtaPOC(3)   = SD_Theta_POC_Ine
        SD_kdiaPOP(1)   = SD_MinRate_POP_Lab
		SD_ThtaPOP(1)   = SD_Theta_POP_Lab
		SD_kdiaPOP(2)   = SD_MinRate_POP_Ref
		SD_ThtaPOP(2)   = SD_Theta_POP_Ref
		SD_kdiaPOP(3)   = SD_MinRate_POP_Ine
		SD_ThtaPOP(3)   = SD_Theta_POP_Ine
		!
		!Compute input fluxes 
		Do iTemp = 1, 3 
			SD_JPOC(iTemp) = SD_Jcin * SD_fpoc(iTemp) 
			SD_JPON(iTemp) = SD_Jnin * SD_fpon(iTemp) 
			SD_JPOP(iTemp) = SD_Jpin * SD_fpop(iTemp) 
        End Do 
		!
        ! computing resuspended POM
        SD_POM=SD_POCT2/ORGC(JW)
        if(CEMA_POM_Resuspension_Processes)then
          if(SD_POMResuspMethod == 0)then
            Call CEMAWindInducedSedimentResuspension
          else
            Call CEMABottomScourResuspension
          end if
          
          SD_EPOC = 0.0
          SD_EPON = 0.0
          SD_EPOP = 0.0
          Do iTemp = 1, 3
            if( SD_POM > NONZERO .and. SD_POCT2 > Nonzero) SD_EPOC(itemp)= SD_POCT2/SD_POM * SD_POC2(itemp)/SD_POCT2 * SD_E            
            if( SD_POM > NONZERO .and. SD_PONT2 > Nonzero) SD_EPON(itemp)= SD_PONT2/SD_POM * SD_PON2(itemp)/SD_PONT2 * SD_E
            if( SD_POM > NONZERO .and. SD_POPT2 > Nonzero) SD_EPOP(itemp)= SD_POPT2/SD_POM * SD_POP2(itemp)/SD_POPT2 * SD_E            
          end do                        
        end if
        
		!Compute particulate organic forms 
		SD_POCT2 = 0.0
		SD_PONT2 = 0.0
		SD_POPT2 = 0.0
		
		!Equation 13.32 DiToro. See also Equation 12.2 - also includes resuspension
		Do iTemp = 1, 3
			SD_POC2(iTemp) = (SD_POC2(iTemp) + SD_JPOC(iTemp) * SD_tc / SD_H2 - SD_EPOC(iTemp) * SD_tc / SD_H2 ) /    &
                (1. + SD_kdiaPOC(iTemp) * SD_ThtaPOC(iTemp) ** (SD_T(2) - 20.) * SD_tc)
			SD_PON2(iTemp) = (SD_PON2(iTemp) + SD_JPON(iTemp) * SD_tc / SD_H2 - SD_EPON(iTemp) * SD_tc / SD_H2 ) /   &
                (1. + SD_kdiaPON(iTemp) * SD_ThtaPON(iTemp) ** (SD_T(2) - 20.) * SD_tc)
			SD_POP2(iTemp) = (SD_POP2(iTemp) + SD_JPOP(iTemp) * SD_tc / SD_H2 - SD_EPOP(iTemp) * SD_tc / SD_H2 ) /    &
                (1. + SD_kdiaPOP(iTemp) * SD_ThtaPOP(iTemp) ** (SD_T(2) - 20.) * SD_tc)
			SD_POCT2 = SD_POCT2 + SD_POC2(iTemp)
			SD_PONT2 = SD_PONT2 + SD_PON2(iTemp)
			SD_POPT2 = SD_POPT2 + SD_POP2(iTemp)
		End Do
	
		!
		!Compute diagenesis fluxes 
		!Equation 13.31 diagenesis term only. See also Equation 12.2, 12.5 and 12.6
		SD_Jc = 0
		SD_Jn = 0
		SD_Jp = 0 
		Do iTemp = 1,3 
			SD_Jc = SD_Jc + SD_kdiaPOC(iTemp) * SD_ThtaPOC(iTemp)**(SD_T(2) - 20.)*SD_POC2(iTemp)*SD_H2 
			SD_Jn = SD_Jn + SD_kdiaPON(iTemp) * SD_ThtaPON(iTemp)**(SD_T(2) - 20.)*SD_PON2(iTemp)*SD_H2 
			SD_Jp = SD_Jp + SD_kdiaPOP(iTemp) * SD_ThtaPOP(iTemp)**(SD_T(2) - 20.)*SD_POP2(iTemp)*SD_H2 
        End Do 
        SD_Jctest=SD_Jc

		maxit = 1000 
		SD_es = 0.001 
		
		!Calculation of final benthic particle mixing velocity Equation 13.6 in DoToro (2001)
		!See also Equation 2.56 and 4.48
		SD_KL12 = SD_PW_DiffCoeff * SD_Theta_PW ** (SD_T(2) - 20.) / (SD_H2 / 2.)
				
		If (FirstTimeinCEMAMFTSedDiag) Then
		    SD_SODold = SD_Jc + 1.714 * SD_Jn
		Else
		    SD_SODold = SD_JSOD
		End If
		
		iter = 0 
		!Saturation conc. of methane in oxygen equivalent units (Equation 10.51) 
		SD_CH4SAT = 100.0D+00*(1.0D+00 + SD_depth/10.0D+00)*(1.024**(20.0D+00 - SD_T(2)))        ![gmO*/m3] 

		Do While (.true.) 
			!CEMA Additional Code
			!Initialize variables to zero
			CO2ProducedSrc1L1 = 0.d0
			SO4ProducedSrc1L1 = 0.d0
			SO4ConsumedSnk1L2 = 0.d0
			CO2ProducedSrc2L2 = 0.d0
			!End CEMA Additional Code
			
			!Equation 3.15 or 10.36
            if(SD_SODold > nonzero)then
			  SD_s = SD_SODold/SD_O20 
            else
              SD_SODold = SD_Jc + 1.714 * SD_Jn   ! needed when segments are added
              SD_s = SD_SODold/SD_O20
            end if
            if(SD_s < nonzero) SD_s=nonzero
			!Equation 4.51
			SD_NH3toNO3 = SD_Ae_NH3_NO3**2.0*SD_Theta_NH3_NO3**(SD_T(1)-20.)/SD_s*SD_Ae_HS_NH4_Nit/(SD_Ae_HS_NH4_Nit + SD_NH3(1))*SD_O20/(2.0*SD_Ae_HS_O2_Nit + SD_O20)  ! [m/d] = [m2/d2] / [m/d] * [-] * [-] 
			! cb note - SD_Ae_HS_NH4_Nit in above line may need to be temperature corrected
			!Calculate dissolved and particulate (sorbed) fractions 
			SD_fd1 = (1.0D+00/(1.0d+00 + SD_m1*SD_KdNH3))                ! cb note -SD_KdNH3 is not defined anywhere
			SD_fp1 = 1.0D+00 - SD_fd1               != ((m1*KdNH3)/(1 + m1*KdNH3)) 
			SD_fd2 = (1.0D+00/(1.0D+00 + SD_m2*SD_KdNH3))               ! cb note -SD_KdNH3 is not defined anywhere
			SD_fp2 = 1.0D+00 - SD_fd2               != ((m2*KdNH3)/(1 + m2*KdNH3)) 
			! cb note -  since SD_KdNH3 not defined (=0), no nh3 is sorbed to particualtes
			!Write linear system of equations around NH3T 
			!Equation 5.1
			!
			!Layer 1
			SD_a11 = -SD_fd1*SD_KL12 - SD_fp1*SD_w12 - SD_fd1*SD_NH3toNO3 - SD_fd1*SD_s - SD_w2 
            !SD_a11 = -SD_H1 / SD_tc - SD_fd1*SD_KL12 - SD_fp1*SD_w12 - SD_fd1*SD_NH3toNO3 - SD_fd1*SD_s - SD_w2 
			SD_a12 = SD_fd2*SD_KL12 + SD_fp2*SD_w12 
			SD_b1 = -SD_s*SD_NH30                   ![m/d]*[mg/m3] 
			!
			!Layer 2
			SD_a21 = SD_fd1*SD_KL12 + SD_fp1*SD_w12 + SD_w2 
			SD_a22 = -SD_fd2 * SD_KL12 - SD_fp2 * SD_w12 - SD_w2 - SD_H2 / SD_tc
			SD_b2 = -SD_Jn - SD_H2 / SD_tc * SD_NH3Tp(2)
							
			Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_NH3T(1), SD_NH3T(2), NFLog, NFCle) 
			!
			!Dissolved Concentrations 
			SD_NH3(1) = SD_fd1*SD_NH3T(1) 
			SD_NH3(2) = SD_fd2*SD_NH3T(2) 
			
			!CEMA Additional Code
			!Calculate NH4+ and NH3 concentrations
			!
			!SD1_Ammonia = SD_NH3(1)/(1. + 10.**(-SD_pHValue(SegNumI))/10.**(-NH4_NH3_Eqb_Const))
            SD1_Ammonia = SD_NH3(1)/(1. + 10.**(-SD_PH(1))/10.**(-NH4_NH3_Eqb_Const))
			SD1_Ammonium =  SD_NH3(1) - SD1_Ammonia
			!SD2_Ammonia = SD_NH3(2)/(1. + 10.**(-SD_pHValue(SegNumI))/10.**(-NH4_NH3_Eqb_Const))
            SD2_Ammonia = SD_NH3(2)/(1. + 10.**(-SD_PH(2))/10.**(-NH4_NH3_Eqb_Const))
			SD2_Ammonium =  SD_NH3(2) - SD2_Ammonia
			!Calculate dissolved <--> gaseous phase distribution
			SedTemp1 = SD_T(1) + 273.15
            SedTemp2 = SD_T(2) + 273.15
			MW_Constituent = 17 !N = 14, H = 1
            !gp estimated thickness of the aerobic sediment layer 1 (DiToro Appendix B Page 576)  MOVED from below, cb 9/6/13
		   SD_H1 = SD_KL12 * SD_H2 / SD_s
		   If(SD_H1 > SD_H2)Then   
        	  SD_H1 = SD_H2
		   End If
		   SD_AerLayerThick(SegNumI) = SD_H1
			!Layer 1
			VolWater = CellArea(SegNumI)*SD_H1
			Call CEMADisGasPhaseDistribution(SD1_Ammonia, MW_Constituent, SedTemp1, HenryConst_NH3, VolWater, AmmoniaG_SD1, AmmoniaD_SD1)
			!Layer 2
			VolWater = CellArea(SegNumI)*SD_H2
			Call CEMADisGasPhaseDistribution(SD2_Ammonia, MW_Constituent, SedTemp2, HenryConst_NH3, VolWater, AmmoniaG_SD2, AmmoniaD_SD2)
			!End CEMA Additional Code
			
			!	
			!Oxygen Flux due to NH3->NO2, see Equation 23.3 in Chapra (1997) 
			SD_nh3conv = SD_NH3toNO3*SD_NH3(1) 
			SD_NSOD = 2*(32.0/14.0)*SD_NH3toNO3*SD_NH3(1) 
			![gmO/m2-d] = [mol O2/mol N]*[gm O2/mol O2]/[gm N/mol N]* 
			![gm/1000mg] * [m/day] * [mgN/m3] 
			!':::::::::::::::::::::::::::: BEGIN Nitrate:::::::::::::::::::::::::::: 
			!
			!Denitrification in layers 1 and 2 (Equation 4.55) 
			SD_Denit(1) = (SD_Ae_NO3_N2**2*SD_Theta_NO3_N2**(SD_T(1) - 20.)/SD_s) 
			SD_Denit(2) = SD_An_NO3_N2*SD_Theta_NO3_N2**(SD_T(2) - 20.) 
			!
			!Layer 1 
			SD_a11 = -SD_KL12 - SD_Denit(1) - SD_s - SD_w2 
            !SD_a11 = -SD_H1 / SD_tc - SD_KL12 - SD_Denit(1) - SD_s - SD_w2 
			SD_a12 = SD_KL12 
			SD_b1 = -SD_s*SD_NO30 - SD_NH3toNO3*SD_NH3(1) 
			!Layer 2 
			SD_a21 = SD_KL12 + SD_w2 
			SD_a22 = -SD_KL12 - SD_Denit(2) - SD_w2 - SD_H2 / SD_tc
			SD_b2 = -SD_H2 / SD_tc * SD_NO3p(2)
			
			!
			Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_NO3(1), SD_NO3(2), NFLog, NFCle) 
			!
			!Nitrate Flux to water column 
			SD_JNO3 = SD_s*(SD_NO3(1) - SD_NO30) 
			!
			!Denitrification Flux [mgN/m2d] 
			SD_JDenit(1) = SD_Denit(1) * SD_NO3(1) 
			SD_JDenit(2) = SD_Denit(2) * SD_NO3(2) 
			SD_JDenitT = SD_JDenit(1) + SD_JDenit(2) 
			!    
			!Methane consumption due to denitrification (Equation 9.16) 
			!
			!Layer 1 
			!2 for 1/2 N2, 10/8 for eqn balance, 16/12 for O2 --> C and 12/14 for N --> C
			SD_JO2NO3(1) = 2.0D+00 * (16.0D+00 / 12.0D+00) * (10.0D+00 / 8.0D+00) * (12.0D+00 / 14.0D+00) * SD_JDenit(1) 
			![gmO*/m2-d] = [molO/molC]*[(gmO/molO)/(gmC/molC)]*[molC/MolN]* 
			![(gmC/molC)/(gmN/molN)] * [mg/m2d] * [gm/1000mg] 
			!where 2*16/12 is the ubiquitous 32/12 (= 2.67) for oxidation of carbon 
			!Layer 2 
			SD_JO2NO3(2) = (32.0D+00 / 12.0D+00) * (10.0D+00 / 8.0D+00) * (12.0D+00 / 14.0D+00) * SD_JDenit(2) 
			!
			!Sum 
			SD_JO2NO3T = SD_JO2NO3(1) + SD_JO2NO3(2) 
			!
			!Calculate methane flux in oxygen equivalent units, adjusted for 
			!the methane consumed in denitrification                            
			!gp also used if sulfide is produced 
			SD_JC_O2equiv = SD_Jc - SD_JO2NO3T 
			If (SD_JC_O2equiv < 0) Then 
				SD_JC_O2equiv = 0.0D+00
			End If
			!
			!***** Methane/sulfide in O2 equivalents 
			!:::::::::::::::::::::::::::: BEGIN methane/sulfide:::::::::::::::::::::::::::: 
			!
			!gp select for freshwater/saltwater methane/sulfide production from C diagenesis 
			!If (SD_SO4 <= SD_Sulfate_CH4_H2S) Then 
            If (SD_SO4(2) <= SD_Sulfate_CH4_H2S) Then          ! cb 7/26/18
				!gp freshwater methane production, no changes to original code 
				!CSODMAX Equations 10.28 and 10.30 
				!SD_CSODmax = DMin1((2.0D+00 * SD_KL12 * SD_CH4SAT * SD_JC_O2equiv)**2.0D+00, SD_JC_O2equiv)    ![gmO*/m2-d] = sqr([m/d] * [gmO*/m3] * [gmO*/m2-d])   ! SW 10/10/2017 MAJOR ERROR
                SD_CSODmax = DMin1((2.0D+00 * SD_KL12 * SD_CH4SAT * SD_JC_O2equiv)**0.5D+00, SD_JC_O2equiv)    ![gmO*/m2-d] = sqr([m/d] * [gmO*/m3] * [gmO*/m2-d])   ! SW 10/10/2017
				If(SD_CH4CompMethod == 0) Then
					!***********************************************************************
					!Analytical solution for methane
					SD_SECH_ARG = (SD_Ae_CH4_CO2 * SD_Theta_CH4_CO2 ** ((SD_T(1) - 20.) / 2.0)) / SD_s
					!CSOD Equation 10.35
					!The hyperbolic secant is defined as HSec(X) = 2 / (Exp(X) + Exp(-X))
					If (SD_SECH_ARG < 400.0) Then !This is the usual case
						SD_CSOD = SD_CSODmax * (1.0 - (2.0 / (Exp(SD_SECH_ARG) + Exp(-SD_SECH_ARG))))
					Else !HSec(SECH_ARG) < 3.8E-174 ~ 0
						SD_CSOD = SD_CSODmax
					End If
					!***********************************************************************
				Else	
					!
					!NumericalSolution for methane
					SD_CH4toCO2 = (SD_Ae_CH4_CO2 ** 2.0D+00 * SD_Theta_CH4_CO2**((SD_T(1) - 20.0D+00) / 2.0D+00)) / SD_s 
					SD_CH4(1) = (SD_CSODmax + SD_s * SD_CH40) / (SD_CH4toCO2 + SD_s) 
					SD_CSOD = SD_CH4toCO2 * SD_CH4(1) 
				End IF
				
				!0.5CH4 + O2 --> 0.5 CO2 + H2O
				!CO2 Produced = 0.5*(12+32)/32 = 0.6875
				CO2ProducedSrc1L1 = SD_CSOD*0.6875    !g CO2/m²/d
				CO2ProducedCon1L1 = CO2ProducedSrc1L1*SD_tc/SD_H1   !g CO2/m²/d*d/m = g CO2/m³
				
			Else 
				!
				!gp saltwater sulfide production by C diagenesis based on DiToro (2001) Appendix B 
				!***** Calculate dissolved and particulate (sorbed) fractions for sulfide 
				SD_fd1 = (1.0D+00 / (1.0D+00 + SD_m1 * SD_KdH2S1)) 
				SD_fp1 = 1.0D+00 - SD_fd1                                    != ((m1*KdH2S1)/(1 + m1*KdH2S1)) 
				SD_fd2 = (1.0D+00 / (1.0D+00 + SD_m2* SD_KdH2S2)) 
				SD_fp2 = 1.0D+00 - SD_fd2                                    != ((m2*KdH2S2)/(1 + m2*KdH2S2)) 
				!
				!***** Temperature adjusted reaction velocities 
				SD_xappd1 = SD_Ae_H2S_SO4 * SD_Theta_H2S_SO4 **((SD_T(1) - 20.) / 2.0D+00) 
				SD_xappp1 = SD_KappaH2Sp1 * SD_Theta_H2S_SO4 ** ((SD_T(1)- 20.) / 2.0D+00) 
				!
				!***** Transport and Decay terms 
				!Equation B.19
				SD_k1h1d = SD_xappd1**2.0D+00/SD_s*(SD_O20/SD_NormConst_H2S_SO4) + SD_s 
				SD_k1h1p = SD_xappp1**2.0D+00/SD_s*(SD_O20/SD_NormConst_H2S_SO4) 
				SD_k2h2d = 0.0d+00 
				SD_k2h2p = 0.0d+00
				SD_F12 = SD_w12 * SD_fp1 + SD_KL12 * SD_fd1 
				SD_F21 = SD_w12 * SD_fp2 + SD_KL12 * SD_fd2 
				SD_xk1 = SD_k1h1d * SD_fd1 + SD_k1h1p * SD_fp1 
				SD_xk2 = SD_k2h2d * SD_fd2 + SD_k2h2p * SD_fp2 
				!
				!***** Matrix and forcing function 
				SD_a11 = -SD_F12 - SD_xk1 - SD_w2                             !note: -fd1 * s is included in DiToro's -xk1 term 
                !SD_a11 = -SD_H1 / SD_tc - SD_F12 - SD_xk1 - SD_w2                             !note: -fd1 * s is included in DiToro's -xk1 term 
				SD_a21 = SD_F12 + SD_w2 
				SD_a12 = SD_F21 
				SD_b1 = 0.0D+00 
				
				SD_a22 = -SD_F21 - SD_xk2 - SD_w2 - SD_H2 / SD_tc
				SD_b2 = -SD_JC_O2equiv - SD_H2 / SD_tc * SD_HSTp(2)
				
				Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_HST(1), SD_HST(2), NFLog, NFCle) 
				!
				!***** dissolved concentrations 
				SD_HS(1) = SD_fd1*SD_HST(1) 
				SD_HS(2) = SD_fd2*SD_HST(2) 
				SD_CSOD = (SD_xappd1**2.0D+00/SD_s*SD_fd1 + SD_xappp1**2.0D+00/SD_s*SD_fp1)*(SD_O20/SD_NormConst_H2S_SO4)*SD_HST(1) 
				
				!H2S + 2 O2 --> 2 H+ + SO42-
				!SO42- Produced = (32+16*4)/(2*32) = 1.5
				SO4ProducedSrc1L1 = SD_CSOD*1.5    !g SO42-/m²/d
				
				!CH2O + 2 H+ + SO42- --> 2 CO2 + H2S + 2 H2O    !CH2O is represented in O2 equivalent in SD_JC_O2equiv
				!SO42- Consumed = (32+16*4)/(2*16) = 3.0
				!SO4ConsumedSnk1L2 = SD_JC_O2equiv*3.0    !g SO42-/m²/d
                SO4ConsumedSnk1L2 = -SD_JC_O2equiv*3.0    !g SO42-/m²/d   ! cb 7/26/18
				
				!CH2O + 2 H+ + SO42- --> 2 CO2 + H2S + 2 H2O    !CH2O is represented in O2 equivalent in SD_JC_O2equiv
				!CO2 Produced = 2*(12+16*2)/(2*16) = 2.75
				CO2ProducedSrc2L2 = SD_JC_O2equiv*2.75    !g CO2/m²/d
				CO2ProducedCon2L2 = CO2ProducedSrc2L2*SD_tc/SD_H2   !g CO2/m²/d*d/m = g CO2/m³
            End If 
            
            ! Metals Start
          if(includeIron)then
            !calculating dissolved and particulate forms of ferrous iron Fe(II) (Chapra, eqn. 25.89)
            SD_fd1=1.0/(SD_Porosity+KdFe1*SD_Rho*(1.0-SD_Porosity))
            SD_fp1=(KdFe1*SD_Rho*(1.0-SD_Porosity))/(SD_Porosity+KdFe1*SD_Rho*(1.0-SD_Porosity))
            SD_fd2=1.0/(SD_Porosity+KdFe2*SD_Rho*(1.0-SD_Porosity))
            SD_fp2=(KdFe2*SD_Rho*(1.0-SD_Porosity))/(SD_Porosity+KdFe2*SD_Rho*(1.0-SD_Porosity))
            
            SD_Fe2toFeOOH=kfe_oxid*SD_O20*10**(2.0*(sd_ph(1)-7.0))*SD_fd1*SD_Fe2T(1)
            SD_CSOD = SD_CSOD + (0.25*2.0*16.0/55.845)*SD_Fe2toFeOOH  ! lumping in DO consumed by Fe(II)>FeOOH into CSOD
            SD_FeOOHtoFe2=kfe_red*SD_FeOOH(1)
            !FeOOH + 0.25CH2O + 2H+ > Fe(II) + 0.25CO2 + 1.75H20
            ! CO2 Produced = 0.25 * (12 + 2*16) / 55.845 = 0.197         !g CO2/m²/d
            CO2ProducedSrc1L2 = CO2ProducedSrc1L2 + SD_FeOOHtoFe2 * 0.197
            
            !Write linear system of equations around total ferrous iron SD_Fe2T
			!Equation 5.1
			!
			!Layer 1
			SD_a11 = -SD_fd1*SD_KL12 - SD_fp1*SD_w12 - SD_fd1*SD_s - SD_w2 - SD_Fe2toFeOOH*SD_H1
			SD_a12 = SD_fd2*SD_KL12 + SD_fp2*SD_w12 
			SD_b1 = -SD_s*SD_Fe20
			!
			!Layer 2
			SD_a21 = SD_fd1*SD_KL12 + SD_fp1*SD_w12 + SD_w2 
			SD_a22 = -SD_fd2 * SD_KL12 - SD_fp2 * SD_w12 - SD_w2 - SD_H2 / SD_tc
			SD_b2 = - SD_H2 / SD_tc * SD_Fe2T(2) - SD_H2 * SD_FeOOHtoFe2
							
			Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_Fe2T(1), SD_Fe2T(2), NFLog, NFCle) 
            
            ! dissolved Fe(2)
            SD_Fe2(1)=SD_fd1*SD_Fe2T(1)
            SD_Fe2(2)=SD_fd2*SD_Fe2T(2)
		
            SD_JFe2= SD_s * (SD_Fe2(1)-SD_Fe20)        
            
            !Write linear system of equations around total ferrous iron SD_FeOOH
			!Equation 5.1
			!
			!Layer 1
			SD_a11 = -SD_w12 -  SD_w2 
			SD_a12 = SD_w12 
			SD_b1 = -SD_JFeOOHin- SD_Fe2toFeOOH*SD_H1
			!
			!Layer 2
			SD_a21 = SD_w12 + SD_w2 
			SD_a22 = - SD_w12 - SD_w2 - SD_H2*SD_FeOOHtoFe2  - SD_H2 / SD_tc
			SD_b2 = - SD_H2 / SD_tc * SD_FeOOH(2)
							
			Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_FeOOH(1), SD_FeOOH(2), NFLog, NFCle)                         		
          end if
            
          if(includeManganese)then
            !calculating dissolved and particulate forms of Mn(II) (Chapra, eqn. 25.89)
            SD_fd1=1.0/(SD_Porosity+KdMn1*SD_Rho*(1.0-SD_Porosity))
            SD_fp1=(KdMn1*SD_Rho*(1.0-SD_Porosity))/(SD_Porosity+KdMn1*SD_Rho*(1.0-SD_Porosity))
            SD_fd2=1.0/(SD_Porosity+KdMn2*SD_Rho*(1.0-SD_Porosity))
            SD_fp2=(KdMn2*SD_Rho*(1.0-SD_Porosity))/(SD_Porosity+KdMn2*SD_Rho*(1.0-SD_Porosity))
            
            SD_Mn2toMnO2=kMn_oxid*SD_O20*10**(2.0*(sd_ph(1)-7.0))*SD_fd1*SD_Mn2T(1)
            SD_CSOD = SD_CSOD + (16.0/54.94)*SD_Mn2toMnO2  ! lumping in DO consumed by Mn(II)>MnO2 into CSOD
            SD_MnO2toMn2=kMn_red*SD_MnO2(1)
            !MnO2 + 0.5CH2O + 2H+ > Mn(II) + 0.5CO2 + 1.5H20
            ! CO2 Produced = 0.5 * (12 + 2*16) / 54.94 = 0.400         !g CO2/m²/d
            CO2ProducedSrc1L2 = CO2ProducedSrc1L2 + SD_MnO2toMn2 * 0.400
            
            !Write linear system of equations around total Mn(II) SD_Mn2T
			!Equation 5.1
			!
			!Layer 1
			SD_a11 = -SD_fd1*SD_KL12 - SD_fp1*SD_w12 - SD_fd1*SD_s - SD_w2 - SD_Mn2toMnO2*SD_H1
			SD_a12 = SD_fd2*SD_KL12 + SD_fp2*SD_w12 
			SD_b1 = -SD_s*SD_Mn20
			!
			!Layer 2
			SD_a21 = SD_fd1*SD_KL12 + SD_fp1*SD_w12 + SD_w2 
			SD_a22 = -SD_fd2 * SD_KL12 - SD_fp2 * SD_w12 - SD_w2 - SD_H2 / SD_tc
			SD_b2 = - SD_H2 / SD_tc * SD_Mn2T(2) - SD_H2 * SD_MnO2toMn2
							
			Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_Mn2T(1), SD_Mn2T(2), NFLog, NFCle) 
            
            ! dissolved Mn(II)
            SD_Mn2(1)=SD_fd1*SD_Mn2T(1)
            SD_Mn2(2)=SD_fd2*SD_Mn2T(2)
		
            SD_JMn2= SD_s * (SD_Mn2(1)-SD_Mn20)        
            
            !Write linear system of equations around manganese dioxide SD_MnO2
			!Equation 5.1
			!
			!Layer 1
			SD_a11 = -SD_w12 -  SD_w2 
			SD_a12 = SD_w12 
			SD_b1 = -SD_JMnO2in- SD_Mn2toMnO2*SD_H1
			!
			!Layer 2
			SD_a21 = SD_w12 + SD_w2 
			SD_a22 = - SD_w12 - SD_w2 - SD_H2*SD_MnO2toMn2  - SD_H2 / SD_tc
			SD_b2 = - SD_H2 / SD_tc * SD_MnO2(2)
							
			Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_MnO2(1), SD_MnO2(2), NFLog, NFCle)                         		
          end if
            ! Metals End
            
			SD_SOD = (SD_SODold + SD_CSOD + SD_NSOD)/2 
			iter = iter + 1 
			SD_ea = Abs((SD_SOD - SD_SODold)/SD_SOD)*100.0D+00 
			!SOD = (SODold + CSOD + NSOD) / 2 
			!iter = iter + 1 
			!ea = Abs((SOD - SODold) / SOD) * 100 
			If (SD_ea <= SD_es) Then 
				Exit  
			ElseIf (iter >= maxit) Then 
				!MsgBox "SOD iterations exceeded" 
				!Write(NFlog,*) 'SOD iterations exceeded'
				Write(*,*) 'SOD iterations exceeded'
				Exit
			End If
			SD_SODold = SD_SOD 
		End Do
		SD_s = SD_SOD/SD_O20
		!
		!Flux to water column
		SD_JSOD = SD_SOD 
		SD_JNH4 = SD_s * (SD_NH3(1) - SD_NH30) 
		SD_JNO3 = SD_s * (SD_NO3(1) - SD_NO30) 
        
        !CEMA Additional Code
        !If (SD_SO4 <= SD_Sulfate_CH4_H2S) Then
        If (SD_SO4(2) <= SD_Sulfate_CH4_H2S) Then         ! 7/26/18
            
        Else
            !Calculate H2S and HS- concentrations
			!
			!SD1_SulfiMinus = SD_HS(1)/(1 + 10**(-SD_pHValue(SegNumI))/10**(-HS_H2S_Eqb_Const))
            SD1_SulfiMinus = SD_HS(1)/(1 + 10**(-SD_PH(1))/10**(-HS_H2S_Eqb_Const))
			SD1_Sulfide =  SD_HS(1) - SD1_SulfiMinus
			!SD2_SulfiMinus = SD_HS(2)/(1 + 10**(-SD_pHValue(SegNumI))/10**(-HS_H2S_Eqb_Const))
            SD2_SulfiMinus = SD_HS(2)/(1 + 10**(-SD_PH(2))/10**(-HS_H2S_Eqb_Const))
			SD2_Sulfide =  SD_HS(2) - SD2_SulfiMinus
			
        End If
        !End CEMA Additional Code
        
		!
		!gp   methane or sulfide fluxes produced from C diagenesis 
		!If (SD_SO4 < SD_Sulfate_CH4_H2S) Then
        If (SD_SO4(2) <= SD_Sulfate_CH4_H2S) Then                 ! cb 7/26/18
			!gp freshwater sediment fluxes - methane
			If (SD_CH4CompMethod == 0) Then
				!***********************************************************************
				!Aqueous methane flux to water column
				SD_SJCH4 = SD_CSODmax - SD_CSOD
				SD_JCH4 = SD_SJCH4
				!Gaseous methane flux to water column
				SD_JCH4g = SD_JC_O2equiv - SD_JCH4 - SD_CSOD           ! flux in gO2/m2/day
				!***********************************************************************
				SD_CH4(1) = SD_JCH4/SD_s + SD_CH40
			Else
				!***********************************************************************
				!numerical solution for methane
				SD_JCH4 = SD_s * (SD_CH4(1) - SD_CH40)
				SD_JCH4g = SD_JC_O2equiv - SD_JCH4 - SD_CSOD           ! flux in gO2/m2/day
				!***********************************************************************
			End If
			SD_JHS = 0     !gp
		Else !gp
			!gp marine sediment fluxes - sulfide
			SD_JCH4 = 0
			SD_JCH4g = 0
			SD_JHS = SD_s * SD_HS(1)
		End If  !gp
		
		!CEMA Additional Code
		!SD_JSO4 = SD_s * SD_SO4                
		!End CEMA Additional Code
		!
		!gp estimated thickness of the aerobic sediment layer 1 (DiToro Appendix B Page 576)  MOVED ABOVE
		!SD_H1 = SD_KL12 * SD_H2 / SD_s
		!If(SD_H1 > SD_H2)Then
        !		    SD_H1 = SD_H2
		!End If
		!SD_AerLayerThick(SegNumI) = SD_H1
		
		!Bubbles formation
		!Aerobic layer is thin and so ignore gas formation from aerobic layer
		  !H2S
		!Calculate dissolved <--> gaseous phase distribution
		SedTemp1 = SD_T(1) + 273.15
		MW_Constituent = 36. !S = 34, H = 1
		!Calculate dissolved <--> gaseous phase distribution
		!VolWater = CellArea(SegNumI)*SD_H2
		SedTemp2 = SD_T(2) + 273.15
		!Layer 1
    	VolWater = CellArea(SegNumI)*SD_H1
		Call CEMADisGasPhaseDistribution(SD1_Sulfide, MW_Constituent, SedTemp1, HenryConst_H2S, VolWater, SulfideG_SD1, SulfideD_SD1)
		!Layer 2
		VolWater = CellArea(SegNumI)*SD_H2
		Call CEMADisGasPhaseDistribution(SD2_Sulfide, MW_Constituent, SedTemp2, HenryConst_H2S, VolWater, SulfideG_SD2, SulfideD_SD2)
		!Old TConc(1,SegNumI) = SD2_Sulfide
		!Old SConc(1,SegNumI) = (TConc(1,SegNumI) - TConcP(1,SegNumI))/dlt    !gm/m³/s
		TConc(1,SegNumI) = SulfideG_SD2*BubbAccFraction + TConcP(1,SegNumI)
		SConc(1,SegNumI) = (SulfideG_SD2*BubbAccFraction)/dlt    !gm/m³/s
		!Old TConcP(1,SegNumI) = SD2_Sulfide
		TConcP(1,SegNumI) = TConc(1,SegNumI)
		
		!CH4
		TConc(2,SegNumI) = SD_CH4(1)*BubbAccFraction + TConcP(2,SegNumI)
		SConc(2,SegNumI) = (SD_CH4(1)*BubbAccFraction)/dlt    !gm/m³/s
		!Old TConcP(1,SegNumI) = SD2_Sulfide
		TConcP(2,SegNumI) = TConc(2,SegNumI)
		
		!NH3
		!Old TConc(3,SegNumI) = SD2_Ammonia
		!Old SConc(3,SegNumI) = (TConc(3,SegNumI) - TConcP(3,SegNumI))/dlt    !gm/m³/s
		TConc(3,SegNumI) = SD_NH3(2) !SD2_Ammonia*BubbAccFraction + TConcP(3,SegNumI)
		SConc(3,SegNumI) = (SD_NH3(2)-TConcP(3,SegNumI))/dlt    !gm/m³/s
		!Old TConcP(3,SegNumI) = SD2_Ammonia
		TConcP(3,SegNumI) = TConc(3,SegNumI)
		
		!CO2
		!Old TConc(4,SegNumI) = CO2ProducedCon2L2
		!Old SConc(4,SegNumI) = (TConc(4,SegNumI) - TConcP(4,SegNumI))/dlt    !gm/m³/s
		TConc(4,SegNumI) = CO2ProducedCon2L2*BubbAccFraction + TConcP(4,SegNumI)
		SConc(4,SegNumI) = (CO2ProducedCon2L2*BubbAccFraction)/dlt    !gm/m³/s
		!Old TConcP(4,SegNumI) = CO2ProducedCon2L2
		TConcP(4,SegNumI) = TConc(4,SegNumI)
		
		!H2S
		DissolvedGasSediments(1, SegNumI) = SulfideD_SD2
	    !CH4
		DissolvedGasSediments(2, SegNumI) = SD_CH4(1)
		!NH4
		DissolvedGasSediments(3, SegNumI) = AmmoniaD_SD2
		!CO2
		DissolvedGasSediments(4, SegNumI) = 0.d0
		
		IF(BUBBLES_CALCULATION)Call GasBubblesFormation(BubbleRadiusSed(SegNumI), dlt, VolWater)    
		!Update Sediment Bed Concentrations - additional mass balance (due to PW release)
		!Update concentrations of SO4, NH3, NO3, H2S
		!Sulfate        
		!SD_SO4 = SD_SO4 + (SO4ProducedSrc1L1 - SO4ConsumedSnk1L2)*SD_tc/(SD_H1 + SD_H2)
        !if(SD_SO4 < 0.0)SD_SO4=0.0                
		
         !Write linear system of equations around sulfate   cb 7/26/18
			!Equation 5.1
			!
			!Layer 1
			SD_a11 = -SD_KL12 -SD_s            
			SD_a12 = SD_KL12
          ! SO4ProducedSrc1L1 converted from  g SO4/m²/d to g S/m²/d ; 32/(4*16+32)=0.33333
            !H2S + 2 O2 --> 2 H+ + SO42-            
			SD_b1 = -SD_s*SD_so40 - SO4ProducedSrc1L1  * 0.3333333
			!
			!Layer 2
			SD_a21 = SD_KL12
			SD_a22 = -SD_KL12 -  SD_H2 / SD_tc
          ! SO4ConsumedSnk1L2 converted from  g SO4/m²/d to g S/m²/d ; 32/(4*16+32)=0.33333
            !CH2O + 2 H+ + SO42- --> 2 CO2 + H2S + 2 H2O
			SD_b2 = - SD_H2 / SD_tc * SD_so4(2) - SO4ConsumedSnk1L2 * 0.33333333
							
			Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_so4(1), SD_so4(2), NFLog, NFCle)
            
                        
            ! calculating diffusive flux between layer 1 and water column
            SD_Jso4= SD_s * (SD_so4(1)-SD_so40)
		
        !Write linear system of equations around Total inorganic carbon
			!Equation 5.1
			!
			!Layer 1
			SD_a11 = -SD_KL12 -SD_s            
			SD_a12 = SD_KL12
          ! COCO2ProducedSrc1L1 converted from  g CO2/m²/d to g C/m²/d ; 12/(2*16+12)=0.272
			SD_b1 = -SD_s*SD_TIC0 - CO2ProducedSrc1L1  * 0.272
			!
			!Layer 2
			SD_a21 = SD_KL12
			SD_a22 = -SD_KL12 -  SD_H2 / SD_tc
          ! COCO2ProducedSrc1L2 converted from  g CO2/m²/d to g C/m²/d ; 12/(2*16+12)=0.272
			SD_b2 = - SD_H2 / SD_tc * SD_tic(2) - CO2ProducedSrc2L2 * 0.272
							
			Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_tic(1), SD_tic(2), NFLog, NFCle)
            
            ! calculating diffusive flux between layer 1 and water column
            SD_JTIC= SD_s * (SD_tic(1)-SD_TIC0)
        
            if(IncludeAlkalinity)then
        !Write linear system of equations around Total Alkalinity
        ! Nitrification of ammonium results in an alkalinity decrease: 2 eq. alk per 1 mole ammonium
        ! Denitrification of nitrate (to nitrogen gas) results in an alkalinity increase: 1 eq. alk per 1 mole nitrate
			!Equation 5.1
			!
			!Layer 1
			SD_a11 = -SD_KL12 -SD_s
            !SD_a11 = -SD_H1/SD_tc -SD_KL12 -SD_s
			SD_a12 = SD_KL12
			SD_b1 = -SD_s*SD_ALK0 +2.0*SD_NH3toNO3*SD_NH3(1) - SD_NO3(1)*SD_Denit(1)                   ![m/d]*[mg/m3] 
			!
			!Layer 2
			SD_a21 = SD_KL12
			SD_a22 = -SD_KL12 -  SD_H2 / SD_tc
			SD_b2 = - SD_H2 / SD_tc * SD_alk(2) - SD_NO3(2)*SD_Denit(2)
							
			Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_alk(1), SD_alk(2), NFLog, NFCle)
            
            ! calculating diffusive flux between layer 1 and water column
            SD_Jalk= SD_s * (SD_alk(1)-SD_alk0)
        
            end if
            
            !calculating dissolved and particulate forms of phosphorus (Chapra, eqn. 25.89)
            SD_fd1=1.0/(SD_Porosity+Kdp1*SD_Rho*(1.0-SD_Porosity))
            SD_fp1=(Kdp1*SD_Rho*(1.0-SD_Porosity))/(SD_Porosity+Kdp1*SD_Rho*(1.0-SD_Porosity))
            SD_fd2=1.0/(SD_Porosity+Kdp2*SD_Rho*(1.0-SD_Porosity))
            SD_fp2=(Kdp2*SD_Rho*(1.0-SD_Porosity))/(SD_Porosity+Kdp2*SD_Rho*(1.0-SD_Porosity))
            
            !Write linear system of equations around total phosphate SD_PO4T
			!Equation 5.1
			!
			!Layer 1
			SD_a11 = -SD_fd1*SD_KL12 - SD_fp1*SD_w12 - SD_fd1*SD_s - SD_w2 
            !SD_a11 = -SD_H1/SD_tc - SD_fd1*SD_KL12 - SD_fp1*SD_w12 - SD_fd1*SD_s - SD_w2 
			SD_a12 = SD_fd2*SD_KL12 + SD_fp2*SD_w12 
			SD_b1 = -SD_s*SD_PO40                   ![m/d]*[mg/m3] 
			!
			!Layer 2
			SD_a21 = SD_fd1*SD_KL12 + SD_fp1*SD_w12 + SD_w2 
			SD_a22 = -SD_fd2 * SD_KL12 - SD_fp2 * SD_w12 - SD_w2 - SD_H2 / SD_tc
			SD_b2 = -SD_Jp - SD_H2 / SD_tc * SD_PO4T(2)
							
			Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_PO4T(1), SD_PO4T(2), NFLog, NFCle) 
            
            ! dissolved PO4
            SD_PO4(1)=SD_fd1*SD_PO4T(1)
            SD_PO4(2)=SD_fd2*SD_PO4T(2)
		
            SD_Jp= SD_s * (SD_PO4(1)-SD_PO40)
        
            
        !Write linear system of equations for temperature
            SD_rhowcp=4.186*1.0e6   ! units J g-1 C-1 * g m-3= J C-1 m-3
			!Equation 5.1
			!
			!Layer 1
			SD_a11 = -SD_KL12 -SD_s            
			SD_a12 = SD_KL12          
			SD_b1 = -SD_s*SD_Tw 
			!
			!Layer 2
			SD_a21 = SD_KL12 * SD_rhowcp
			SD_a22 = -SD_KL12 * SD_rhowcp -  SD_H2 / SD_tc * SD_rhowcp   + SD_Ksw
			SD_b2 = - SD_H2 * SD_rhowcp / SD_tc * SD_T(2) - SD_Ksw * SD_Tsed
							
			Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_T(1), SD_T(2), NFLog, NFCle)
            
            ! calculating diffusive heat flux between layer 1 and water column
            SD_JT= SD_s * SD_rhowcp* (SD_T(1)-SD_Tw)
            
        if(IncludeDynamicpH)then
        ! Sediment pH for layers 1 and 2
        SD_TDS=0.0  !sediments not simulating tds at the moment        
        SD_POCT1=0.0 ! POC not predicted for aerobic layer.
        call PH_SEDIMENTS(SD_T(1),sd_tic(1),sd_alk(1),sd_nh3(1),sd_po4(1),SD_POCT1,SD_TDS(1),SD_PH(1))  !layer 1
        
        call PH_SEDIMENTS(SD_T(2),sd_tic(2),sd_alk(2),sd_nh3(2),sd_po4(2),SD_POCT2,SD_TDS(2),SD_PH(2))  !layer 2
        end if
            
		!Transfer to variables that will be used to hook up the model to W2 water quality model
        
        ! SW compute CO2 in layer 1 and layer 2 10/20/2017
            !HION = 10.0**SD_PH(1)
            !SD_CO2(1) = SD_TIC(1)/(1.0+K1/HION+K1*K2/(HION*HION))    ! K1 AND K2 FROM PH_SEDIMENTS SUBROUTINE
            !
            !
            !HION = 10.0**SD_PH(2)
            !SD_CO2(2) = SD_TIC(2)/(1.0+K1/HION+K1*K2/(HION*HION))
            !
            
		
		!1
		!Bubbles release
		!Bubbles of CH4g, NH3g and H2Sg
		
		
		!2
		!Porewater release
		!Flux of CH4d, NH3d + NH4d, H2Sd + HSd, SO42-d, NO3d, CO2d 
		!Volume of porewater
		VolWater = CellArea(SegNumI)*SD_H2
        
        !Aerobic Layer
        PW_RelRate1 = PorewaterRelRate(SegNumI)*SD_H1/(SD_H1 + SD_H2)*86400.0  !m³/d
        !Anaerobic Layer
        PW_RelRate2 = PorewaterRelRate(SegNumI)*SD_H2/(SD_H1 + SD_H2)*86400.0  !m³/d
		
		!Dissolved_CH4_Src = SD_JCH4*SD_tc/(SD_H1 + SD_H2)   !CHECK CODE This is incorrect original code since SD_tc is in units of s*d/s or d ==> g/m3 Why didvide by SD_H1 and SD_H2 inncorrect
        !Dissolved_CH4_Src = SD_JCH4*SD_tc/(SD_H1 + SD_H2)   !SW 10/10/2017    gO2/m2/d ==> gO2/m3/d
		
		!Dissolved_NH3_Src = (AmmoniaD_SD1*SD_H1 + AmmoniaD_SD2*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
		!Dissolved_NH3_Src = Dissolved_NH3_Src + (SD1_Ammonium*SD_H1 + SD2_Ammonium*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))  !gm/m³/s
        Dissolved_NH3_Src = (AmmoniaD_SD1*PW_RelRate1 + AmmoniaD_SD2*PW_RelRate2)/(CellThickness*CellArea(SegNumI))                      !gm/m³/d
		Dissolved_NH3_Src = Dissolved_NH3_Src + (SD1_Ammonium*PW_RelRate1 + SD2_Ammonium*PW_RelRate2)/(CellThickness*CellArea(SegNumI))  !gm/m³/d
		
        !Dissolved_H2S_Src = (SulfideD_SD1*SD_H1 + SulfideD_SD2*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
		!Dissolved_H2S_Src = Dissolved_H2S_Src + (SD1_SulfiMinus*SD_H1 + SD2_SulfiMinus*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))  !gm/m³/s
        Dissolved_H2S_Src = (SulfideD_SD1*PW_RelRate1 + SulfideD_SD2*PW_RelRate2)/(CellThickness*CellArea(SegNumI))                      !gm/m³/d
		Dissolved_H2S_Src = Dissolved_H2S_Src + (SD1_SulfiMinus*PW_RelRate1 + SD2_SulfiMinus*PW_RelRate2)/(CellThickness*CellArea(SegNumI))  !gm/m³/d
		
		!Dissolved_NO3_Src = (SD_NO3(1)*SD_H1 + SD_NO3(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
        Dissolved_NO3_Src = (SD_NO3(1)*PW_RelRate1 + SD_NO3(2)*PW_RelRate2)/(CellThickness*CellArea(SegNumI))                      !gm/m³/d
		
		!Dissolved_SO4_Src = SD_SO4 * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
        !Dissolved_SO4_Src = SD_SO4 * (PW_RelRate1+PW_RelRate2)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
        Dissolved_SO4_Src = (SD_SO4(1)*PW_RelRate1 + SD_SO4(2)*PW_RelRate2)/(CellThickness*CellArea(SegNumI))     ! cb 7/26/18
        		
		!Dissolved_CO2_Src = (CO2ProducedCon1L1*SD_H1 + CO2ProducedCon2L2*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
        Dissolved_CO2_Src = (CO2ProducedCon1L1*PW_RelRate1 + CO2ProducedCon2L2*PW_RelRate2)/(CellThickness*CellArea(SegNumI))                      !gm/m³/d
        !Dissolved_ALK_Src = (SD_ALK(1)*SD_H1 + SD_ALK(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
        if(IncludeAlkalinity)then
        Dissolved_ALK_Src = (SD_ALK(1)*PW_RelRate1 + SD_ALK(2)*PW_RelRate2)/(CellThickness*CellArea(SegNumI))                      !gm/m³/d
        else
          Dissolved_ALK_Src = 0.0
        end if
        !Dissolved_PO4_Src = (SD_PO4(1)*SD_H1 + SD_PO4(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
        Dissolved_PO4_Src = (SD_PO4(1)*PW_RelRate1 + SD_PO4(2)*PW_RelRate2)/(CellThickness*CellArea(SegNumI))                      !gm/m³/d
        !Dissolved_Fe2_Src = (SD_Fe2(1)*SD_H1 + SD_Fe2(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
        if(IncludeIron)then
        Dissolved_Fe2_Src = (SD_Fe2(1)*PW_RelRate1 + SD_Fe2(2)*PW_RelRate2)/(CellThickness*CellArea(SegNumI))                      !gm/m³/d
        else
          Dissolved_Fe2_Src = 0.0
        end if
        !Dissolved_Mn2_Src = (SD_Mn2(1)*SD_H1 + SD_Mn2(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m³/s
        if(IncludeManganese)then
        Dissolved_Mn2_Src = (SD_Mn2(1)*PW_RelRate1 + SD_Mn2(2)*PW_RelRate2)/(CellThickness*CellArea(SegNumI))                      !gm/m³/d
        else
          Dissolved_Mn2_Src = 0.0
        end if
        !Sediment_Heat_Src = SD_rhowcp*(SD_T(1)*SD_H1 + SD_T(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !J/m³/s
        Sediment_Heat_Src = SD_rhowcp*(SD_T(1)*PW_RelRate1 + SD_T(2)*PW_RelRate2)/(CellThickness*CellArea(SegNumI))                      !J/m³/d
		
        !3
        !Diffusive flux of NH3, NO3, CH4, SO4, H2S, CO2
        Dissolved_NH3_Src = Dissolved_NH3_Src + SD_JNH4/CellThickness   !g/m³/d SD_JNH4/CellThickness = g/m²/d/m = g/m³/d
        Dissolved_NO3_Src = Dissolved_NO3_Src + SD_JNO3/CellThickness   !g/m³/d SD_JNO3/CellThickness = g/m²/d/m = g/m³/d
        !Dissolved_CH4_Src = Dissolved_CH4_Src + SD_JCH4/CellThickness   !g/m³/d SD_JCH4/CellThickness = g/m²/d/m = g/m³/d commented out because double counting... 
        Dissolved_CH4_Src = Dissolved_CH4_Src + SD_JCH4/CellThickness   !g/m³/d SD_JCH4/CellThickness = g/m²/d/m = g/m³/d commented out because double counting... ! SW 10/10/2017 added back because eliminated the problem above
        Dissolved_SO4_Src = Dissolved_SO4_Src + SD_JSO4/CellThickness   !g/m³/d SD_JSO4/CellThickness = g/m²/d/m = g/m³/d	
        Dissolved_H2S_Src = Dissolved_H2S_Src + SD_JHS/CellThickness
        Dissolved_CO2_Src = Dissolved_CO2_Src + SD_JTIC/CellThickness
        Dissolved_Alk_Src = Dissolved_ALK_Src + SD_JALK/CellThickness
        Dissolved_PO4_Src = Dissolved_PO4_Src + SD_Jp/CellThickness
        Dissolved_Fe2_Src = Dissolved_Fe2_Src + SD_JFe2/CellThickness
        Dissolved_Mn2_Src = Dissolved_Mn2_Src + SD_JMn2/CellThickness
        
        GasReleaseCH4=GasReleaseCH4+SD_JCH4g*B(KB(SegNumI),SegNumI)*DLX(SegNumI)*SD_tc/2.67                    ! Flux of gas CH4 in gC   SW 10/10/2017
        !GasReleaseCO2=GasReleaseCO2+
        
        Sediment_Heat_Src = Sediment_Heat_Src + SD_JT/CellThickness  ! Heat J/m3/d  SD_JT/CellThickness = J/m3/d/m= J/m3/d

		!Flux to CSOD and NSOD
		Dissolved_O2_Snk = SD_JSOD/CellThickness   !g/m³/d SD_JSOD/CellThickness = g/m²/d/m = g/m³/d
        !BOD sediment feature
		if(IncludeCEMAGenBODConstituents)Call ComputeGenBODFate
				
		!Convert all source/sink to g/m³/s from g/m³/d
		!NH3, NO3, CH4, SO4, DO, CO2, H2S
		Dissolved_NH3_Src   = Dissolved_NH3_Src/86400.0     !g/m³/d --> g/m³/s
        Dissolved_NO3_Src   = Dissolved_NO3_Src/86400.0     !g/m³/d --> g/m³/s
        Dissolved_CH4_Src   = Dissolved_CH4_Src/86400.0     !g/m³/d --> g/m³/s
        Dissolved_SO4_Src   = Dissolved_SO4_Src/86400.0     !g/m³/d --> g/m³/s
		Dissolved_CO2_Src   = Dissolved_CO2_Src/86400.0     !g/m³/d --> g/m³/s
		Dissolved_H2S_Src   = Dissolved_H2S_Src/86400.0     !g/m³/d --> g/m³/s
		Dissolved_O2_Snk    = Dissolved_O2_Snk/86400.0      !g/m³/d --> g/m³/s
        Dissolved_ALK_Src   = Dissolved_ALK_Src/86400.0      !g/m³/d --> g/m³/s
        Dissolved_PO4_Src   = Dissolved_PO4_Src/86400.0      !g/m³/d --> g/m³/s
        Dissolved_Fe2_Src   = Dissolved_Fe2_Src/86400.0      !g/m³/d --> g/m³/s
        Dissolved_Mn2_Src   = Dissolved_Mn2_Src/86400.0      !g/m³/d --> g/m³/s
        if(IncludeCEMAGenBODConstituents)then
		  Do GenBODNum = 1, NumGenBODConstituents
		    SedGenBODSS(GenBODNum) = SedGenBODSS(GenBODNum)/CellThickness   !g/m³/d SD_JSO4/CellThickness = g/m²/d/m = g/m³/d
		    SedGenBODSS(GenBODNum) = SedGenBODSS(GenBODNum)/86400.0         !g/m³/d --> g/m³/s
          End Do !GenBODNum                   
        end if
       Sediment_Heat_Src   = Sediment_Heat_Src/86400.0      !J/m³/d --> J/m³/s
       
       ! resuspension of POM, POP, and PON
       Do iTemp = 1, 2  ! only labile and refractory, not including inert for now                             
         LPOM_Resuspension = LPOM_Resuspension + SD_EPOC(itemp)/CellThickness/ORGC(JW)
         RPOM_Resuspension = RPOM_Resuspension + SD_EPOC(itemp)/CellThickness/ORGC(JW)
         LPOMN_Resuspension = LPOMN_Resuspension + SD_EPON(itemp)/CellThickness
         RPOMN_Resuspension = RPOMN_Resuspension + SD_EPON(itemp)/CellThickness
         LPOMP_Resuspension = LPOMP_Resuspension + SD_EPOP(itemp)/CellThickness
         RPOMP_Resuspension = RPOMP_Resuspension + SD_EPOP(itemp)/CellThickness
      end do
       
		
		!Update Sediment Bed Concentrations - additional mass balance for SO4
		!NH3, NO3, H2S don't need to be updated as the diffusive flux is included in the mass balance equation
		!Sulfate
		!SD_SO4 = SD_SO4 - SD_JSO4*SD_tc/(SD_H1 + SD_H2) !Only SO4 has this term as NH3, NO3 and H2S equations include diffusion
        !if(SD_SO4 < 0.0)SD_SO4=0.0        
		!PW release loss only reduces total mass but not the concentration. Lost mass is lost along with water so the concentration remains the same
		
		Do iTemp = 1, 3
			SD_POC22(iTemp) = SD_POC2(iTemp)
			SD_PON22(iTemp) = SD_PON2(iTemp)
			SD_POP22(iTemp) = SD_POP2(iTemp)
		End Do

		Do iTemp = 1, 2
			SD_NH3Tp2(iTemp) = SD_NH3T(iTemp)
			SD_NO3p2(iTemp)  = SD_NO3(iTemp)
			!SD_PO4Tp2(iTemp) = SD_PO4T(iTemp)
			SD_NH3p2(iTemp)  = SD_NH3(iTemp)
			SD_PO4p2(iTemp)  = SD_PO4(iTemp)
			SD_CH4p2(iTemp)  = SD_CH4(iTemp)     !gp
			SD_HSTp2(iTemp)  = SD_HST(iTemp)     !gp
			SD_HSp2(iTemp)   = SD_HS(iTemp)       !gp
		End Do
		SD_BEN_STRp2 = SD_BEN_STR
		
	RETURN
	
	Entry ComputeGenBODFate
	
    
    ! BOD debug
    !if(BODtestout)then
    !  open(8765,file='Sed_bod_debug.opt',status='unknown')
    !  write(8765,'("        JDAY   BOD_S2_L1   BOD_S3_L1   BOD_S4_L1   BOD_S2_L2   BOD_S3_L2   BOD_S4_L2")')
    !  BODtestout=.false.
    !end if
   
    
	SD_BODD = 0.d0
    Dissolved_BOD_Src = 0.d0
    SedGenBODSS = 0.d00
    
    !Decay and oxygen consumption
    If(SD_O20 > SD_Ox_Threshol)Then !Only if aerobic conditions exist
        Do GenBODNum = 1, NumGenBODConstituents
	        Do iTemp = 1, 1 !Only Aerobic layer consumption and decay
    	         
	             SedGenBODDecayRate(GenBODNum,SegNumI) = SedGenBODConsRate(GenBODNum,SegNumI) * SedGenBODConsTCoeff(GenBODNum,SegNumI)**(SD_Tw - 20.)
	             SedGenBODConc(GenBODNum,SegNumI,iTemp) = (1-SedGenBODDecayRate(GenBODNum,SegNumI)*SD_tc)*SedGenBODConc(GenBODNum,SegNumI,iTemp)
	             Dissolved_O2_Snk = Dissolved_O2_Snk + SedGenBODDecayRate(GenBODNum,SegNumI)*SedGenBODConc(GenBODNum,SegNumI,iTemp)
                        	         
	        End Do !iTemp
        End Do !GenBODNum                   
    End If
    
    !Limit the concentrations to minimum of zero
    Do GenBODNum = 1, NumGenBODConstituents
        Do iTemp = 1, 1 !Only Aerobic layer consumption and decay
            If(SedGenBODConc(GenBODNum,SegNumI,iTemp) < 0.d00)SedGenBODConc(GenBODNum,SegNumI,iTemp) = 0.d00   
        End Do !iTemp
    End Do !GenBODNum
    
    !Diffusion
    Do GenBODNum = 1, NumGenBODConstituents
        
        ! Diffusion between layers
        SD_BOD = SD_KL12 * (SedGenBODConc(GenBODNum,SegNumI,2) - SedGenBODConc(GenBODNum,SegNumI,1))    
        SD_BODD(GenBODNum,1) =  1.d0*SD_BOD !If SD_BOD > 0 then Layer 1 receives    gm/m²/d
        SD_BODD(GenBODNum,2) = -1.d0*SD_BOD !If SD_BOD < 0 then Layer 2 receives    gm/m²/d
        
        !Diffusion to overlying water column
        SD_BOD = 0.d0
        !Only if there are aerobic conditions
        !If(SD_O20 > 1.0e-10)SD_BOD = SD_s * (SedGenBODConc(GenBODNum,SegNumI,1) - SedGenBODConc0(GenBODNum))
        If(SD_O20 > SD_Ox_Threshol)SD_BOD = SD_s * (SedGenBODConc(GenBODNum,SegNumI,1) - SedGenBODConc0(GenBODNum))
        SD_BODD(GenBODNum,1) =  SD_BODD(GenBODNum,1) - 1.d0*SD_BOD !If SD_BOD < 0 then Layer 1 receives    gm/m²/d
        SedGenBODSS(GenBODNum) = 1.d0*SD_BOD
        
    End Do !GenBODNum
    
    !Porewater release
    Do GenBODNum = 1, NumGenBODConstituents
    
        !Aerobic Layer
        PW_RelRate = PorewaterRelRate(SegNumI)*SD_H1/(SD_H1 + SD_H2)*86400.  !m³/d
        Dissolved_BOD_Src(GenBODNum,1) = -1.d0*SedGenBODConc(GenBODNum,SegNumI,1) * PW_RelRate/CellArea(SegNumI)   !gm/m²/d    
        
        !Anaerobic Layer
        PW_RelRate = PorewaterRelRate(SegNumI)*SD_H2/(SD_H1 + SD_H2)*86400.  !m³/d
        Dissolved_BOD_Src(GenBODNum,2) = -1.d0*SedGenBODConc(GenBODNum,SegNumI,2) * PW_RelRate/CellArea(SegNumI)   !gm/m²/d    
        
        SedGenBODSS(GenBODNum) = SedGenBODSS(GenBODNum) - Dissolved_BOD_Src(GenBODNum,1) - Dissolved_BOD_Src(GenBODNum,2)
    
    End Do !GenBODNum
 
    !Resulting BOD concentrations in each layer
    Do GenBODNum = 1, NumGenBODConstituents
        !Aerobic Layer
        SedGenBODConc(GenBODNum,SegNumI,1) = SedGenBODConc(GenBODNum,SegNumI,1) + ( SD_BODD(GenBODNum,1) +  Dissolved_BOD_Src(GenBODNum,1) )/SD_H1 *SD_tc
        !Anaerobic Layer
        SedGenBODConc(GenBODNum,SegNumI,2) = SedGenBODConc(GenBODNum,SegNumI,2) + ( SD_BODD(GenBODNum,2) +  Dissolved_BOD_Src(GenBODNum,2) )/SD_H2 *SD_tc
    End Do !GenBODNum
    
    !Limit the concentrations to minimum of zero
    Do GenBODNum = 1, NumGenBODConstituents
        Do iTemp = 1, 1 !Only Aerobic layer consumption and decay
            If(SedGenBODConc(GenBODNum,SegNumI,iTemp) < 0.d00)SedGenBODConc(GenBODNum,SegNumI,iTemp) = 0.d00   
        End Do !iTemp
    End Do !GenBODNum
    
    !if(segnumi == 4)then
    !  write(8765,'(f12.6,100f10.6)')jday,(SedGenBODConc(1,ii,1), ii=2,4),(SedGenBODConc(1,ii,2), ii=2,4)
    !end if
        	
	RETURN
	
	
	Entry EndSedimentFluxVariables
		    
        MFTSedFlxVars(SegNumI,1)		    =   0.d00
        MFTSedFlxVars(SegNumI,2)		    =   0.d00
        MFTSedFlxVars(SegNumI,3)		    =   0.d00
        MFTSedFlxVars(SegNumI,4)		    =   0.d00
        MFTSedFlxVars(SegNumI,5)		    =   0.d00
        MFTSedFlxVars(SegNumI,6)		    =   0.d00
        MFTSedFlxVars(SegNumI,7)		    =   0.d00
        MFTSedFlxVars(SegNumI,8)		    =   0.d00
        MFTSedFlxVars(SegNumI,9)		    =   0.d00
        MFTSedFlxVars(SegNumI,10)	        =   0.d00
        MFTSedFlxVars(SegNumI,11)	        =   0.d00
        MFTSedFlxVars(SegNumI,12)	        =   0.d00
        MFTSedFlxVars(SegNumI,13)	        =   0.d00
        MFTSedFlxVars(SegNumI,14)	        =   0.d00
        MFTSedFlxVars(SegNumI,15)	        =   0.d00

        MFTSedFlxVars(SegNumI,16)           =   0.d00
        MFTSedFlxVars(SegNumI,17)           =   0.d00
        MFTSedFlxVars(SegNumI,18)           =   0.d00
        MFTSedFlxVars(SegNumI,19)           =   0.d00
        MFTSedFlxVars(SegNumI,20)           =   0.d00
        MFTSedFlxVars(SegNumI,21)           =   0.d00
        MFTSedFlxVars(SegNumI,22)           =   0.d00
        MFTSedFlxVars(SegNumI,23)           =   0.d00
        MFTSedFlxVars(SegNumI,24)           =   0.d00
        MFTSedFlxVars(SegNumI,25)           =   0.d00
    	
        MFTSedFlxVars(SegNumI,26)           =   0.d00
        MFTSedFlxVars(SegNumI,27)           =   0.d00
        MFTSedFlxVars(SegNumI,28)           =   0.d00
    	
        MFTSedFlxVars(SegNumI,29)           =   0.d00
        MFTSedFlxVars(SegNumI,30)           =   0.d00
        
        MFTSedFlxVars(SegNumI,31)           =   0.d00
        MFTSedFlxVars(SegNumI,32)           =   0.d00
        MFTSedFlxVars(SegNumI,33)           =   0.d00
        MFTSedFlxVars(SegNumI,34)           =   0.d00
        MFTSedFlxVars(SegNumI,35)           =   0.d00
        MFTSedFlxVars(SegNumI,36)           =   0.d00
        MFTSedFlxVars(SegNumI,37)           =   0.d00
        MFTSedFlxVars(SegNumI,38)           =   0.d00
        MFTSedFlxVars(SegNumI,39)           =   0.d00
        MFTSedFlxVars(SegNumI,40)           =   0.d00
        MFTSedFlxVars(SegNumI,41)           =   0.d00
        MFTSedFlxVars(SegNumI,42)           =   0.d00
        MFTSedFlxVars(SegNumI,43)           =   0.d00
        MFTSedFlxVars(SegNumI,44)           =   0.d00
        MFTSedFlxVars(SegNumI,45)           =   0.d00
        MFTSedFlxVars(SegNumI,46)           =   0.d00
        MFTSedFlxVars(SegNumI,47)           =   0.d00
                	        
		        
	Return
	
	!Entry ComputeCEMADiagenesisSourceSinks  SW 6/27/2017
 !       
 !       Do JW=1, NWB
 !           KT = KTWB(JW)
 !           Do JB=BS(JW),BE(JW)
 !               IU = CUS(JB)
 !               ID = DS(JB)
 !               Do SegNumI = IU, ID
 !               
 !                                       
 !               End Do !SegNumI
 !           End Do !JB
 !       End Do !JW
 !  Return	
        
    ENTRY PH_SEDIMENTS(t1sed,ticsed,alksed,nh4sed,po4sed,pocsed,tdssed,phsed) ! Enhancements added for buffering by ammonia, phosphate, and OM ! SR 01/01/12
! pH and carbonate species
 
     T1K = t1sed+273.15
     CART = ticsed/12011. ! SR 01/01/12
     ALKT = alksed/50044. ! SR 01/01/12
     AMMT = NH4sed/14006.74 ! SR 01/01/12
     PHOST = PO4sed/30973.762 ! SR 01/01/12
     !OMCT = (LDOM(K,I)+RDOM(K,I))*ORGC(JW)/12011. ! moles carbon per liter from DOM ! SR 01/01/12
     omct=0.0  ! DOM is not simulated in the sediments yet...
     !IF (POM_BUFFERING) OMCT = OMCT + (LPOM(K,I)+RPOM(K,I))*ORGC(JW)/12011. ! SR 01/01/12
     IF (POM_BUFFERING) OMCT = OMCT + pocsed/12011. ! SR 01/01/12
     omct=0.0    ! 
 !**** Ionic strength
     IF (FRESH_WATER(JW)) S2 = 2.5E-05*TDSsed
     IF (SALT_WATER(JW)) S2 = 1.47E-3+1.9885E-2*TDSsed+3.8E-5*TDSsed*TDSsed
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
     KP3 = 10.0**(-12.38) *HPO4T/PO4T ! Dean (1985 )! SR 01/01/12
!**** pH evaluation
     !PHT = -PH(K,I)-2.1
     !IF (PH(K,I) <= 0.0) PHT = -14.0
     PHT = -phsed-2.1
     IF (phsed <= 0.0) PHT = -14.0
     INCR1 = 10.0
     DO N=1,3
        F1 = 1.0
        INCR1 = INCR1/10.0
        ITER1 = 0
        DO WHILE (F1 > 0.0 .AND. ITER1 < 12)
          PHT = PHT+INCR1
          HION = 10.0**PHT
          F1 = CART*K1*(HION+2.0*K2)/(HION*HION+K1*HION+K1*K2)+KW/HION-ALKT-HION/HT ! SR 01/01/12
          IF (AMMONIA_BUFFERING) THEN ! SR 01/01/12
            F1 = F1 + AMMT*KAMM/(HION+KAMM) ! SR 01/01/12
          END IF ! SR 01/01/12
          IF (PHOSPHATE_BUFFERING) THEN ! SR 01/01/12
            F1 = F1+ PHOST*( KP1*KP2*HION + 2*KP1*KP2*KP3 - HION*HION*HION ) &
                /( HION*HION*HION + KP1*HION*HION + KP1*KP2*HION + KP1*KP2*KP3) ! SR 01/01/12
          END IF ! SR 01/01/12
          IF (OM_BUFFERING) THEN ! SR 01/01/12
            DO JA=1,NAG ! SR 01/01/12
              F1 = F1 + OMCT*SDEN(JA)*( 1.0/(1.0+HION*(10.0**PK(JA))) - 1.0/(1.0+(10.0**(PK(JA)-4.5))) ) ! SR 01/01/12
            END DO ! SR 01/01/12
          END IF ! SR 01/01/12
          ITER1 = ITER1+1
        END DO
        PHT = PHT-INCR1
     END DO
!**** pH, carbon dioxide, bicarbonate, and carbonate concentrations
     HION = 10.0**PHT
     PHsed = -PHT
     CO2sed = TICsed/(1.0+K1/HION+K1*K2/(HION*HION))
     HCO3sed = TICsed/(1.0+HION/K1+K2/HION)
     CO3sed = TICsed/((HION*HION)/(K1*K2)+HION/K2+1.0)
 
RETURN


Entry CEMAWindInducedSedimentResuspension

  epsilon=0.0
  FETCHW = FETCHD(SegNumI,JB)
  IF (COS(PHI(JW)-PHI0(SegNumI)) < 0.0) FETCHW = FETCHU(SegNumI,JB)
  FETCHW = MAX(FETCHW,BI(KT,SegNumI),DLX(SegNumI))
  U2    = WIND(JW)*WSC(SegNumI)*WIND(JW)*WSC(SegNumI)+NONZERO
  COEF1 = 0.53  *(G*DEPTHB(kb(SegNumI),SegNumI)/U2)**0.75
  COEF2 = 0.0125*(G*FETCHW/U2)**0.42
  COEF3 = 0.833* (G*DEPTHB(kb(SegNumI),SegNumI)/U2)**0.375
  COEF4 = 0.077* (G*FETCHW/U2)**0.25
  HS    = 0.283 *U2/G*0.283*TANH(COEF1)*TANH(COEF2/TANH(COEF1))
  !TS    = 2.0*PI*U2/G*1.2*  TANH(COEF3)*TANH(COEF4/TANH(COEF3))
  TS    = 2.0*PI*sqrt(U2)/G*1.2*  TANH(COEF3)*TANH(COEF4/TANH(COEF3))   ! cb 5/9/14
  LW0    = G*TS*TS/(2.0*PI)  
 
  LW1 = LW0
  LW  = LW0*TANH(2.0*PI*DEPTHB(Kb(SegNumI),SegNumI)/LW1)
  DO WHILE (ABS(LW-LW1) > 0.001)
    LW1 = LW
    LW  = LW0*TANH(2.0*PI*DEPTHB(Kb(SegNumI),SegNumI)/LW1)
  END DO
  COEF = MIN(710.0,2.0*PI*DEPTHB(Kb(SegNumI),SegNumI)/LW)
  UORB = PI*HS/TS*100.0/SINH(COEF)
  TAU  = 0.003*UORB*UORB
  IF (TAU-TAUCRPOM > 0.0) EPSILON = MAX(0.0,0.008/49.0*(TAU-TAUCRPOM)**3*10000.0/DLT)						        
  SD_E = EPSILON*DLX(SegNumI)*BI(KB(SegNumI),SegNumI)/VOL(KB(SegNumI),SegNumI)  ! SD_E: g/m^2/s

Return
        
Entry CEMABottomScourResuspension

  if(cao_method)then
    reyn_resusp = dia_POM * sqrt(spgrav_POM * g * dia_POM)
    if(reyn_resusp < 6.0)then
      crshields = 0.1414 * reyn_resusp ** (-0.2306)
    else if(reyn_resusp >= 6.0 .and. reyn_resusp <= 282.8)then
      crshields = ( 1.0 + (0.0223 * reyn_resusp)**2.8358)**0.3542 / (3.0946 * reyn_resusp**0.6769)
    else if(reyn_resusp > 282.8)then
      crshields = 0.045
    end if
  end if  
  molvisc_h2o = 1.79e-6 * exp(0.0266 * SD_T(1))
  shields = SD_taubot / (g*(spgrav_POM-1.0) * dia_POM)
  
  Vscour = 0.00033*(shields/crshields - 1.0)*(spgrav_POM-1.0)**0.6 * g**0.6 * dia_POM**0.8/molvisc_h2o
  c_bottom = (c2(kb(SegNumI),SegNumI,NLPOM) + c2(kb(SegNumI),SegNumI,NRPOM)) * dexp( poms(jw) * h(kb(SegNumI),jw) / DZ(kb(SegNumI)-1, SegNumI))
    
  if(spgrav_POM < 1.2)then
    c_bottom2=1.0
  else if(spgrav_POM >= 1.2 .and. spgrav_POM < 1.8)then
     c_bottom2= 1.0 * (1.8 - spgrav_POM)/0.6 + 3.0 * (spgrav_POM - 1.2)/0.6
  else if(spgrav_POM >= 1.8 .and. spgrav_POM <= 2.2)then
      c_bottom2= 3.0 * (2.2 - spgrav_POM)/0.4 + 5.0 * (spgrav_POM - 1.8)/0.4
  else if(spgrav_POM > 2.2)then
    c_bottom2=5.0
  end if 
  
  c_bottom=dmin1(c_bottom,c_bottom2)

  if(Vscour > 0.0)then
  SD_E = c_bottom*Vscour
  else
    SD_E = 0.0
  end if


Return

End Subroutine
    
!
! 
Subroutine Lin_Sys(a11, a12, a21, a22, b1, b2, x1, x2,NFLog,NFCle) 
	Real(8) a11, a12, a21, a22, b1, b2, x1, x2
	Byte NFLog, NFCle
	!from 03-Nov-2003 version of Q2KMaster 
	!This subroutine solves a linear system of 2 equations and 2 unknowns 
	If (a11 * a22 - a12 * a21 == 0) Then 
		!MsgBox "The sediment flux solution matrix is singular: " & a11 & ", " & a12 & ", " & a21 & ", " & a22 
		Write(NFLog,'(a)') 'The sediment flux solution matrix is singular: '
		Write(NFLog,*) 'a11  == ', a11, 'a12 = ', a12, 'a21 = ', a21, 'a22 = ',a22 
		Write(NFCle,'(a)') 'The sediment flux solution matrix is singular: '
		Write(NFCle,*) 'a11  == ', a11, 'a12 = ', a12, 'a21 = ', a21, 'a22 = ',a22 
		Write(NFCle,*) 'Error in the solution of linear system of equations used in sediment diagenesis model'
		Write(NFLog,*) 'Error in the solution of linear system of equations used in sediment diagenesis model'
        write(w2err,*)'Error in Sediment Diagenesis - see log files: Please review the sediment diagensis parameters'
		Stop 'Please review the sediment diagensis parameters' 
	End If 
	x1 = (a22 * b1 - a12 * b2) / (a11 * a22 - a12 * a21) 
	x2 = (a11 * b2 - a21 * b1) / (a11 * a22 - a12 * a21) 

End Subroutine 


Subroutine CEMADisGasPhaseDistribution(CTotal, MolWt, GasTemp, HenryConst, Vwtr, CGasPh, CLiqPh)
    Use CEMAVars
    
    Real(8) CTotal, MolWt, GasTemp, HenryConst
    Real(8) Vwtr, CGasPh, CLiqPh
    
    CGasPh  =   CTotal/(1 + GasConst_R*GasTemp/HenryConst)
    CLiqPh  =   CTotal*(GasConst_R*GasTemp/HenryConst)/(1 + GasConst_R*GasTemp/HenryConst)

End Subroutine

