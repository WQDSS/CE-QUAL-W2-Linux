Subroutine CEMA_W2_Input
    Use MAIN 
    Use GLOBAL
    Use GEOMC
    Use CEMAVars
    
    ! Type declarations
    IMPLICIT NONE
    
    Logical SkipLoop, file_exists
    Character(256) MessageTemp
    
    ! start CEMA testing variables
    integer monzz,dayzz,yearzz,ninp

     SD_global=.false.
     IncludeIron=.false.
     IncludeManganese=.false.
     IncludeDynamicpH=.false.
     IncludeAlkalinity=.false.    
     ngh2s=0
     ngch4=0
     ngso4=0
     ngturb=0
     ngFe2=0
     ngFeOOH=0
     ngMn2=0
     ngMnO2=0
     ngmft=0
    
    ! end CEMA testing variables
 
 	INQUIRE(FILE="W2_CEMA_Input.npt", EXIST=file_exists)   ! file_exists will be TRUE if the file
    
	If(.NOT. file_exists) Then
	    CEMARelatedCode = .FALSE.
	    IncludeBedConsolidation = .FALSE.
	    Return
	End If
	CEMARelatedCode = .TRUE. 
	
    CEMAFilN=NUNIT ; NUNIT=NUNIT+1   ! SW 5/26/15
	Open(CEMAFilN, File = "W2_CEMA_Input.npt")
	Open(CEMALogFilN, File = "CEMALogFile.opt")
	
	!Read Header
	SkipLoop = .FALSE.
	Do While(.NOT. SkipLoop)
		Read(CEMAFilN,'(a)')MessageTemp
		If(index(MessageTemp, "$") == 0)SkipLoop = .TRUE.
	End Do

	BackSpace(CEMAFilN)
    Read(CEMAFilN,*)MessageTemp, SD_global
    If(.NOT. SD_global) Then         
	    CEMARelatedCode = .FALSE.
	    IncludeBedConsolidation = .FALSE.
	    Return
	End If
	Read(CEMAFilN,*)MessageTemp, nh2s     ! cb 2/18/13  reading constituent #'s for h2s, ch4, so4, co2, turbidity, and mft
    Read(CEMAFilN,*)MessageTemp, nch4
    Read(CEMAFilN,*)MessageTemp, nso4    
    Read(CEMAFilN,*)MessageTemp, nturb
    Read(CEMAFilN,*)MessageTemp, nFe2
    Read(CEMAFilN,*)MessageTemp, nFeOOH
    Read(CEMAFilN,*)MessageTemp, nMn2
    Read(CEMAFilN,*)MessageTemp, nMnO2
    Read(CEMAFilN,*)MessageTemp, nmft
    ngh2s=nh2s-NGCS+1
    ngch4=nch4-NGCS+1
    ngso4=nso4-NGCS+1
    ngturb=nturb-NGCS+1
    ngFe2=nFe2-NGCS+1
    ngFeOOH=nFeOOH-NGCS+1
    ngMn2=nMn2-NGCS+1
    ngMnO2=nMnO2-NGCS+1
    ngmft=nmft-NGCS+1
    
    Read(CEMAFilN,*)MessageTemp, IncludeBedConsolidation
	Read(CEMAFilN,*)MessageTemp, LayerAddThkFrac
	Read(CEMAFilN,*)MessageTemp, NumConsolidRegns
	
	Allocate(ConsolidationType(NumConsolidRegns),ConstConsolidRate(NumConsolidRegns))
	Allocate(ConstPoreWtrRate(NumConsolidRegns),ConsolidRateRegnFil(NumConsolidRegns))
	Allocate(ConsRegSegSt(NumConsolidRegns), ConsRegSegEn(NumConsolidRegns))
	
	Read(CEMAFilN,*)MessageTemp, (ConsRegSegSt(i), i = 1, NumConsolidRegns)
	Read(CEMAFilN,*)MessageTemp, (ConsRegSegEn(i), i = 1, NumConsolidRegns)
	Read(CEMAFilN,*)MessageTemp, (ConsolidationType(i), i = 1, NumConsolidRegns)
	Read(CEMAFilN,*)MessageTemp, (ConstConsolidRate(i), i = 1, NumConsolidRegns)
	Read(CEMAFilN,*)MessageTemp, (ConsolidRateRegnFil(i), i = 1, NumConsolidRegns)
	
	Read(CEMAFilN,*)MessageTemp, BedElevationInit
	Read(CEMAFilN,*)MessageTemp, BedPorosityInit
	Read(CEMAFilN,*)MessageTemp, CEMAPWpH
    Read(CEMAFilN,*)MessageTemp, IncludeDynamicpH
    Read(CEMAFilN,*)MessageTemp, IncludeAlkalinity
	Read(CEMAFilN,*)MessageTemp, CEMAParticleSize
	CEMAParticleSize = 1.d-6*CEMAParticleSize   !Microns to m
	Read(CEMAFilN,*)MessageTemp, CEMASedimentType
	Read(CEMAFilN,*)MessageTemp, CEMASedimentDensity
	Read(CEMAFilN,*)MessageTemp, CEMASedimentSVelocity
	CEMASedimentSVelocity = CEMASedimentSVelocity/86400.d0  !m/d to m/s
	Read(CEMAFilN,*)MessageTemp, CEMASedimentProcessesInc
	
    Allocate(BedElevation(IMX), BedElevationLayer(IMX), BedPorosity(IMX))
    !allocate (sedcellwidth(imx))
    Allocate(ConsolidRegnNum(IMX), BedConsolidRate(IMX), PorewaterRelRate(IMX))
    Allocate(CEMASedConc(IMX,KMX))
    Allocate(CEMACumPWRelease(IMX), CEMALayerAdded(IMX), CEMASSApplied(IMX))
    Allocate(CEMACumPWToRelease(IMX),CEMACumPWReleased(IMX))
    Allocate(NumCEMAPWInst(IMX))
    Allocate(ApplyCEMAPWRelease(IMX))
    Allocate(CEMACumPWReleaseRate(IMX))
    !Allocate(EndBedConsolidation(IMX))
    Allocate(EndBedConsolidation(IMX),BedConsolidationSeg(IMX))   ! cb 6/28/18
    Allocate(CEMATSSCopy(KMX,IMX))
    Allocate(VOLCEMA(NBR))
    
    BedElevation = BedElevationInit
    BedElevationLayer = 0.d00
    !sedcellwidth=0.d00
    BedPorosity = BedPorosityInit
    BedConsolidRate = 0.d00
    PorewaterRelRate = 0.d00
    CEMASedConc = 0.d00
    CEMACumPWRelease = 0.d00
    CEMACumPWReleaseRate = 0.d00
    CEMACumPWToRelease = 0.d00
    CEMACumPWReleased = 0.d00
    EndBedConsolidation = .FALSE.
    BedConsolidationSeg = .FALSE.  ! cb 6/28/18
    VOLCEMA = 0.d00
    NumCEMAPWInst = 0
    ApplyCEMAPWRelease = .FALSE.
    
    Read(CEMAFilN,*)MessageTemp, WriteBESnp
    Read(CEMAFilN,*)MessageTemp, WritePWSnp
    
    Read(CEMAFilN,*)MessageTemp, IncludeFFTLayer
    Read(CEMAFilN,*)MessageTemp, NumFFTActivePrds
    If(IncludeFFTLayer)FirstTimeInFFTCode = .TRUE.
    Allocate(FFTActPrdSt(NumFFTActivePrds), FFTActPrdEn(NumFFTActivePrds))
    Allocate(FFTLayConc(IMX))
    Read(CEMAFilN,*)MessageTemp, (FFTActPrdSt(i), i = 1, NumFFTActivePrds) 
    Read(CEMAFilN,*)MessageTemp, (FFTActPrdEn(i), i = 1, NumFFTActivePrds)
    Read(CEMAFilN,*)MessageTemp, InitFFTLayerConc
    Read(CEMAFilN,*)MessageTemp, FFTLayerSettVel
    FFTLayerSettVel = FFTLayerSettVel/86400.d00 !m/d --> m/s
    FFTLayConc = 0.d00
    FFTActPrd = 1
    MoveFFTLayerDown = .FALSE.
    Read(CEMAFilN,*)MessageTemp, MoveFFTLayerDown
    
    Read(CEMAFilN,*)MessageTemp, IncludeCEMASedDiagenesis
    If(IncludeCEMASedDiagenesis)sediment_diagenesis=.true.
    !If(.NOT.IncludeCEMASedDiagenesis)Return  ! commented out cb 6/29/18
    If(IncludeCEMASedDiagenesis)FirstTimeinCEMAMFTSedDiag = .TRUE.
    IncludeCEMAGenBODConstituents = .FALSE.
    Read(CEMAFilN,*)MessageTemp, NumRegnsSedimentBedComposition
    
    Allocate(SDRegnPOC_T(NumRegnsSedimentBedComposition), SDRegnPON_T(NumRegnsSedimentBedComposition), SDRegnSul_T(NumRegnsSedimentBedComposition))
    Allocate(SDRegnPOP_T(NumRegnsSedimentBedComposition))
    Allocate(SDRegnH2S_T(NumRegnsSedimentBedComposition), SDRegnNH3_T(NumRegnsSedimentBedComposition), SDRegnCH4_T(NumRegnsSedimentBedComposition))
    Allocate(SDRegnTIC_T(NumRegnsSedimentBedComposition),SDRegnALK_T(NumRegnsSedimentBedComposition), SDRegnPO4_T(NumRegnsSedimentBedComposition))
    Allocate(SDRegnFe2_T(NumRegnsSedimentBedComposition),SDRegnFeOOH_T(NumRegnsSedimentBedComposition))
    Allocate(SDRegnMn2_T(NumRegnsSedimentBedComposition),SDRegnMnO2_T(NumRegnsSedimentBedComposition))
    Allocate(SDRegnT_T(NumRegnsSedimentBedComposition))
    Allocate(SedBedInitRegSegSt(NumRegnsSedimentBedComposition), SedBedInitRegSegEn(NumRegnsSedimentBedComposition))
	
	Read(CEMAFilN,*)MessageTemp, (SedBedInitRegSegSt(i), i = 1, NumRegnsSedimentBedComposition)
	Read(CEMAFilN,*)MessageTemp, (SedBedInitRegSegEn(i), i = 1, NumRegnsSedimentBedComposition)
    
    Read(CEMAFilN,*)MessageTemp, (SDRegnT_T(i), i = 1, NumRegnsSedimentBedComposition)
    Read(CEMAFilN,*)MessageTemp, (SDRegnPOC_T(i), i = 1, NumRegnsSedimentBedComposition)
    Read(CEMAFilN,*)MessageTemp, (SDRegnPON_T(i), i = 1, NumRegnsSedimentBedComposition)
    Read(CEMAFilN,*)MessageTemp, (SDRegnPOP_T(i), i = 1, NumRegnsSedimentBedComposition)
    Read(CEMAFilN,*)MessageTemp, (SDRegnSul_T(i), i = 1, NumRegnsSedimentBedComposition)
    Read(CEMAFilN,*)MessageTemp, (SDRegnNH3_T(i), i = 1, NumRegnsSedimentBedComposition)
    Read(CEMAFilN,*)MessageTemp, (SDRegnPO4_T(i), i = 1, NumRegnsSedimentBedComposition)
    Read(CEMAFilN,*)MessageTemp, (SDRegnH2S_T(i), i = 1, NumRegnsSedimentBedComposition)
    Read(CEMAFilN,*)MessageTemp, (SDRegnCH4_T(i), i = 1, NumRegnsSedimentBedComposition)
    Read(CEMAFilN,*)MessageTemp, (SDRegnTIC_T(i), i = 1, NumRegnsSedimentBedComposition)
    Read(CEMAFilN,*)MessageTemp, (SDRegnALK_T(i), i = 1, NumRegnsSedimentBedComposition)
    Read(CEMAFilN,*)MessageTemp, (SDRegnFe2_T(i), i = 1, NumRegnsSedimentBedComposition)
    Read(CEMAFilN,*)MessageTemp, (SDRegnFeOOH_T(i), i = 1, NumRegnsSedimentBedComposition)
    Read(CEMAFilN,*)MessageTemp, (SDRegnMn2_T(i), i = 1, NumRegnsSedimentBedComposition)
    Read(CEMAFilN,*)MessageTemp, (SDRegnMnO2_T(i), i = 1, NumRegnsSedimentBedComposition)
    
    Read(CEMAFilN,*)MessageTemp, IncludeCEMAGenBODConstituents    
    !If(IncludeCEMAGenBODConstituents)Then
        Read(CEMAFilN,*)MessageTemp, NumGenBODConstituents
        Allocate(SedGenBODName(NumGenBODConstituents))
        Read(CEMAFilN,*)MessageTemp, (SedGenBODName(i), i = 1, NumGenBODConstituents)

        Read(CEMAFilN,*)MessageTemp, NumGenBODInitRegns
        Allocate(SedGenBODInit(NumGenBODConstituents, NumGenBODInitRegns))
        Allocate(SedGenBODRegSegSt(NumGenBODInitRegns), SedGenBODRegSegEn(NumGenBODInitRegns))
        Read(CEMAFilN,*)MessageTemp, (SedGenBODRegSegSt(i), i = 1, NumGenBODInitRegns)
	    Read(CEMAFilN,*)MessageTemp, (SedGenBODRegSegEn(i), i = 1, NumGenBODInitRegns)
    	
	    Do i = 1, NumGenBODConstituents	       
            Read(CEMAFilN,*)MessageTemp, (SedGenBODInit(i,j), j = 1, NumGenBODInitRegns)
	    End Do
        
        Read(CEMAFilN,*)MessageTemp, NumGenBODConsumptionRegions
        Allocate(SedGenBODConsRegSegSt(NumGenBODConsumptionRegions), SedGenBODConsRegSegEn(NumGenBODConsumptionRegions))
        Allocate(SedGenBODRegnRate(NumGenBODConstituents, NumGenBODConsumptionRegions), SedGenBODRegnTCoeff(NumGenBODConstituents, NumGenBODConsumptionRegions))
        Read(CEMAFilN,*)MessageTemp, (SedGenBODConsRegSegSt(i), i = 1, NumGenBODConsumptionRegions)
	    Read(CEMAFilN,*)MessageTemp, (SedGenBODConsRegSegEn(i), i = 1, NumGenBODConsumptionRegions)
        
        Do i = 1, NumGenBODConstituents
	        Read(CEMAFilN,*)MessageTemp, (SedGenBODRegnRate(i,j), j = 1, NumGenBODConsumptionRegions)
	    End Do
	    Do i = 1, NumGenBODConstituents
	        Read(CEMAFilN,*)MessageTemp, (SedGenBODRegnTCoeff(i,j), j = 1, NumGenBODConsumptionRegions)
	    End Do
    !End If
	
    Read(CEMAFilN,*)MessageTemp, NumRegnsSedimentDiagenesis
    
    Allocate(SDRegnPOC_L_Fr(NumRegnsSedimentDiagenesis), SDRegnPOC_R_Fr(NumRegnsSedimentDiagenesis), SDRegnPON_L_Fr(NumRegnsSedimentDiagenesis))
    Allocate(SDRegnPON_R_Fr(NumRegnsSedimentDiagenesis), SDRegnPW_DiffCoeff(NumRegnsSedimentDiagenesis), SDRegnOx_Threshold(NumRegnsSedimentDiagenesis))
    Allocate(SDRegnPOP_L_Fr(NumRegnsSedimentDiagenesis), SDRegnPOP_R_Fr(NumRegnsSedimentDiagenesis))
    Allocate(SDRegnAe_NH3_NO3_L(NumRegnsSedimentDiagenesis), SDRegnAe_NH3_NO3_H(NumRegnsSedimentDiagenesis), SDRegnAe_NO3_N2_L(NumRegnsSedimentDiagenesis))
    Allocate(SDRegnAe_NO3_N2_H(NumRegnsSedimentDiagenesis), SDRegnAn_NO3_N2(NumRegnsSedimentDiagenesis), SDRegnAe_CH4_CO2(NumRegnsSedimentDiagenesis))
    Allocate(SDRegnAe_HS_NH4_Nit(NumRegnsSedimentDiagenesis), SDRegnAe_HS_O2_Nit(NumRegnsSedimentDiagenesis), SDRegn_Theta_PW(NumRegnsSedimentDiagenesis))
    Allocate(SDRegn_Theta_NH3_NO3(NumRegnsSedimentDiagenesis), SDRegn_Theta_NO3_N2(NumRegnsSedimentDiagenesis), SDRegn_Theta_CH4_CO2(NumRegnsSedimentDiagenesis))
    Allocate(SDRegn_Sulfate_CH4_H2S(NumRegnsSedimentDiagenesis), SDRegnAe_H2S_SO4(NumRegnsSedimentDiagenesis), SDRegn_Theta_H2S_SO4(NumRegnsSedimentDiagenesis))
    Allocate(SDRegn_NormConst_H2S_SO4(NumRegnsSedimentDiagenesis), SDRegn_MinRate_PON_Lab(NumRegnsSedimentDiagenesis), SDRegn_MinRate_PON_Ref(NumRegnsSedimentDiagenesis))
    Allocate(SDRegn_MinRate_PON_Ine(NumRegnsSedimentDiagenesis), SDRegn_MinRate_POC_Lab(NumRegnsSedimentDiagenesis), SDRegn_MinRate_POC_Ref(NumRegnsSedimentDiagenesis))
    Allocate(SDRegn_MinRate_POC_Ine(NumRegnsSedimentDiagenesis), SDRegn_Theta_PON_Lab(NumRegnsSedimentDiagenesis), SDRegn_Theta_PON_Ref(NumRegnsSedimentDiagenesis))
    Allocate(SDRegn_Theta_PON_Ine(NumRegnsSedimentDiagenesis), SDRegn_Theta_POC_Lab(NumRegnsSedimentDiagenesis), SDRegn_Theta_POC_Ref(NumRegnsSedimentDiagenesis))
    Allocate(SDRegn_Theta_POC_Ine(NumRegnsSedimentDiagenesis), SDRegn_CH4CompMethod(NumRegnsSedimentDiagenesis), SDRegn_POMResuspMethod(NumRegnsSedimentDiagenesis))
    Allocate(SDRegn_Theta_POP_Lab(NumRegnsSedimentDiagenesis), SDRegn_Theta_POP_Ref(NumRegnsSedimentDiagenesis), SDRegn_Theta_POP_Ine(NumRegnsSedimentDiagenesis))
    Allocate(SDRegn_MinRate_POP_Lab(NumRegnsSedimentDiagenesis), SDRegn_MinRate_POP_Ref(NumRegnsSedimentDiagenesis), SDRegn_MinRate_POP_Ine(NumRegnsSedimentDiagenesis))
    Allocate(SedBedDiaRCRegSegSt(NumRegnsSedimentDiagenesis), SedBedDiaRCRegSegEn(NumRegnsSedimentDiagenesis))
	
	Read(CEMAFilN,*)MessageTemp, (SedBedDiaRCRegSegSt(i), i = 1, NumRegnsSedimentDiagenesis)
	Read(CEMAFilN,*)MessageTemp, (SedBedDiaRCRegSegEn(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegnPOC_L_Fr(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegnPOC_R_Fr(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegnPON_L_Fr(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegnPON_R_Fr(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegnPOP_L_Fr(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegnPOP_R_Fr(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegnPW_DiffCoeff(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegnOx_Threshold(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegnAe_NH3_NO3_L(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegnAe_NH3_NO3_H(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegnAe_NO3_N2_L(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegnAe_NO3_N2_H(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegnAn_NO3_N2(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegnAe_CH4_CO2(i), i = 1, NumRegnsSedimentDiagenesis)   !Eq. 10.35
    Read(CEMAFilN,*)MessageTemp, (SDRegnAe_HS_NH4_Nit(i), i = 1, NumRegnsSedimentDiagenesis)   !Eq. 3.3
    Read(CEMAFilN,*)MessageTemp, (SDRegnAe_HS_O2_Nit(i), i = 1, NumRegnsSedimentDiagenesis)   !Eq. 3.3
    Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_PW(i), i = 1, NumRegnsSedimentDiagenesis)   
    Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_NH3_NO3(i), i = 1, NumRegnsSedimentDiagenesis)   
    Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_NO3_N2(i), i = 1, NumRegnsSedimentDiagenesis)  
    Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_CH4_CO2(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_Sulfate_CH4_H2S(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegnAe_H2S_SO4(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_H2S_SO4(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_NormConst_H2S_SO4(i), i = 1, NumRegnsSedimentDiagenesis)   !Eq. 9.6
    Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_PON_Lab(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_PON_Ref(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_PON_Ine(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_POC_Lab(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_POC_Ref(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_POC_Ine(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_POP_Lab(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_POP_Ref(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_POP_Ine(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_PON_Lab(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_PON_Ref(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_PON_Ine(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_POC_Lab(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_POC_Ref(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_POC_Ine(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_POP_Lab(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_POP_Ref(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_POP_Ine(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, o2h2s          ! cb 5/22/15   moved from below
    Read(CEMAFilN,*)MessageTemp, o2ch4          ! cb 5/22/15
    Read(CEMAFilN,*)MessageTemp, Kdp1        ! cb 5/22/15
    Read(CEMAFilN,*)MessageTemp, Kdp2        ! cb 5/22/15
    Read(CEMAFilN,*)MessageTemp, (SDRegn_CH4CompMethod(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, NH4_NH3_Eqb_Const
    Read(CEMAFilN,*)MessageTemp, HS_H2S_Eqb_Const
    Read(CEMAFilN,*)MessageTemp, HenryConst_NH3
    Read(CEMAFilN,*)MessageTemp, HenryConst_CH4
    Read(CEMAFilN,*)MessageTemp, HenryConst_H2S
    Read(CEMAFilN,*)MessageTemp, HenryConst_CO2
    Read(CEMAFilN,*)MessageTemp, GasDiff_Sed    ! in m^2/s
    Read(CEMAFilN,*)MessageTemp, CalibParam_R1
    Read(CEMAFilN,*)MessageTemp, YoungModulus
    Read(CEMAFilN,*)MessageTemp, CritStressIF
    Read(CEMAFilN,*)MessageTemp, BubbRelScale
    Read(CEMAFilN,*)MessageTemp, CrackCloseFraction
    Read(CEMAFilN,*)MessageTemp, LimBubbSize
    Read(CEMAFilN,*)MessageTemp, MaxBubbRad
    Read(CEMAFilN,*)MessageTemp, UseReleaseFraction
    Read(CEMAFilN,*)MessageTemp, BubbRelFraction
    Read(CEMAFilN,*)MessageTemp, BubbAccFraction
    Read(CEMAFilN,*)MessageTemp, NumBubRelArr
    Read(CEMAFilN,*)MessageTemp, BubbRelFractionAtm
    Read(CEMAFilN,*)MessageTemp, BubbWatGasExchRate
    Read(CEMAFilN,*)MessageTemp, ApplyBubbTurb
    Read(CEMAFilN,*)MessageTemp, CEMATurbulenceScaling
    Read(CEMAFilN,*)MessageTemp, CoeffA_Turb
    Read(CEMAFilN,*)MessageTemp, CoeffB_Turb
    Read(CEMAFilN,*)MessageTemp, WriteCEMAMFTSedFlx
    Read(CEMAFilN,*)MessageTemp, PartMixVel
    Read(CEMAFilN,*)MessageTemp, BurialVel
    Read(CEMAFilN,*)MessageTemp, IncludeIron        ! cb 5/22/15  moved from above
    Read(CEMAFilN,*)MessageTemp, IncludeManganese   ! cb 5/22/15
    Read(CEMAFilN,*)MessageTemp, KFeOOH_HalfSat
    Read(CEMAFilN,*)MessageTemp, kfe_red
    Read(CEMAFilN,*)MessageTemp, kfe_oxid
    Read(CEMAFilN,*)MessageTemp, FeSetVel
    Read(CEMAFilN,*)MessageTemp, KdFe1
    Read(CEMAFilN,*)MessageTemp, KdFe2
    Read(CEMAFilN,*)MessageTemp, o2fe2
    Read(CEMAFilN,*)MessageTemp, KMnO2_HalfSat
    Read(CEMAFilN,*)MessageTemp, kMn_red
    Read(CEMAFilN,*)MessageTemp, kMn_oxid
    Read(CEMAFilN,*)MessageTemp, MnSetVel
    Read(CEMAFilN,*)MessageTemp, KdMn1
    Read(CEMAFilN,*)MessageTemp, KdMn2
    Read(CEMAFilN,*)MessageTemp, o2Mn2
    Read(CEMAFilN,*)MessageTemp, CEMA_POM_Resuspension_Processes
    Read(CEMAFilN,*)MessageTemp, (SDRegn_POMResuspMethod(i), i = 1, NumRegnsSedimentDiagenesis)
    Read(CEMAFilN,*)MessageTemp, TAUCRPOM
    Read(CEMAFilN,*)MessageTemp, crshields
    Read(CEMAFilN,*)MessageTemp, cao_method
    Read(CEMAFilN,*)MessageTemp, spgrav_POM
    Read(CEMAFilN,*)MessageTemp, dia_POM
    READ(CEMAFilN,*,END=100)MessageTemp, SEDIAGFREQ   ! FREQUENCY OF OUTPUT SW 5/25/2017
    READ(CEMAFilN,*,END=100)MessageTemp, BUBBLES_CALCULATION   ! BUBBLES_CALCULATION SW 5/25/2017
    GO TO 101
100 CONTINUE
    ! IF INPUT FILE DOES NOT HAVE NEW LAST LINE - ERROR TRAPPING FOR NEW CODE
    SEDIAGFREQ=7.0
    BUBBLES_CALCULATION=.FALSE.
101 close(CEMAFilN)
    
    !Allocate other variables
    
    !Allocate(CEMAMFT_RandC_RegN(IMX), CEMAMFT_InCond_RegN(IMX), MFTSedFlxVars(IMX,47), CEMA_SD_Vars(IMX,21))
    Allocate(CEMAMFT_RandC_RegN(IMX), CEMAMFT_InCond_RegN(IMX), MFTSedFlxVars(IMX,48), CEMA_SD_Vars(IMX,22))  ! cb 7/26/18
    Allocate(SD_NO3p2(2), SD_NH3p2(2), SD_NH3Tp2(2), SD_CH4p2(2), SD_PO4p2(2), SD_PO4Tp2(2))
	Allocate(SD_HSp2(2), SD_HSTp2(2), SD_POC22(3), SD_PON22(3), SD_POP22(3))
	Allocate(SD_poc2(3), SD_pon2(3), SD_pop2(3), SD_NH3Tp(2), SD_NO3p(2), SD_PO4Tp(2), SD_HSTp(2))
	Allocate(SD_fpon(3), SD_fpoc(3), SD_kdiaPON(3), SD_ThtaPON(3), SD_kdiaPOC(3), SD_ThtaPOC(3))
	Allocate(SD_JPOC(3), SD_JPON(3), SD_JPOP(3), SD_NH3(2),SD_TIC(2),SD_ALK(2),SD_PH(2),SD_TDS(2))
    Allocate(SD_EPOC(3), SD_EPON(3), SD_EPOP(3))
    Allocate(SD_T(2))
	!Allocate(SD_Denit(2), SD_JDenit(2), SD_JO2NO3(2), SD_CH4(2), SD_HS(2), SD_PO4T(2))
    Allocate(SD_Denit(2), SD_JDenit(2), SD_JO2NO3(2), SD_CH4(2), SD_HS(2), SD_PO4T(2), sd_so4(2))   ! cb 7/26/18
    Allocate(SD_Fe2T(2), SD_Fe2(2),SD_FeOOH(2),SD_Mn2T(2), SD_Mn2(2),SD_MnO2(2))
	Allocate(SD_kdiaPOP(3), SD_ThtaPOP(3), SD_NH3T(2), SD_NO3(2), SD_HST(2), SD_PO4(2), SD_FPOP(3))
	!Allocate(SD_SO4Conc(IMX), SD_pHValue(IMX))
    Allocate(SD_pHValue(IMX))   ! cb 7/26/18
	Allocate(SD_AerLayerThick(IMX))
	Allocate(H2SDis(IMX), H2SGas(IMX), CH4Dis(IMX), CH4Gas(IMX))
	Allocate(NH4Dis(IMX), NH4Gas(IMX), CO2Dis(IMX), CO2Gas(IMX))
	Allocate(BubbleRadiusSed(IMX),PresBubbSed(IMX), PresCritSed(IMX))
	Allocate(CgSed(IMX), C0Sed(IMX), CtSed(IMX))
	Allocate(TConc(NumGas,IMX), TConcP(NumGas,IMX), SConc(NumGas,IMX))
	Allocate(DissolvedGasSediments(NumGas,IMX))
	Allocate(CrackOpen(IMX), MFTBubbReleased(IMX), LastDiffVolume(IMX))
	Allocate(BubblesCarried(IMX,NumBubRelArr), BubblesRadius(IMX,NumBubRelArr))
	Allocate(BubblesLNumber(IMX,NumBubRelArr),BubblesStatus(IMX,NumBubRelArr))
	Allocate(BubblesRiseV(IMX,NumBubRelArr))
	Allocate(BubblesGasConc(IMX,NumBubRelArr,NumGas))
	Allocate(BRVoluAGas(IMX,NumBubRelArr,NumGas), BRRateAGas(IMX,NumBubRelArr,NumGas))
	Allocate(FirstBubblesRelease(IMX,NumBubRelArr), BubblesReleaseAllValue(IMX,NumBubRelArr),BubbleRelWB(NWB,NumGas))    ! SW 7/1/2017
	Allocate(BRRateAGasNet(IMX, NumGas))
	Allocate(BubblesAtSurface(IMX,NumBubRelArr))
	Allocate(BottomTurbulence(IMX))
	!Allocate(IceQSS(IMX)) Put in main code sw 6/17/15
	
	!Sediment Generic BOD
	Allocate(SedGenBODConc(NumGenBODConstituents,IMX,2), SedGenBODDecayRate(NumGenBODConstituents,IMX))
	Allocate(SedGenBODConsRate(NumGenBODConstituents,IMX), SedGenBODConsTCoeff(NumGenBODConstituents,IMX))
	Allocate(SedGenBODSS(NumGenBODConstituents),SDPFLUX(NWB),SDNH4FLUX(NWB),SDNO3FLUX(NWB))
		
    SD_NO3p2 = 0.d00; SD_NH3p2 = 0.d00; SD_NH3Tp2 = 0.d00; SD_CH4p2 = 0.d00; BubbleRelWB=0.0   ! SW 7/1/2017
    SD_PO4p2 = 0.d00; SD_PO4Tp2 = 0.d00; SD_HSp2 = 0.d00; SD_HSTp2 = 0.d00; 
    SD_POC22 = 0.d00; SD_PON22 = 0.d00; SD_POP22 = 0.d00; 
    SD_poc2 = 0.d00; SD_pon2 = 0.d00; SD_pop2 = 0.d00; SD_NH3Tp = 0.d00; 
    SD_NO3p = 0.d00; SD_PO4Tp = 0.d00; SD_HSTp = 0.d00; 
    SD_fpon = 0.d00; SD_fpoc = 0.d00; SD_kdiaPON = 0.d00; SD_ThtaPON = 0.d00; 
    SD_kdiaPOC = 0.d00; SD_ThtaPOC = 0.d00; 
    SD_JPOC = 0.d00; SD_JPON = 0.d00; SD_JPOP = 0.d00; SD_NH3 = 0.d00; 
    SD_EPOC = 0.d00; SD_EPON = 0.d00; SD_EPOP = 0.d00
    SD_Denit = 0.d00; SD_JDenit = 0.d00; SD_JO2NO3 = 0.d00; SD_CH4 = 0.d00; 
    SD_HST = 0.d00; SD_PO4 = 0.d00; SD_FPOP = 0.d00; SD_HS = 0.d00; SD_PO4T = 0.d00; 
    SD_Fe2 = 0.d00; SD_Fe2T = 0.d00; SD_FeOOH = 0.d00
    SD_Mn2 = 0.d00; SD_Mn2T = 0.d00; SD_MnO2 = 0.d00; SD_T = 0.d00
    SD_kdiaPOP = 0.d00; SD_ThtaPOP = 0.d00; SD_NH3T = 0.d00; SD_NO3 = 0.d00;
    SD_AerLayerThick = 0.d00
    H2SDis = 0.d00; H2SGas = 0.d00; CH4Dis = 0.d00; CH4Gas = 0.d00; 
    NH4Dis = 0.d00; NH4Gas = 0.d00; CO2Dis = 0.d00; CO2Gas = 0.d00; 
    BubbleRadiusSed = 0.d00; PresBubbSed = 0.d00; PresCritSed = 0.d00; 
    CgSed = 0.d00; C0Sed = 0.d00; CtSed = 0.d00; 
    TConc = 0.d00; TConcP = 0.d00; SConc = 0.d00
    CrackOpen = .FALSE.
    MFTBubbReleased = 0; LastDiffVolume = 0.d00
    BubblesCarried = 0; BubblesLNumber = 0; BubblesStatus = 0
    BubblesRadius = 0.d00; BubblesRiseV = 0.d00; BubblesGasConc = 0.d00
    CEMAMFT_RandC_RegN = 0
    MFTSedFlxVars = 0.d00
    CEMA_SD_Vars = 0.d00
    BubblesReleaseAllValue = 0.d00
    !SD_SO4Conc = 0.d0
    BRVoluAGas = 0.d00; BRRateAGas = 0.d00; BRRateAGasNet = 0.d00
    BottomTurbulence = 0.d00
    !IceQSS = 0.d00  Put in main code SW 6/17/15
    SedGenBODConc = 0.d00; SedGenBODConsRate = 0.d00
    SedGenBODDecayRate = 0.d00; SedGenBODConsTCoeff = 0.d00 
    SedGenBODSS = 0.d00;
    DissolvedGasSediments = 0.d00
    FeSetVel=FeSetVel/86400.0;  kfe_red=kfe_red/86400.0; kfe_oxid=kfe_red/86400.0
    MnSetVel=MnSetVel/86400.0;  kMn_red=kfe_red/86400.0; kMn_oxid=kfe_red/86400.0
    SDPFLUX=0.0
    SDNH4FLUX=0.0
    SDNO3FLUX=0.0
        
    FirstTimeInBubbles = .TRUE. 
    FirstBubblesRelease = .TRUE. 
    BubblesAtSurface = .FALSE.
    
    BODtestout=.TRUE.  ! BOD testing
! special input for CEMA code testing - start

!open(1399,file= 'seddatin.txt',status='old', form='formatted ' )  !Data input 

!read(1399,*)
!read(1399,*) ninp 

!do i = 1,ninp 
!read(1399,*)
!read(1399,*) monzz,dayzz,yearzz,sodzz,jnh4zz,jno3zz,jpo4zz,jsizz,o2zz,nh4zz,no3zz,po4zz,sizz,tempzz   !(B) ! NOTE: sodzz, jno3zz, jpo4zz, and jsizz NOT USED
!end do


!open(1199 ,file='sedin.txt' ,status='old',form='formatted')   !Parameter input 

!read(1199,*)
!read(1199,*)
!read(1199,*) m1,m2,Dp,w2,Dd ,thtaDp,thtaDd                     !(A) 
!read(1199,*)
!read(1199,*)
!read(1199,*) kappnh4 ,pienh4 ,thtanh4 ,kmnh4,thtakmnh4 ,kmnh4o2 
!read(1199,*)
!read(1199,*)
!read(1199,*) kapp1no3,k2no3,thtano3 
!read(1199,*)
!read(1199,*)
!read(1199,*) kappd1,kappp1,pie1s,pie2s,thtapd1,kmhso2 
!read(1199,*)
!read(1199,*)
!read(1199,*) ksi,csisat,dpie1si,pie2si 
!read(1199,*)
!read(1199,*)
!read(1199,*) h2ss,thtasi,kmpsi,o2critsi 
!read(1199,*)
!read(1199,*)
!read(1199,*) dpie1po4,pie2po4,o2crit,kmo2Dp 
!read(1199,*)
!read(1199,*)
!read(1199,*) frpon1,kpon1,thtapon1 
!read(1199,*)
!read(1199,*)
!read(1199,*) frpon2,kpon2,thtapon2 
!read(1199,*)
!read(1199,*)
!read(1199,*) frpon3,kpon3,thtapon3 
!read(1199,*)
!read(1199,*)
!read(1199,*) frpoc1,kpoc1,thtapoc1 
!read(1199,*)
!read(1199,*)
!read(1199,*) frpoc2,kpoc2,thtapoc2 
!read(1199,*)
!read(1199,*)
!read(1199,*) frpoc3,kpoc3,thtapoc3 
!read(1199,*)
!read(1199,*)
!read(1199,*) frpop1,kpop1,thtapop1 
!read(1199,*)
!read(1199,*)
!read(1199,*) frpop2,kpop2,thtapop2 
!read(1199,*)
!read(1199,*)
!read(1199,*) frpop3,kpop3,thtapop3 
!read(1199,*)
!read(1199,*)
!read(1199,*) ratiocn,ratiocp,ratiocsi

!  CEMA testing end
   
    Return

    End Subroutine
    
    Subroutine Deallocate_CEMA
        Use CEMAVars
        IMPLICIT NONE
    
    DEALLOCATE(ConsolidationType,ConstConsolidRate)
	DEALLOCATE(ConstPoreWtrRate,ConsolidRateRegnFil)
	DEALLOCATE(ConsRegSegSt, ConsRegSegEn,BubbleRelWB)     ! SW 7/1/2017

    DEALLOCATE(BedElevation, BedElevationLayer, BedPorosity)
    !DEALLOCATE (sedcellwidth)
    DEALLOCATE(ConsolidRegnNum, BedConsolidRate, PorewaterRelRate)
    DEALLOCATE(CEMASedConc)
    DEALLOCATE(CEMACumPWRelease, CEMALayerAdded, CEMASSApplied)
    DEALLOCATE(CEMACumPWToRelease,CEMACumPWReleased)
    DEALLOCATE(NumCEMAPWInst)
    DEALLOCATE(ApplyCEMAPWRelease)
    DEALLOCATE(CEMACumPWReleaseRate)
    DEALLOCATE(EndBedConsolidation)
    DEALLOCATE(BedConsolidationSeg)  ! cb 6/28/18
    DEALLOCATE(CEMATSSCopy)
    DEALLOCATE(VOLCEMA)

    DEALLOCATE(FFTActPrdSt, FFTActPrdEn)
    DEALLOCATE(FFTLayConc)

    DEALLOCATE(SDRegnPOC_T, SDRegnPON_T, SDRegnSul_T)
    DEALLOCATE(SDRegnPOP_T)
    DEALLOCATE(SDRegnH2S_T, SDRegnNH3_T, SDRegnCH4_T)
    DEALLOCATE(SDRegnTIC_T,SDRegnALK_T, SDRegnPO4_T)
    DEALLOCATE(SDRegnFe2_T,SDRegnFeOOH_T)
    DEALLOCATE(SDRegnMn2_T,SDRegnMnO2_T)
    DEALLOCATE(SDRegnT_T)
    DEALLOCATE(SedBedInitRegSegSt, SedBedInitRegSegEn)
        DEALLOCATE(SedGenBODName)
        DEALLOCATE(SedGenBODInit)
        DEALLOCATE(SedGenBODRegSegSt, SedGenBODRegSegEn)      
        DEALLOCATE(SedGenBODConsRegSegSt, SedGenBODConsRegSegEn)
        DEALLOCATE(SedGenBODRegnRate, SedGenBODRegnTCoeff)
    
    DEALLOCATE(SDRegnPOC_L_Fr, SDRegnPOC_R_Fr, SDRegnPON_L_Fr)
    DEALLOCATE(SDRegnPON_R_Fr, SDRegnPW_DiffCoeff, SDRegnOx_Threshold)
    DEALLOCATE(SDRegnPOP_L_Fr, SDRegnPOP_R_Fr)
    DEALLOCATE(SDRegnAe_NH3_NO3_L, SDRegnAe_NH3_NO3_H, SDRegnAe_NO3_N2_L)
    DEALLOCATE(SDRegnAe_NO3_N2_H, SDRegnAn_NO3_N2, SDRegnAe_CH4_CO2)
    DEALLOCATE(SDRegnAe_HS_NH4_Nit, SDRegnAe_HS_O2_Nit, SDRegn_Theta_PW)
    DEALLOCATE(SDRegn_Theta_NH3_NO3, SDRegn_Theta_NO3_N2, SDRegn_Theta_CH4_CO2)
    DEALLOCATE(SDRegn_Sulfate_CH4_H2S, SDRegnAe_H2S_SO4, SDRegn_Theta_H2S_SO4)
    DEALLOCATE(SDRegn_NormConst_H2S_SO4, SDRegn_MinRate_PON_Lab, SDRegn_MinRate_PON_Ref)
    DEALLOCATE(SDRegn_MinRate_PON_Ine, SDRegn_MinRate_POC_Lab, SDRegn_MinRate_POC_Ref)
    DEALLOCATE(SDRegn_MinRate_POC_Ine, SDRegn_Theta_PON_Lab, SDRegn_Theta_PON_Ref)
    DEALLOCATE(SDRegn_Theta_PON_Ine, SDRegn_Theta_POC_Lab, SDRegn_Theta_POC_Ref)
    DEALLOCATE(SDRegn_Theta_POC_Ine, SDRegn_CH4CompMethod, SDRegn_POMResuspMethod)
    DEALLOCATE(SDRegn_Theta_POP_Lab, SDRegn_Theta_POP_Ref, SDRegn_Theta_POP_Ine)
    DEALLOCATE(SDRegn_MinRate_POP_Lab, SDRegn_MinRate_POP_Ref, SDRegn_MinRate_POP_Ine)
    DEALLOCATE(SedBedDiaRCRegSegSt, SedBedDiaRCRegSegEn)

    DEALLOCATE(CEMAMFT_RandC_RegN, CEMAMFT_InCond_RegN, MFTSedFlxVars, CEMA_SD_Vars)
    DEALLOCATE(SD_NO3p2, SD_NH3p2, SD_NH3Tp2, SD_CH4p2, SD_PO4p2, SD_PO4Tp2)
	DEALLOCATE(SD_HSp2, SD_HSTp2, SD_POC22, SD_PON22, SD_POP22)
	DEALLOCATE(SD_poc2, SD_pon2, SD_pop2, SD_NH3Tp, SD_NO3p, SD_PO4Tp, SD_HSTp)
	DEALLOCATE(SD_fpon, SD_fpoc, SD_kdiaPON, SD_ThtaPON, SD_kdiaPOC, SD_ThtaPOC)
	DEALLOCATE(SD_JPOC, SD_JPON, SD_JPOP, SD_NH3,SD_TIC,SD_ALK,SD_PH,SD_TDS)
    DEALLOCATE(SD_EPOC, SD_EPON, SD_EPOP)
    DEALLOCATE(SD_T)
	!DEALLOCATE(SD_Denit, SD_JDenit, SD_JO2NO3, SD_CH4, SD_HS, SD_PO4T)
    DEALLOCATE(SD_Denit, SD_JDenit, SD_JO2NO3, SD_CH4, SD_HS, SD_PO4T, sd_so4)   ! cb 7/26/18
    DEALLOCATE(SD_Fe2T, SD_Fe2,SD_FeOOH,SD_Mn2T, SD_Mn2,SD_MnO2)
	DEALLOCATE(SD_kdiaPOP, SD_ThtaPOP, SD_NH3T, SD_NO3, SD_HST, SD_PO4, SD_FPOP)
	!DEALLOCATE(SD_SO4Conc, SD_pHValue)
    DEALLOCATE(SD_pHValue)   ! cb 7/26/18
	DEALLOCATE(SD_AerLayerThick)
	DEALLOCATE(H2SDis, H2SGas, CH4Dis, CH4Gas)
	DEALLOCATE(NH4Dis, NH4Gas, CO2Dis, CO2Gas)
	DEALLOCATE(BubbleRadiusSed,PresBubbSed, PresCritSed)
	DEALLOCATE(CgSed, C0Sed, CtSed)
	DEALLOCATE(TConc, TConcP, SConc)
	DEALLOCATE(DissolvedGasSediments)
	DEALLOCATE(CrackOpen, MFTBubbReleased, LastDiffVolume)
	DEALLOCATE(BubblesCarried, BubblesRadius)
	DEALLOCATE(BubblesLNumber,BubblesStatus)
	DEALLOCATE(BubblesRiseV)
	DEALLOCATE(BubblesGasConc)
	DEALLOCATE(BRVoluAGas, BRRateAGas)
	DEALLOCATE(FirstBubblesRelease, BubblesReleaseAllValue)
	DEALLOCATE(BRRateAGasNet)
	DEALLOCATE(BubblesAtSurface)
	DEALLOCATE(BottomTurbulence)
	!Sediment Generic BOD
	DEALLOCATE(SedGenBODConc, SedGenBODDecayRate)
	DEALLOCATE(SedGenBODConsRate, SedGenBODConsTCoeff)
	DEALLOCATE(SedGenBODSS,SDPFLUX,SDNH4FLUX,SDNO3FLUX)
    
    
    RETURN
    End Subroutine Deallocate_CEMA