Subroutine InitializeBedConsolidationFiles(TempFilNum, TempFilName)
    Use MAIN 
    Use GLOBAL
    Use SCREENC
    Use CEMAVars
    
    Implicit None
    
    Logical SkipLoop
    Integer(4) TempFilNum
    Character(256) MessageTemp, TempFilName
    
    !Open File	
	Open(TempFilNum, File = TempFilName(1:len_trim(TempFilName)-1))
	
	!Read Header
	SkipLoop = .FALSE.
	Do While(.NOT. SkipLoop)
		Read(TempFilNum,'(a)')MessageTemp
		If(index(MessageTemp, "$") == 0)SkipLoop = .TRUE.
	End Do
			
End Subroutine

Subroutine ReadBedConsolidationFiles(TempFilNum, ConsolidRateTemp)
    Use MAIN 
    Use GLOBAL
    Use SCREENC
    Use CEMAVars
    
    Logical SkipLoop
    Integer(4) TempFilNum
    Real(8) TimeJD1, ConsolidRateTemp1
    Real(8) TimeJD2, ConsolidRateTemp2
    Real(8) FactorInterp, ConsolidRateTemp
    
    !Read Data
	SkipLoop = .FALSE.
	Do While(.NOT. SkipLoop .or. EOF(TempFilNum))
		
	    Read(TempFilNum,'(2F8.0)')TimeJD1, ConsolidRateTemp1
	    Read(TempFilNum,'(2F8.0)')TimeJD2, ConsolidRateTemp2

	    If(JDay >= TimeJD1 .and. JDay <= TimeJD2)then
		    SkipLoop = .TRUE.
	    End If
	    BackSpace(TempFilNum)

	End Do
	
	BackSpace(TempFilNum)
	
	FactorInterp = (JDay - TimeJD1)/(TimeJD2 - TimeJD1)

	ConsolidRateTemp = ConsolidRateTemp1*(1-FactorInterp) + ConsolidRateTemp2*FactorInterp
	
End Subroutine