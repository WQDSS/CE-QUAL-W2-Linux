Subroutine CEMAComputeTurbidity

Use GLOBAL
Use SCREENC
Use CEMAVars

! Type declarations
IMPLICIT NONE

Real(8) CellTSSValue
Integer(2) K

    Do JW=1, NWB
        KT = KTWB(JW)
        Do JB=BS(JW),BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            Do SegNumI = IU, ID
                Do K = KT, KB(SegNumI)
                    !CellTSSValue = C1(K,SegNumI,6)
                    CellTSSValue = C1(K,SegNumI,nturb)      ! cb 2/18/13
                    !C1(K,SegNumI,6) = exp(CoeffA_Turb*log(CellTSSValue) + CoeffB_Turb)
                    !C2(K,SegNumI,6) = exp(CoeffA_Turb*log(CellTSSValue) + CoeffB_Turb)
                    C1(K,SegNumI,nturb) = exp(CoeffA_Turb*log(CellTSSValue) + CoeffB_Turb)       ! cb 2/18/13
                    C2(K,SegNumI,nturb) = exp(CoeffA_Turb*log(CellTSSValue) + CoeffB_Turb)
                End Do !K                 
            End Do                    
        End Do                
    End Do              


Return
End Subroutine