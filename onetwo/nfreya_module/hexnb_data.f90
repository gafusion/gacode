    MODULE hexnb_data 

    USE nrtype,                                                  ONLY : DP,I4B

      INTEGER(I4B),PARAMETER :: ms = 21, mc = 35
      INTEGER(I4B),PARAMETER :: mz =  1, mi =  2
      INTEGER(I4B),PARAMETER ::   mzz=4, mt1=  4, mt2= 3 , mt3=2



      REAL(DP)       er0, v0, te, ti, ami(mi), deni(mi), amz(mz), &
                    denz(mz), zcor(mz), zsqcor(mz), dene            ! b0
    
      INTEGER(I4B)  ns, nc, numi, numz, iz(mz), izstrp(mz)          ! b0
 
      INTEGER(I4B)  kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl   ! b1

      INTEGER(I4B)  nouthx,istart,ihxbug                            ! b2

      REAL(DP)      f(ms,mc),ar(ms+1,ms+1)                          ! b3

      REAL(DP)       en(mc+1),dg(mc+1),ae(ms,mc),be(ms,mc),       &
                    de1(ms,mc),de2(ms,mc),ge1(ms,mc),ge2(ms,mc)     ! b4

      REAL(DP)      al(ms+1)                                        ! b5

      REAL(DP)      cii(ms+1,mi),cei(ms+1,ms+2,mi),ciz(ms+1,mz),  &
                    cez(ms+1,ms+2,mz),cie(ms+1),cee(ms+1,ms+2),   &
                    ccxi(ms+1,mi),ccxz(ms+1,mz)                     ! b6

      REAL(DP)      q(ms+1,ms+1)                                    ! b7

      INTEGER(I4B)  iexcit,ilorent,mstate,ncont                     ! b8

      REAL(DP)      eigvr(ms+1,ms+1), eigvl(ms+1),cj(ms+1)          ! b9
      
      REAL(DP) xfrac                                                ! b10

      REAL(DP) arad(8,15,mz)                                        ! locrad

      REAL(DP)           A1(mt1,mt2,mt3), A2(mt1,mt2,mt3,mzz)

      INTEGER(I4B)       nth1, nth2, nth3, ntz1(mzz), ntz2(mzz), ntz3(mzz)

    END     MODULE hexnb_data
 
