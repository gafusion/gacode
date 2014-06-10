 

     SUBROUTINE copy_injector_attributes(orig_injector_no,duplicate_injector_no)
!---------------------------------------------------------------------------------------
! -- create a pseudo injector, identified as duplicate_injector_no, that is
! -- a copy of a real injector,identified as orig_injector_no
!---------------------------------------------------------------------------------------
     USE nrtype,                                 ONLY : DP,I4B

     USE neutral_beams,                          ONLY  :                                         &
                                                         anglev, angleh,nashape, aheigh,         &
                                                         awidth, bcur, bptor, blenp, nbshape,    &
                                                         bleni,bheigh, bwidth, bhfoc, bvfoc,     &
                                                         bhdiv, bvdiv, ebkev, fbcur, nbeams,     &
                                                         naptr, alen, bvofset, bhofset, nsourc,  &
                                                         sfrac1, npart, npskip, rpivot,          &
                                                         zpivot,pseudo_injectors  
                                                         
                                       
                                       



     IMPLICIT NONE 

     INTEGER(I4B) orig_injector_no,duplicate_injector_no


     !--------------------------------------------------------------------------- 
     ! we assign the same power and current to the duplicate as the original
     ! the duplicate injector is treated as an independent statistical sample of
     ! the original ijector when results are collected
     !---------------------------------------------------------------------------
     bcur(duplicate_injector_no)       =       bcur(orig_injector_no)     
     bptor(duplicate_injector_no)      =       bptor(orig_injector_no) 
   
     
     anglev(duplicate_injector_no)     =       anglev(orig_injector_no)   
     angleh(duplicate_injector_no)     =       angleh(orig_injector_no)   
     bvofset(duplicate_injector_no)    =       bvofset(orig_injector_no)  
     bhofset(duplicate_injector_no)    =       bhofset(orig_injector_no)  
     bleni(duplicate_injector_no)      =       bleni(orig_injector_no) 
     nbshape(duplicate_injector_no)    =       nbshape(orig_injector_no)  
     bheigh(duplicate_injector_no)     =       bheigh(orig_injector_no)   
     bwidth(duplicate_injector_no)     =       bwidth(orig_injector_no)   
     bhdiv(duplicate_injector_no)      =       bhdiv(orig_injector_no)    
     bvdiv(duplicate_injector_no)      =       bvdiv(orig_injector_no)    
     fbcur(1,duplicate_injector_no)    =       fbcur(1,orig_injector_no)  
     fbcur(2,duplicate_injector_no)    =       fbcur(2,orig_injector_no)  
     fbcur(3,duplicate_injector_no)    =       fbcur(3,orig_injector_no)  
     bhfoc(duplicate_injector_no)      =       bhfoc(orig_injector_no)    
     bvfoc(duplicate_injector_no)      =       bvfoc(orig_injector_no)    
     ebkev(duplicate_injector_no)      =       ebkev(orig_injector_no)    
     ebkev(duplicate_injector_no)      =       ebkev(orig_injector_no)    
     sfrac1(duplicate_injector_no)     =       sfrac1(orig_injector_no)   
     nashape(1,duplicate_injector_no)  =       nashape(1,orig_injector_no)
     nashape(2,duplicate_injector_no)  =       nashape(2,orig_injector_no)
     nashape(3,duplicate_injector_no)  =       nashape(3,orig_injector_no)
     nashape(4,duplicate_injector_no)  =       nashape(4,orig_injector_no)
     aheigh(1,duplicate_injector_no)   =       aheigh(1,orig_injector_no) 
     awidth(1,duplicate_injector_no)   =       awidth(1,orig_injector_no) 
     alen(1,duplicate_injector_no)     =       alen(1,orig_injector_no)   
     aheigh(2,duplicate_injector_no)   =       aheigh(2,orig_injector_no) 
     awidth(2,duplicate_injector_no)   =       awidth(2,orig_injector_no) 
     alen(2,duplicate_injector_no)     =       alen(2,orig_injector_no)   
     alen(3,duplicate_injector_no)     =       alen(3,orig_injector_no)   
     awidth(4,duplicate_injector_no)   =       awidth(4,orig_injector_no) 
     alen(4,duplicate_injector_no)     =       alen(4,orig_injector_no)   
     blenp(duplicate_injector_no)      =       blenp(orig_injector_no)    
     rpivot(duplicate_injector_no)     =       rpivot(orig_injector_no)       
     zpivot(duplicate_injector_no)     =       zpivot(orig_injector_no)       
     

     RETURN
     END SUBROUTINE copy_injector_attributes
