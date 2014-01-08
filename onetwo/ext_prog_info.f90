       MODULE ext_prog_info
!          
             IMPLICIT NONE
             !number of different architectures, ncpu-arch:
!             INTEGER, PARAMETER :: ncpu_arch = 40 !jmp.ibm
!             INTEGER, PARAMETER :: ncpu_arch = 41  !add benten
!             INTEGER, PARAMETER :: ncpu_arch  = 44  !add lohan5,6,7
!             INTEGER, PARAMETER :: ncpu_arch  = 49  !add lohan8,9,10,11,12
!             INTEGER, PARAMETER :: ncpu_arch  = 50  !add lohan1
!             INTEGER, PARAMETER :: ncpu_arch  = 56  !add star cluster
              INTEGER, PARAMETER :: ncpu_arch  = 57  !add vidar

             INTEGER, PARAMETER :: ncpu_arch_32 = 7
             !number of grid size variants available:
             INTEGER, PARAMETER :: grid_types = 3
             INTEGER, PARAMETER :: grid_types_prplt = 3
             !number of toray verisons available:
             INTEGER, PARAMETER :: n_versions = 3
             REAL  versions(n_versions)
             INTEGER echin_save,nw_rf,nh_rf,kj_rf,nj_rf,nchars_12
             !INTEGER,EXTERNAL :: ihost_name !jmp.ibm.par
             CHARACTER (LEN = 256) :: host_name_ext !jmp.ibm.par
             INTEGER nchr, host_index,nchars_host
             REAL toray_version,toray_version_switch
             REAL, PARAMETER :: vertol = 1.e-3  !version tolerance. toray 
                                  !versions are
                                  !assumed equal if difference
                                  !is less than vertol
             LOGICAL is_32_bit
             CHARACTER (len = 256),DIMENSION(ncpu_arch,n_versions,grid_types) ::  &
                   toray_paths
             CHARACTER *12 storesp
             CHARACTER (len = 256),DIMENSION(ncpu_arch,grid_types) ::  &
                   preplt_paths
             CHARACTER (len = 256) :: toray_path,onetwo_xsct
             CHARACTER (len = 256) :: curray_path
             CHARACTER (len = 256) :: genray_path
             CHARACTER (len = 256) :: toq_path
             CHARACTER (len = 256) :: pedestal_path
             CHARACTER (len = 256) :: nubeam_path,nubeam_setup_ext
             CHARACTER (len = 256),DIMENSION(ncpu_arch,grid_types) :: curray_paths
             CHARACTER (len = 256),DIMENSION(ncpu_arch,grid_types) :: toq_paths
             CHARACTER (len = 256),DIMENSION(ncpu_arch) :: toq_base
             CHARACTER (len = 256),DIMENSION(ncpu_arch) :: nubeam_paths
             CHARACTER (len = 256),DIMENSION(ncpu_arch) :: nubeam_setup
             CHARACTER (len = 256) :: preplt_path
             CHARACTER (len = 256),DIMENSION(ncpu_arch) :: pedestal_paths
             CHARACTER (len = 32) ,DIMENSION(ncpu_arch) :: host_names
             CHARACTER (len = 32) ,DIMENSION(ncpu_arch_32) :: host_names_32
             CHARACTER (len = 16) ::  host_name
             CHARACTER (len = 64)  :: root_str(ncpu_arch)
             CHARACTER (len = 64)  :: root_str_prplt(ncpu_arch)
             CHARACTER (len = 128) :: fastcd_paths( ncpu_arch,grid_types)
             CHARACTER (len = 256) :: fastcd_path  !user specified path
!            plotcod3eid is keyed     on PREPLT and TRPLOT:
!            character (len =8)    :: plotcodeid ='04feb03v'
!            CHARACTER (len =8)    :: plotcodeid ='04feb04v'
!            CHARACTER (len =8)    ::  plotcodeid ='01apr04v'
!             CHARACTER (len =8)   ::  plotcodeid ='14apr06v'
             CHARACTER (len =8)    ::  plotcodeid  ='10dec08v'
             CHARACTER (len =16)   :: gcnmp_hosts(16)
             CHARACTER (len =128)  :: gcnmp_paths(16)
             CHARACTER (len =16)   :: P_Nfreya_hosts(22)
             CHARACTER (len =128)  :: P_Nfreya_paths(22)

             CHARACTER *256  preplt_py ,preplt_to_run
             CHARACTER *256 pedestal_to_run
             DATA host_names                                 &
                  /'VENUS','DELPHI2','LOHAN2','HYDRA ',     &
                   'DELPHI','CARDEA','KATZE ', 'NEMSIS',     &
                   'ULAM','IRENIC','HERA','PHOBOS',          &
                   'USCWS8','HERMES',                        &
                   'AETNA','LOHAN3','RANIER','ZEUS',         &
                   'USCWSD','EOS','PHOEBE','HESTIA','ELCAP', &
                   'AJAX','BOB5','PLUTO','LOHAN4',           &
                   'NFRC','JAGUAR','RFPLASMA','SNUIBM',      & !jmp.ibm
                   'ISIS1','ISIS2',                          & !jmp.ibm
                   'BENTEN','LOHAN5','HEAD','NODE01',        & ! head = 36
                   'NODE02','NODE03','NODE04','NODE05',      &
                   'NODE06','LOHAN1',                        &
                   'STAR','STAR1','STAR2','STAR3','STAR4',   & ! star = 44
                   'STAR5','STAR6','STAR7','STAR8','STAR9',  &
                   'STAR10','STAR11','STAR12', 'VIDAR'/     

             DATA host_names_32                              &
                  / 'LOHAN2','DELPHI','LOHAN3',              &
                  'ZEUS','EOS','PHOEBE','HESTIA'/


             DATA gcnmp_hosts                                &
                  /'VENUS','DELPHI2','LOHAN2','LOHAN3',      &
                   'LOHAN4','LOHAN5','LOHAN6'  ,'LOHAN7',    &
                   'BENTEN','SCYLD','CHARM',                 &
                   'LOHAN8','LOHAN9','LOHAN10','LOHAN11',    &
                   'LOHAN12'                                        /
             ! select 32 or 64 bit version:
             DATA gcnmp_paths(1) / '/p/linux/GCNMP/source/gcnmp'      /  ! VENUS
             DATA gcnmp_paths(2) / '/p/linux/GCNMP/source/32/gcnmp'   /  ! DELPHI2
             DATA gcnmp_paths(3) / '/p/linux/GCNMP/source/32/gcnmp'   /  ! LOHAN2
             DATA gcnmp_paths(4) / '/p/linux/GCNMP/source/32/gcnmp'   /  ! LOHAN3
             DATA gcnmp_paths(5) / '/p/linux/GCNMP/source/gcnmp'      /  ! LOHAN4
             DATA gcnmp_paths(6) / '/p/linux/GCNMP/source/gcnmp'      /  ! LOHAN5
             DATA gcnmp_paths(7) / '/p/linux/GCNMP/source/gcnmp'      /  ! LOHAN6, or HEAD
             DATA gcnmp_paths(8) / '/p/linux/GCNMP/source/gcnmp'      /  ! LOHAN7, or NODE01
             DATA gcnmp_paths(9) / '/p/linux/GCNMP/source/gcnmp'      /  ! BENTEN
             DATA gcnmp_paths(10) / '/home/stjohn/GCNMP/source/gcnmp' /  ! SCYLD
             DATA gcnmp_paths(11) / '/home/stjohn/GCNMP/source/gcnmp' /  ! CHARM
             DATA gcnmp_paths(12) / '/p/linux/GCNMP/source/gcnmp'     /  ! LOHAN8,NODE02 
             DATA gcnmp_paths(13) / '/p/linux/GCNMP/source/gcnmp'     /  ! LOHAN9,NODE03 
             DATA gcnmp_paths(14) / '/p/linux/GCNMP/source/gcnmp'     /  ! LOHAN10,NODE04 
             DATA gcnmp_paths(15) / '/p/linux/GCNMP/source/gcnmp'     /  ! LOHAN11,NODE05 
             DATA gcnmp_paths(16) / '/p/linux/GCNMP/source/gcnmp'     /  ! LOHAN12,NODE06   
             
             DATA P_Nfreya_hosts                              & ! Note only lohan1
                  /'VENUSA',  'LOHAN1', 'LOHAN2',              & ! lohan6 ,venus and  star* are 
                   'DELPHI2','LOHAN3', 'LOHAN4',              & ! parallel nodes
                   'LOHAN5', 'LOHAN6',                        & 
                   'STAR','STAR1','STAR2','STAR3','STAR4',    &
                   'STAR5','STAR6','STAR7','STAR8','STAR9',   &
                   'STAR10','STAR11','STAR12','VIDAR'/
     
             DATA P_Nfreya_paths                                                      &
                  / '/p/linux/onetwo/pgf90/P_Nfreya_ven_star',  &  ! VENUS
                    '/p/linux/onetwo/pgf90/P_Nfreya_l1',        &  ! LOHAN!
                    '/p/linux/onetwo/pgf90/P_Nfreya_l2',        &  ! not defined on l2
                    '/p/linux/onetwo/pgf90/P_Nfreya_l2',        &  ! not defined on delphi
                    '/p/linux/onetwo/pgf90/P_Nfreya_l3',        &  ! not defined on l3
                    '/p/linux/onetwo/pgf90/P_Nfreya_l4',        &  ! not defined on l4
                    '/p/linux/onetwo/pgf90/P_Nfreya_l5',        &  ! not defined on l5
                    '/p/linux/onetwo/pgf90/P_Nfreya_l6',        &  ! LOHAN6
                    '/p/linux/onetwo/pgf90/P_Nfreya_ven_star',  &  ! STAR
                    '/p/linux/onetwo/pgf90/P_Nfreya_ven_star',  &  ! STAR1
                    '/p/linux/onetwo/pgf90/P_Nfreya_ven_star',  &  ! STAR2
                    '/p/linux/onetwo/pgf90/P_Nfreya_ven_star',  &  ! STAR3
                    '/p/linux/onetwo/pgf90/P_Nfreya_ven_star',  &  ! STAR4
                    '/p/linux/onetwo/pgf90/P_Nfreya_ven_star',  &  ! STAR5
                    '/p/linux/onetwo/pgf90/P_Nfreya_ven_star',  &  ! STAR6
                    '/p/linux/onetwo/pgf90/P_Nfreya_ven_star',  &  ! STAR7
                    '/p/linux/onetwo/pgf90/P_Nfreya_ven_star',  &  ! STAR8
                    '/p/linux/onetwo/pgf90/P_Nfreya_ven_star',  &  ! STAR9
                    '/p/linux/onetwo/pgf90/P_Nfreya_ven_star',  &  ! STAR10
                    '/p/linux/onetwo/pgf90/P_Nfreya_ven_star',  &  ! STAR11
                    '/p/linux/onetwo/pgf90/P_Nfreya_ven_star',  &  ! STAR12
                    '/p/linux/onetwo/pgf90/P_Nfreya_l1'         &  ! VIDAR
                  /




!                  For NFRC,SNUIBM,JAGUAR,RFPLASMA,ISIS1,ISIS2
!                  ONLY onetwo_xsct and preplt_paths are defined 

             DATA toray_version_switch /1.80/

!             DATA genray_path       / '' / !user specified path for genray 
             DATA genray_path       / '/task/imd/genray/xgenray' / ! HSJ 2/14/2013

                   



             DATA fastcd_paths(1,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'          / !VENUS
             DATA fastcd_paths(2,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'          / !DELPHI2
             DATA fastcd_paths(3,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'          / !LOHAN2
             DATA fastcd_paths(4,1) / '/u/stjohn/fastcd/hp/v1.0/129_101/fastcd'     / !HYDRA
             DATA fastcd_paths(5,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'          / !DELPHI
             DATA fastcd_paths(6,1) / '/u/stjohn/fastcd/osf/v1.0/129_101/fastcd'    / !CARDEA
             DATA fastcd_paths(7,1) / '/u/stjohn/fastcd/osf/v1.0/129_101/fastcd'    / !KATZE
             DATA fastcd_paths(8,1) / '/u/stjohn/fastcd/osf/v1.0/129_101/fastcd'    / !NEMSIS
             DATA fastcd_paths(9,1) / '/u/stjohn/fastcd/osf/v1.0/129_101/fastcd'    / !ULAM
             DATA fastcd_paths(10,1) / '/u/stjohn/fastcd/osf/v1.0/129_101/fastcd'   / !IRENIC
             DATA fastcd_paths(11,1) / '/u/stjohn/fastcd/osf/v1.0/129_101/fastcd'   / !HERA
             DATA fastcd_paths(12,1) / '/u/stjohn/fastcd/osf/v1.0/129_101/fastcd'   / !PHOBOS
             DATA fastcd_paths(13,1) / '/u/stjohn/fastcd/osf/v1.0/129_101/fastcd'   / !USCWS8
             DATA fastcd_paths(14,1) / '/u/stjohn/fastcd/osf/v1.0/129_101/fastcd'   / !HERMES
             DATA fastcd_paths(15,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !AETNA
             DATA fastcd_paths(15,2) / '/p/linux/fastcd/v1.0/lf95/129_101/fastcd'   / !AETNA
             DATA fastcd_paths(16,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !LOHAN3
             DATA fastcd_paths(16,2) / '/p/linux/fastcd/v1.0/lf95/129_101/fastcd'   / !LOHAN3
             DATA fastcd_paths(17,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !RANIER
             DATA fastcd_paths(17,2) / '/p/linux/fastcd/v1.0/lf95/129_101/fastcd'   / !RANIER
             DATA fastcd_paths(18,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !ZEUS
             DATA fastcd_paths(19,1) / '/u/stjohn/fastcd/hp/v1.0/129_101/fastcd'    / !USCWSD
             DATA fastcd_paths(20,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !EOS
             DATA fastcd_paths(21,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !PHOEBE
             DATA fastcd_paths(22,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !HESTIA
             DATA fastcd_paths(23,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !ELCAP
             DATA fastcd_paths(24,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !AJAX
             DATA fastcd_paths(25,1) / 'not yet defined'                            / !BOB5
             DATA fastcd_paths(26,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !PLUTO
             DATA fastcd_paths(27,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !LOHAN4
             DATA fastcd_paths(28,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !NFRC
             DATA fastcd_paths(29,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !JAGUAR
             DATA fastcd_paths(30,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !RFPLASMA
             DATA fastcd_paths(31,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !SNUIBM
             DATA fastcd_paths(32,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !ISIS1
             DATA fastcd_paths(33,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !ISIS2
             DATA fastcd_paths(34,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !BENTEN
             DATA fastcd_paths(35,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !LOHAN5
             DATA fastcd_paths(36,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !HEAD
             DATA fastcd_paths(37,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !NODE01
             DATA fastcd_paths(38,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !NODE02
             DATA fastcd_paths(39,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !NODE03
             DATA fastcd_paths(40,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !NODE04
             DATA fastcd_paths(41,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !NODE05
             DATA fastcd_paths(42,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !NODE06
             DATA fastcd_paths(43,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !LOHAN1
             DATA fastcd_paths(44,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !STAR
             DATA fastcd_paths(45,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !STAR1
             DATA fastcd_paths(46,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !STAR2
             DATA fastcd_paths(47,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !STAR3
             DATA fastcd_paths(48,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !STAR4
             DATA fastcd_paths(49,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !STAR5
             DATA fastcd_paths(50,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !STAR6
             DATA fastcd_paths(51,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !STAR7
             DATA fastcd_paths(52,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !STAR8
             DATA fastcd_paths(53,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !STAR9
             DATA fastcd_paths(54,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !STAR10
             DATA fastcd_paths(55,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !STAR11
             DATA fastcd_paths(56,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !STAR12
             DATA fastcd_paths(57,1) / '/p/linux/fastcd/v1.0/129_51/fastcd'         / !VIDAR






             DATA fastcd_path       / '' /  !user specified path for fastcd


             DATA curray_paths(1,1)   /     '/p/linux/curray/65_65/xcurray'          / !VENUS
             DATA curray_paths(1,2)   /     '/p/linux/curray/65_65/xcurray'          / !VENUS
             DATA curray_paths(2,1)   /     '/p/linux/curray/129_129/xcurray'        / !DELPHI2
             DATA curray_paths(2,2)   /     '/p/linux/curray/65_65/xcurray'          / !DELPHI2

             DATA curray_paths(3,1)   /     '/p/linux/curray/129_129/xcurray'        / !LOHAN2
             DATA curray_paths(3,2)   /     '/p/linux/curray/65_65/xcurray'          / !LOHAN2

             DATA curray_paths(4,1)   /     '/u/stjohn/curray/hp/129_129/xcurray'    / !HP
             DATA curray_paths(4,2)   /     '/u/stjohn/curray/hp/65_65/xcurray'      / !HP

             DATA curray_paths(5,1)   /     '/p/linux/curray/129_129/xcurray'        / !DELPHI
             DATA curray_paths(5,2)   /     '/p/linux/curray/65_65/xcurray'          / !DELPHI

             DATA curray_paths(6,1)   /     '/u/stjohn/curray/osf/129_129/xcurray'   / !CARDEA
             DATA curray_paths(6,2)   /     '/u/stjohn/curray/osf/65_65/xcurray'     / !CARDEA

             DATA curray_paths(7,1)   /     '/u/stjohn/curray/osf/129_129/xcurray'   / !KATZE
             DATA curray_paths(7,2)   /     '/u/stjohn/curray/osf/65_65/xcurray'     / !KATZE

             DATA curray_paths(8,1)   /     '/u/stjohn/curray/osf/129_129/xcurray'   / !NEMSIS  
             DATA curray_paths(8,2)   /     '/u/stjohn/curray/osf/65_65/xcurray'     / !NEMSIS  

             DATA curray_paths(9,1)   /     '/u/stjohn/curray/osf/129_129/xcurray'   / !ULAM 
             DATA curray_paths(9,2)   /     '/u/stjohn/curray/osf/65_65/xcurray'     / !ULAM
     
             DATA curray_paths(10,1)   /     '/u/stjohn/curray/osf/129_129/xcurray'  /!IRENIC   
             DATA curray_paths(10,2)   /     '/u/stjohn/curray/osf/65_65/xcurray'    /!IRENIC  
  
             DATA curray_paths(11,1)   /     '/u/stjohn/curray/osf/129_129/xcurray'  /!HERA
             DATA curray_paths(11,2)   /     '/u/stjohn/curray/osf/65_65/xcurray'    /!HERA
  
             DATA curray_paths(12,1)   /     '/u/stjohn/curray/osf/129_129/xcurray'  /!PHOBOS
             DATA curray_paths(12,2)   /     '/u/stjohn/curray/osf/65_65/xcurray'    /!PHOBOS

             DATA curray_paths(13,1)   /     '/u/stjohn/curray/osf/129_129/xcurray'  /!USCWS8
             DATA curray_paths(13,2)   /     '/u/stjohn/curray/osf/65_65/xcurray'    /!USCWS8

             DATA curray_paths(14,1)   /     '/u/stjohn/curray/osf/129_129/xcurray'  /!HERMES
             DATA curray_paths(14,2)   /     '/u/stjohn/curray/osf/65_65/xcurray'    /!HERMES

             DATA curray_paths(15,1)   /     '/p/linux/curray/129_129/xcurray'        / !AETNA
             DATA curray_paths(15,2)   /     '/p/linux/curray/65_65/xcurray'          / !AETNA

             DATA curray_paths(16,1)   /     '/p/linux/curray/129_129/xcurray'        / !LOHAN3
             DATA curray_paths(16,2)   /     '/p/linux/curray/65_65/xcurray'          / !LOHAN3

             DATA curray_paths(17,1)   /     '/p/linux/curray/129_129/xcurray'        / !RANIER
             DATA curray_paths(17,2)   /     '/p/linux/curray/65_65/xcurray'          / !RANIER

             DATA curray_paths(18,1)   /     '/p/linux/curray/129_129/xcurray'        / !ZEUS
             DATA curray_paths(18,2)   /     '/p/linux/curray/65_65/xcurray'          / !ZEUS


             DATA curray_paths(19,1)   /     '/u/stjohn/curray/hp/129_129/xcurray'    / !USCWSD
             DATA curray_paths(19,2)   /     '/u/stjohn/curray/hp/65_65/xcurray'      / !USCWSD

             DATA curray_paths(20,1)   /     '/p/linux/curray/129_129/xcurray'        / !EOS
             DATA curray_paths(20,2)   /     '/p/linux/curray/65_65/xcurray'          / !EOS

             DATA curray_paths(21,1)   /     '/p/linux/curray/129_129/xcurray'        / !PHOEBE
             DATA curray_paths(21,2)   /     '/p/linux/curray/65_65/xcurray'          / !PHOEBE

             DATA curray_paths(22,1)   /     '/p/linux/curray/129_129/xcurray'        / !HESTIA
             DATA curray_paths(22,2)   /     '/p/linux/curray/65_65/xcurray'           / !HESTIA


             DATA curray_paths(23,1)   /     '/p/linux/curray/129_129/xcurray'        /  !ELCAP
             DATA curray_paths(23,2)   /     '/p/linux/curray/65_65/xcurray'           / !ELCAP

             DATA curray_paths(24,1)   /     '/p/linux/curray/129_129/xcurray'        /  !AJAX
             DATA curray_paths(24,2)   /     '/p/linux/curray/65_65/xcurray'           / !AJAX
 

             DATA curray_paths(25,1)   /     'not yet defined '                         / !BOB5
             DATA curray_paths(25,2)   /     'not yet defined '                         / !BOB5


             DATA curray_paths(26,1)   /     '/p/linux/curray/129_129/xcurray'        / !PLUTO
             DATA curray_paths(26,2)   /     '/p/linux/curray/65_65/xcurray'          / !PLUTO


             DATA curray_paths(27,1)   /     '/p/linux/curray/129_129/xcurray'        / !LOHAN4
             DATA curray_paths(27,2)   /     '/p/linux/curray/65_65/xcurray'          / !LOHAN4
             DATA curray_paths(28,1)   /     '/p/linux/curray/129_129/xcurray'        / !NFRC
             DATA curray_paths(28,2)   /     '/p/linux/curray/65_65/xcurray'          / !NFRC

             DATA curray_paths(29,1)   /     '/p/linux/curray/129_129/xcurray'        / !JAGUAR
             DATA curray_paths(29,2)   /     '/p/linux/curray/65_65/xcurray'          / !JAGUAR
             DATA curray_paths(30,1)   /     '/p/linux/curray/129_129/xcurray'        / !RFPLASMA
             DATA curray_paths(30,2)   /     '/p/linux/curray/65_65/xcurray'          / !RFPLASMA
             DATA curray_paths(31,1)   /     '/p/linux/curray/129_129/xcurray'        / !SNUIBM
             DATA curray_paths(31,2)   /     '/p/linux/curray/65_65/xcurray'          / !SNUIBM
             DATA curray_paths(32,1)   /     '/p/linux/curray/129_129/xcurray'        / !ISIS1
             DATA curray_paths(32,2)   /     '/p/linux/curray/65_65/xcurray'          / !ISIS1
             DATA curray_paths(33,1)   /     '/p/linux/curray/129_129/xcurray'        / !ISIS2
             DATA curray_paths(33,2)   /     '/p/linux/curray/65_65/xcurray'          / !ISIS2

             DATA curray_paths(34,1)   /     '/p/linux/curray/129_129/xcurray'        / !BENTEN
             DATA curray_paths(34,2)   /     '/p/linux/curray/65_65/xcurray'          / !BENTEN
             DATA curray_paths(35,1)   /     '/p/linux/curray/129_129/xcurray'        / !LOHAN5
             DATA curray_paths(35,2)   /     '/p/linux/curray/65_65/xcurray'          / !LOHAN5

             DATA curray_paths(36,1)   /     '/p/linux/curray/129_129/xcurray'        / !HEAD
             DATA curray_paths(36,2)   /     '/p/linux/curray/65_65/xcurray'          / !HEAD
             DATA curray_paths(37,1)   /     '/p/linux/curray/129_129/xcurray'        / !node01
             DATA curray_paths(37,2)   /     '/p/linux/curray/65_65/xcurray'          / !node01
             DATA curray_paths(38,1)   /     '/p/linux/curray/129_129/xcurray'        / !node02
             DATA curray_paths(38,2)   /     '/p/linux/curray/65_65/xcurray'          / !node02
             DATA curray_paths(39,1)   /     '/p/linux/curray/129_129/xcurray'        / !node03
             DATA curray_paths(39,2)   /     '/p/linux/curray/65_65/xcurray'          / !node03
             DATA curray_paths(40,1)   /     '/p/linux/curray/129_129/xcurray'        / !node04
             DATA curray_paths(40,2)   /     '/p/linux/curray/65_65/xcurray'          / !node04
             DATA curray_paths(41,1)   /     '/p/linux/curray/129_129/xcurray'        / !node05
             DATA curray_paths(41,2)   /     '/p/linux/curray/65_65/xcurray'          / !node05
             DATA curray_paths(42,1)   /     '/p/linux/curray/129_129/xcurray'        / !node06
             DATA curray_paths(42,2)   /     '/p/linux/curray/65_65/xcurray'          / !node06
             DATA curray_paths(43,1)   /     '/p/linux/curray/129_129/xcurray'        / !LOHAN1
             DATA curray_paths(43,2)   /     '/p/linux/curray/65_65/xcurray'          / !LOHAN1

             DATA curray_paths(44,1)   /     '/p/linux/curray/65_65/xcurray'          / !STAR
             DATA curray_paths(44,2)   /     '/p/linux/curray/65_65/xcurray'          / !STAR
             DATA curray_paths(45,1)   /     '/p/linux/curray/65_65/xcurray'          / !STAR1
             DATA curray_paths(45,2)   /     '/p/linux/curray/65_65/xcurray'          / !STAR1
             DATA curray_paths(46,1)   /     '/p/linux/curray/65_65/xcurray'          / !STAR2
             DATA curray_paths(46,2)   /     '/p/linux/curray/65_65/xcurray'          / !STAR2
             DATA curray_paths(47,1)   /     '/p/linux/curray/65_65/xcurray'          / !STAR3
             DATA curray_paths(47,2)   /     '/p/linux/curray/65_65/xcurray'          / !STAR3
             DATA curray_paths(48,1)   /     '/p/linux/curray/65_65/xcurray'          / !STAR4
             DATA curray_paths(48,2)   /     '/p/linux/curray/65_65/xcurray'          / !STAR4
             DATA curray_paths(49,1)   /     '/p/linux/curray/65_65/xcurray'          / !STAR5
             DATA curray_paths(49,2)   /     '/p/linux/curray/65_65/xcurray'          / !STAR5
             DATA curray_paths(50,1)   /     '/p/linux/curray/65_65/xcurray'          / !STAR6
             DATA curray_paths(50,2)   /     '/p/linux/curray/65_65/xcurray'          / !STAR6
             DATA curray_paths(51,1)   /     '/p/linux/curray/65_65/xcurray'          / !STAR7
             DATA curray_paths(51,2)   /     '/p/linux/curray/65_65/xcurray'          / !STAR7
             DATA curray_paths(52,1)   /     '/p/linux/curray/65_65/xcurray'          / !STAR8
             DATA curray_paths(52,2)   /     '/p/linux/curray/65_65/xcurray'          / !STAR8
             DATA curray_paths(53,1)   /     '/p/linux/curray/65_65/xcurray'          / !STAR9
             DATA curray_paths(53,2)   /     '/p/linux/curray/65_65/xcurray'          / !STAR9
             DATA curray_paths(54,1)   /     '/p/linux/curray/65_65/xcurray'          / !STAR10
             DATA curray_paths(54,2)   /     '/p/linux/curray/65_65/xcurray'          / !STAR10
             DATA curray_paths(55,1)   /     '/p/linux/curray/65_65/xcurray'          / !STAR11
             DATA curray_paths(55,2)   /     '/p/linux/curray/65_65/xcurray'          / !STAR11
             DATA curray_paths(56,1)   /     '/p/linux/curray/65_65/xcurray'          / !STAR12
             DATA curray_paths(56,2)   /     '/p/linux/curray/65_65/xcurray'          / !STAR12
             DATA curray_paths(57,2)   /     '/p/linux/curray/65_65/xcurray'          / !VIDAR
 
             DATA curray_path       / '' / !user specified path for curray




             !VENUS   data:
             !DATA  nubeam_paths(1) /'/p/linux/nubeam/pgf90/LINUX/test/nubeam_driver' /
             !DATA  nubeam_setup(1) /'/p/linux/nubeam/.nubeam' /

             !VENUS   data (new nubeam)
             DATA  nubeam_paths(1) /'/p/linux/nubeam/201201/bin/d3d_nubeam_driver' /
             DATA  nubeam_setup(1) /'none' /

             !delphi2 data:
             DATA nubeam_paths(2) /'/p/linux/nubeam/jan2006/nbdrive/bin/nubeam_driver' /
             DATA  nubeam_setup(2) /'/p/linux/nubeam/.nubeam' /

             !lohan2 data:
             DATA nubeam_paths(3) /'/p/linux/nubeam/jan2006/nbdrive/bin/nubeam_driver' / 
             DATA  nubeam_setup(3) /'/p/linux/nubeam/.nubeam' /


             !hydra data:
             DATA nubeam_paths(4) /'none'/
             DATA  nubeam_setup(4) /'none'/


             !delphi data:
             DATA nubeam_paths(5)  /'/p/linux/nubeam/pgf90/LINUX/test/nubeam_driver' / 
             DATA  nubeam_setup(5) /'/p/linux/nubeam/.nubeam' /

             !cardea data:
             DATA nubeam_paths(6)  /'none'/
             DATA  nubeam_setup(6) /'none'/

             !katze data:
             DATA nubeam_paths(7) /'none'/
             DATA  nubeam_setup(7) /'none'/

             !nemsis  data:
             DATA nubeam_paths(8)  /'none'/
             DATA  nubeam_setup(8) /'none'/

             !ulam data:
             DATA nubeam_paths(9)  /'none'/
             DATA  nubeam_setup(9) /'none'/

             !irenic data:
             DATA nubeam_paths(10)  /'none'/
             DATA  nubeam_setup(10) /'none'/

             !hera data:
             DATA nubeam_paths(11)  /'none'/
             DATA  nubeam_setup(11) /'none'/

             !phobos data:
             DATA nubeam_paths(12)   /'none'/
             DATA  nubeam_setup(12)  /'none'/

             !uscws8 data:
             DATA nubeam_paths(13)   /'none'/
             DATA  nubeam_setup(13)  /'none'/


             !hermes data:
             DATA nubeam_paths(14)  /'none'/
             DATA  nubeam_setup(14) /'none'/


             !aetna data:
             DATA nubeam_paths(15)   /'/p/linux/nubeam/nubeam_lf95_rh3/nubeam_lf95/transp_beam_lf95/LINUX/test/nubeam_driver'/
             DATA  nubeam_setup(15) /'/p/linux/nubeam/nubeam_lf95_rh3/nubeam_lf95/transp_beam_lf95/LINUX/test/.nubeam' /

             !lohan3 data:
             DATA nubeam_paths(16)  /'/p/linux/nubeam/jan2006/nbdrive/bin/nubeam_driver'/
             DATA  nubeam_setup(16) /'/p/linux/nubeam/.nubeam' /

             !ranier data:
             DATA nubeam_paths(17)  /'/home/stjohn/transport_codes/transp_beam_pgf90/LINUX/test/nubeam_driver'/
             DATA  nubeam_setup(17) /'none'/

             !ZEUS data:
             DATA nubeam_paths(18) /'/p/linux/nubeam/pgf90/LINUX/test/nubeam_driver' /
             DATA  nubeam_setup(18) /'/p/linux/nubeam/.nubeam' /


             !USCWSD data:
             DATA nubeam_paths(19) /'none'/
             DATA  nubeam_setup(19) /'none'/

             !EOS data:
             DATA nubeam_paths(20) /'/p/linux/nubeam/pgf90/LINUX/test/nubeam_driver' /
             DATA  nubeam_setup(20) /'/p/linux/nubeam/.nubeam' /

             !PHOEBE data:
             DATA nubeam_paths(21) /'/p/linux/nubeam/pgf90/LINUX/test/nubeam_driver' /
             DATA  nubeam_setup(21) /'/p/linux/nubeam/.nubeam' /

             !HESTIA data:
             DATA nubeam_paths(22) /'/p/linux/nubeam/pgf90/LINUX/test/nubeam_driver' /
             DATA  nubeam_setup(22) /'/p/linux/nubeam/.nubeam' /


             !ELCAP  data:
             DATA nubeam_paths(23) /'/p/linux/nubeam/pgf90/LINUX/test/nubeam_driver' /
             DATA  nubeam_setup(23) /'/p/linux/nubeam/.nubeam' /



             !AJAX  data:
             DATA nubeam_paths(24) /'/p/linux/nubeam/pgf90/LINUX/test/nubeam_driver' /
             DATA  nubeam_setup(24) /'/p/linux/nubeam/.nubeam' /

             !BOB5 data:
             DATA nubeam_paths(25) /'not yet defined' /
             DATA  nubeam_setup(25) /'not yet defined' /

             !PLUTO data:
             DATA nubeam_paths(26) /'/p/linux/nubeam/pgf90/LINUX/test/nubeam_driver' /
             DATA  nubeam_setup(26) /'/p/linux/nubeam/.nubeam' /

             !lohan4 data:
             !DATA nubeam_paths(28) /'/p/linux/nubeam/jan2006/pgf90_6.1_64/nubeam/LINUX/test/nubeam_driver
             !Not sure if above 64 bit nubeam works. Use 32 bit version:
             DATA nubeam_paths(27) /'/p/linux/nubeam/jan2006/nbdrive/bin/nubeam_driver' /
             DATA nubeam_setup(27) /'/p/linux/nubeam/.nubeam' /

             !NFRCdata:
             DATA nubeam_paths(28) /'none'/
             DATA  nubeam_setup(28) /'none'/   
             !JAGUAR         
             DATA nubeam_paths(29) /'none'/
             DATA  nubeam_setup(29) /'none'/
             !RFPLASMA     
             DATA nubeam_paths(30) /'none'/
             DATA  nubeam_setup(30) /'none'/
             !SNUIBM     
             DATA nubeam_paths(31) /'none'/
             DATA  nubeam_setup(31) /'none'/ 
             !ISIS1     
             DATA nubeam_paths(32) /'none'/
             DATA  nubeam_setup(32) /'none'/ 
             !ISIS2  
             DATA nubeam_paths(33) /'none'/
             DATA  nubeam_setup(33) /'none'/

             ! following machines are wired to new nubeam 
             !BENTEN
             DATA nubeam_paths(34) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(34) / 'none'/
             !LOHAN5
             DATA nubeam_paths(35) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(35) / 'none'/

             !HEAD
             DATA nubeam_paths(36) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(36) / 'none'/
             !node01
             DATA nubeam_paths(37) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(37) / 'none'/
             !node02
             DATA nubeam_paths(38) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(38) / 'none'/
             !node03
             DATA nubeam_paths(39) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(39) / 'none'/
             !node04
             DATA nubeam_paths(40) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(40) / 'none'/
             !node05
             DATA nubeam_paths(41) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(41) / 'none'/
             !node06
             DATA nubeam_paths(42) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(42) / 'none'/

             !LOHAN1
             DATA nubeam_paths(43) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(43) / 'none'/

             !STAR
             DATA nubeam_paths(44) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(44) / 'none'/
             !STAR1
             DATA nubeam_paths(45) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(45) / 'none'/
             !STAR2
             DATA nubeam_paths(46) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(46) / 'none'/
             !STAR3
             DATA nubeam_paths(47) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(47) / 'none'/
             !STAR4
             DATA nubeam_paths(48) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(48) / 'none'/
             !STAR5
             DATA nubeam_paths(49) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(49) / 'none'/
             !STAR6
             DATA nubeam_paths(50) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(50) / 'none'/
             !STAR7
             DATA nubeam_paths(51) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(51) / 'none'/
             !STAR8
             DATA nubeam_paths(52) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(52) / 'none'/
             !STAR9
             DATA nubeam_paths(53) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(53) / 'none'/
             !STAR10
             DATA nubeam_paths(54) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(54) / 'none'/
             !STAR11
             DATA nubeam_paths(55) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(55) / 'none'/
             !STAR12
             DATA nubeam_paths(56) / '/p/linux/nubeam/201201/bin/d3d_nubeam_driver'/
             DATA nubeam_setup(56) / 'none'/
             !VIDAR   data (new nubeam)
             DATA  nubeam_paths(57) /'/p/linux/nubeam/201201/bin/d3d_nubeam_driver' /
             DATA  nubeam_setup(57) /'none' /





             !NOTE: fixb129xy.sc is stored in toq_base(host_index)
             !The protocol for external file usage in Onetwo
             !was broken by the Toq team. Consequently we have
             !paths to the fixed boundary code that toq uses
             !to generate efit type eqdsks stored in 
             !the script fixb129xy.sc   --- HSJ 03/18/05 ----

             DATA toq_paths(1,1)   /     '/p/linux/toq/129_129/toq.x'        / !VENUS
             DATA toq_paths(1,2)   /     '/p/linux/toq/65_65/toq.x'          / !VENUS
             DATA toq_base(1)      /     '/p/linux/toq/129_129/'             / !VENUS

             DATA toq_paths(2,1)   /     '/p/linux/toq/129_129/toq.x'        / !DELPHI2
             DATA toq_paths(2,2)   /     '/p/linux/toq/65_65/toq.x'          / !DELPHI2
             DATA toq_base(2)      /     '/p/linux/toq/129_129/'             / !DELPHI2

             DATA toq_paths(3,1)   /     '/p/linux/toq/129_129/toq.x'        / !LOHAN2
             DATA toq_paths(3,2)   /     '/p/linux/toq/65_65/toq.x'          / !LOHAN2
             DATA toq_base(3)      /     '/p/linux/toq/129_129/'             / !LOHAN2



             DATA toq_paths(4,1)   /     '/u/stjohn/toq/hp/129_129/toq.x'    / !HP
             DATA toq_paths(4,2)   /     '/u/stjohn/toq/hp/65_65/toq.x'      / !HP
             DATA toq_base(4)      /     'not supported '                    / !HP

             DATA toq_paths(5,1)   /     '/p/linux/toq/129_129/toq.x'        / !DELPHI
             DATA toq_paths(5,2)   /     '/p/linux/toq/65_65/toq.x'          / !DELPHI
             DATA toq_base(5)      /     '/p/linux/toq/129_129/'             / !DELPHI


             DATA toq_paths(6,1)   /     '/u/stjohn/toq/osf/129_129/toq.x'   / !CARDEA
             DATA toq_paths(6,2)   /     '/u/stjohn/toq/osf/65_65/toq.x'     / !CARDEA
             DATA toq_base(6)      /     'not supported '                    / !CARDEA

             DATA toq_paths(7,1)   /     '/u/stjohn/toq/osf/129_129/toq.x'   / !KATZE
             DATA toq_paths(7,2)   /     '/u/stjohn/toq/osf/65_65/toq.x'     / !KATZE
             DATA toq_base(7)      /     'not supported '                    / !KATZE

             DATA toq_paths(8,1)   /     '/u/stjohn/toq/osf/129_129/toq.x'   / !NEMSIS  
             DATA toq_paths(8,2)   /     '/u/stjohn/toq/osf/65_65/toq.x'     / !NEMSIS 
             DATA toq_base(8)      /     'not supported '                    / !NEMSIS

             DATA toq_paths(9,1)   /     '/u/stjohn/toq/osf/129_129/toq.x'   / !ULAM 
             DATA toq_paths(9,2)   /     '/u/stjohn/toq/osf/65_65/toq.x'     / !ULAM
             DATA toq_base(9)      /     'not supported '                    / !ULAM
     
             DATA toq_paths(10,1)   /     '/u/stjohn/toq/osf/129_129/toq.x'  /!IRENIC   
             DATA toq_paths(10,2)   /     '/u/stjohn/toq/osf/65_65/toq.x'    /!IRENIC 
             DATA toq_base(10)      /     'not supported '                    /!IRENIC 
  
             DATA toq_paths(11,1)   /     '/u/stjohn/toq/osf/129_129/toq.x'  /!HERA
             DATA toq_paths(11,2)   /     '/u/stjohn/toq/osf/65_65/toq.x'    /!HERA
             DATA toq_base(11)      /      'not supported '                    /!HERA
  
             DATA toq_paths(12,1)   /     '/u/stjohn/toq/osf/129_129/toq.x'  /!PHOBOS
             DATA toq_paths(12,2)   /     '/u/stjohn/toq/osf/65_65/toq.x'    /!PHOBOS
             DATA toq_base(12)      /     'not supported '                    / !PHOBOS

             DATA toq_paths(13,1)   /     '/u/stjohn/toq/osf/129_129/toq.x'  /!USCWS8
             DATA toq_paths(13,2)   /     '/u/stjohn/toq/osf/65_65/toq.x'    /!USCWS8
             DATA toq_base(13)      /     'not supported '                    / !USCWS8

             DATA toq_paths(14,1)   /     '/u/stjohn/toq/osf/129_129/toq.x'  /!HERMES
             DATA toq_paths(14,2)   /     '/u/stjohn/toq/osf/65_65/toq.x'    /!HERMES
             DATA toq_base(14)      /     'not supported '                    / !HERMES

             DATA toq_paths(15,1)   /     '/p/linux/toq/129_129/toq.x'        / !AETNA
             DATA toq_paths(15,2)   /     '/p/linux/toq/65_65/toq.x'          / !AETNA
             DATA toq_base(15)      /     '/p/linux/toq/129_129/'             / !AETNA


             DATA toq_paths(16,1)   /     '/p/linux/toq/129_129/toq.x'        / !LOHAN3
             DATA toq_paths(16,2)   /     '/p/linux/toq/65_65/toq.x'          / !LOHAN3
             DATA toq_base(16)      /     '/p/linux/toq/129_129/'             / !LOHAN3



             DATA toq_paths(17,1)   /     '/p/linux/toq/129_129/toq.x'        / !RANIER
             DATA toq_paths(17,2)   /     '/p/linux/toq/65_65/toq.x'          / !RANIER
             DATA toq_base(17)      /     '/p/linux/toq/129_129/'             / !RANIER


             DATA toq_paths(18,1)   /     '/p/linux/toq/129_129/toq.x'        / !RANIER
             DATA toq_paths(18,2)   /     '/p/linux/toq/65_65/toq.x'          / !RANIER
             DATA toq_base(18)      /     '/p/linux/toq/129_129/'             / !RANIER


             DATA toq_paths(19,1)   /     '/u/stjohn/toq/hp/129_129/toq.x'     / !USCWSD
             DATA toq_paths(19,2)   /     '/u/stjohn/toq/hp/65_65/toq.x'       / !USCWSD
             DATA toq_base(19)      /     ' not suppoerted  '                  / !USCWSD


             DATA toq_paths(20,1)   /     '/p/linux/toq/129_129/toq.x'        / !EOS
             DATA toq_paths(20,2)   /     '/p/linux/toq/65_65/toq.x'          / !EOS
             DATA toq_base(20)      /     '/p/linux/toq/129_129/'             / !EOS


             DATA toq_paths(21,1)   /     '/p/linux/toq/129_129/toq.x'        / !PHOEBE
             DATA toq_paths(21,2)   /     '/p/linux/toq/65_65/toq.x'          / !PHOEBE
             DATA toq_base(21)      /     '/p/linux/toq/129_129/'             / !PHOEBE


             DATA toq_paths(22,1)   /     '/p/linux/toq/129_129/toq.x'        / !HESTIA
             DATA toq_paths(22,2)   /     '/p/linux/toq/65_65/toq.x'          / !HESTIA
             DATA toq_base(22)      /     '/p/linux/toq/129_129/'             / !HESTIA


             DATA toq_paths(23,1)   /     '/p/linux/toq/129_129/toq.x'        / !ELCAP
             DATA toq_paths(23,2)   /     '/p/linux/toq/65_65/toq.x'          / !ELCAP
             DATA toq_base(23)      /     '/p/linux/toq/129_129/'             / !ELCAP



             DATA toq_paths(24,1)   /     '/p/linux/toq/129_129/toq.x'        / !AJAX
             DATA toq_paths(24,2)   /     '/p/linux/toq/65_65/toq.x'          / !AJAX
             DATA toq_base(24)      /     '/p/linux/toq/129_129/'             / !AJAX


             DATA toq_paths(25,1)   /     'not yet defined'                    / !BOB5
             DATA toq_paths(25,2)   /     'not yet defined'                    / !BOB5
             DATA toq_base(25)      /     'not yet defined'                    / !BOB5


             DATA toq_paths(26,1)   /     '/p/linux/toq/129_129/toq.x'        / !PLUTO
             DATA toq_paths(26,2)   /     '/p/linux/toq/65_65/toq.x'          / !PLUTO
             DATA toq_base(26)      /     '/p/linux/toq/129_129/'             / !PLUTO


             DATA toq_paths(27,1)   /     '/p/linux/toq/129_129/toq.x'        / !LOHAN4
             DATA toq_paths(27,2)   /     '/p/linux/toq/65_65/toq.x'          / !LOHAN4
             DATA toq_base(27)      /     '/p/linux/toq/129_129/'             / !LOHAN4

             DATA toq_paths(28,1)   /     '/p/linux/toq/129_129/toq.x'        / !NFRC
             DATA toq_paths(28,2)   /     '/p/linux/toq/65_65/toq.x'          / !NFRC
             DATA toq_base(28)      /     '/p/linux/toq/129_129/'             / !NFRC

             DATA toq_paths(29,1)   /     '/p/linux/toq/129_129/toq.x'        / !JAGUAR
             DATA toq_paths(29,2)   /     '/p/linux/toq/65_65/toq.x'          / !JAGUAR
             DATA toq_base(29)      /     '/p/linux/toq/129_129/'             / !JAGUAR


             DATA toq_paths(30,1)   /     '/p/linux/toq/129_129/toq.x'        / !RFPLASMA
             DATA toq_paths(30,2)   /     '/p/linux/toq/65_65/toq.x'          / !RFPLASMA
             DATA toq_base(30)      /     '/p/linux/toq/129_129/'             / !RFPLASMA

             DATA toq_paths(31,1)   /     '/p/linux/toq/129_129/toq.x'        / !SNUIBM
             DATA toq_paths(31,2)   /     '/p/linux/toq/65_65/toq.x'          / !SNUIBM
             DATA toq_base(31)      /     '/p/linux/toq/129_129/'             / !SNUIBM

             DATA toq_paths(32,1)   /     '/p/linux/toq/129_129/toq.x'        / !ISIS1
             DATA toq_paths(32,2)   /     '/p/linux/toq/65_65/toq.x'          / !ISIS1
             DATA toq_base(32)      /     '/p/linux/toq/129_129/'             / !ISIS1

             DATA toq_paths(33,1)   /     '/p/linux/toq/129_129/toq.x'        / !ISIS2
             DATA toq_paths(33,2)   /     '/p/linux/toq/65_65/toq.x'          / !ISIS2
             DATA toq_base(33)      /     '/p/linux/toq/129_129/'             / !ISIS2


             DATA toq_paths(34,1)   /     '/p/linux/toq/129_129/toq.x'        / !BENTEN
             DATA toq_paths(34,2)   /     '/p/linux/toq/65_65/toq.x'          / !BENTEN
             DATA toq_base(34)      /     '/p/linux/toq/129_129/'             / !BENTEN


             DATA toq_paths(35,1)   /     '/p/linux/toq/129_129/toq.x'        / !LOHAN5
             DATA toq_paths(35,2)   /     '/p/linux/toq/65_65/toq.x'          / !LOHAN5
             DATA toq_base(35)      /     '/p/linux/toq/129_129/'             / !LOHAN5


             DATA toq_paths(36,1)   /     '/p/linux/toq/129_129/toq.x'        / !HEAD
             DATA toq_paths(36,2)   /     '/p/linux/toq/65_65/toq.x'          / !HEAD
             DATA toq_base(36)      /     '/p/linux/toq/129_129/'             / !HEAD
             DATA toq_paths(37,1)   /     '/p/linux/toq/129_129/toq.x'        / !NODE01
             DATA toq_paths(37,2)   /     '/p/linux/toq/65_65/toq.x'          / !NODE01
             DATA toq_base(37)      /     '/p/linux/toq/129_129/'             / !NODE01
             DATA toq_paths(38,1)   /     '/p/linux/toq/129_129/toq.x'        / !NODE02
             DATA toq_paths(38,2)   /     '/p/linux/toq/65_65/toq.x'          / !NODE02
             DATA toq_base(38)      /     '/p/linux/toq/129_129/'             / !NODE02
             DATA toq_paths(39,1)   /     '/p/linux/toq/129_129/toq.x'        / !NODE03
             DATA toq_paths(39,2)   /     '/p/linux/toq/65_65/toq.x'          / !NODE03
             DATA toq_base(39)      /     '/p/linux/toq/129_129/'             / !NODE03
             DATA toq_paths(40,1)   /     '/p/linux/toq/129_129/toq.x'        / !NODE04
             DATA toq_paths(40,2)   /     '/p/linux/toq/65_65/toq.x'          / !NODE04
             DATA toq_base(40)      /     '/p/linux/toq/129_129/'             / !NODE04
             DATA toq_paths(41,1)   /     '/p/linux/toq/129_129/toq.x'        / !NODE05
             DATA toq_paths(41,2)   /     '/p/linux/toq/65_65/toq.x'          / !NODE05
             DATA toq_base(41)      /     '/p/linux/toq/129_129/'             / !NODE05
             DATA toq_paths(42,1)   /     '/p/linux/toq/129_129/toq.x'        / !NODE06
             DATA toq_paths(42,2)   /     '/p/linux/toq/65_65/toq.x'          / !NODE06
             DATA toq_base(42)      /     '/p/linux/toq/129_129/'             / !NODE06

             DATA toq_paths(43,1)   /     '/p/linux/toq/129_129/toq.x'        / !LOHAN1
             DATA toq_paths(43,2)   /     '/p/linux/toq/65_65/toq.x'          / !LOHAN1
             DATA toq_base(43)      /     '/p/linux/toq/129_129/'             / !LOHAN1


             DATA toq_paths(44,1)   /     '/p/linux/toq/129_129/toq.x'        / !STAR
             DATA toq_paths(44,2)   /     '/p/linux/toq/65_65/toq.x'          / !STAR
             DATA toq_base(44)      /     '/p/linux/toq/129_129/'             / !STAR
             DATA toq_paths(45,1)   /     '/p/linux/toq/129_129/toq.x'        / !STAR1
             DATA toq_paths(45,2)   /     '/p/linux/toq/65_65/toq.x'          / !STAR1
             DATA toq_base(45)      /     '/p/linux/toq/129_129/'             / !STAR1

             DATA toq_paths(46,1)   /     '/p/linux/toq/129_129/toq.x'        / !STAR2
             DATA toq_paths(46,2)   /     '/p/linux/toq/65_65/toq.x'          / !STAR2
             DATA toq_base(46)      /     '/p/linux/toq/129_129/'             / !STAR2

             DATA toq_paths(47,1)   /     '/p/linux/toq/129_129/toq.x'        / !STAR3
             DATA toq_paths(47,2)   /     '/p/linux/toq/65_65/toq.x'          / !STAR3
             DATA toq_base(47)      /     '/p/linux/toq/129_129/'             / !STAR3

             DATA toq_paths(48,1)   /     '/p/linux/toq/129_129/toq.x'        / !STAR4
             DATA toq_paths(48,2)   /     '/p/linux/toq/65_65/toq.x'          / !STAR4
             DATA toq_base(48)      /     '/p/linux/toq/129_129/'             / !STAR4

             DATA toq_paths(49,1)   /     '/p/linux/toq/129_129/toq.x'        / !STAR5
             DATA toq_paths(49,2)   /     '/p/linux/toq/65_65/toq.x'          / !STAR5
             DATA toq_base(49)      /     '/p/linux/toq/129_129/'             / !STAR5

             DATA toq_paths(50,1)   /     '/p/linux/toq/129_129/toq.x'        / !STAR6
             DATA toq_paths(50,2)   /     '/p/linux/toq/65_65/toq.x'          / !STAR6
             DATA toq_base(50)      /     '/p/linux/toq/129_129/'             / !STAR6

             DATA toq_paths(51,1)   /     '/p/linux/toq/129_129/toq.x'        / !STAR7
             DATA toq_paths(51,2)   /     '/p/linux/toq/65_65/toq.x'          / !STAR7
             DATA toq_base(51)      /     '/p/linux/toq/129_129/'             / !STAR7


             DATA toq_paths(52,1)   /     '/p/linux/toq/129_129/toq.x'        / !STAR8
             DATA toq_paths(52,2)   /     '/p/linux/toq/65_65/toq.x'          / !STAR8
             DATA toq_base(52)      /     '/p/linux/toq/129_129/'             / !STAR8


             DATA toq_paths(53,1)   /     '/p/linux/toq/129_129/toq.x'        / !STAR9
             DATA toq_paths(53,2)   /     '/p/linux/toq/65_65/toq.x'          / !STAR9
             DATA toq_base(53)      /     '/p/linux/toq/129_129/'             / !STAR9


             DATA toq_paths(54,1)   /     '/p/linux/toq/129_129/toq.x'        / !STAR10
             DATA toq_paths(54,2)   /     '/p/linux/toq/65_65/toq.x'          / !STAR10
             DATA toq_base(54)      /     '/p/linux/toq/129_129/'             / !STAR10

             DATA toq_paths(55,1)   /     '/p/linux/toq/129_129/toq.x'        / !STAR11
             DATA toq_paths(55,2)   /     '/p/linux/toq/65_65/toq.x'          / !STAR11
             DATA toq_base(55)      /     '/p/linux/toq/129_129/'             / !STAR11

             DATA toq_paths(56,1)   /     '/p/linux/toq/129_129/toq.x'        / !STAR12
             DATA toq_paths(56,2)   /     '/p/linux/toq/65_65/toq.x'          / !STAR12
             DATA toq_base(56)      /     '/p/linux/toq/129_129/'             / !STAR12
             DATA toq_paths(57,1)   /     '/p/linux/toq/129_129/toq.x'        / !VIDAR
             DATA toq_paths(57,2)   /     '/p/linux/toq/65_65/toq.x'          / !VIDAR
             DATA toq_base(57)      /     '/p/linux/toq/129_129/'             / !VIDAR
             DATA toq_path       / '' / !user specified path for toq
             DATA is_32_bit      /.false./

      CONTAINS

       SUBROUTINE set_rf_data
             !------------------------------------------------------------------------
             !it should be fairly obvious  how to modify this
             !the assumed properties (in sub get_toray) are
             !given here:
             ! (1) all hosts must have n_versions  versions available:
              versions = (/ 0.97,1.70,1.80 /)

             !(2) last three fields in toray_paths 
             !(ie fields after root_str) must follow the
             !root_str(j)//'/vx.xx/grid/toray' 
             !pattern otherwise sub get_toray will fail !!!!!!
             !(3) obviously if any of the parameters 
             !ncpu_arch,n_versions,grid_types are 
             !changed then the following data statements must be
             !changed accordingly. 
             !(I do not have  the luxury of creating code generators
             !to do such things - HSJ )



             !VENUS  data 1 :

             toray_paths(1,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(1,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(1,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(1,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(1,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(1,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(1,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(1,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(1,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(1) = '/p/linux/toray' 


              


             !DELPHI2 data 2 :

             toray_paths(2,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(2,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(2,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(2,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(2,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(2,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(2,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(2,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(2,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(2) = '/p/linux/toray' 

             !lohan2 data 3 :
             toray_paths(3,1,1)=   toray_paths(2,1,1)
             toray_paths(3,2,1)=   toray_paths(2,2,1)
             toray_paths(3,3,1)=   toray_paths(2,3,1)
             toray_paths(3,1,2)=   toray_paths(2,1,2)
             toray_paths(3,2,2)=   toray_paths(2,2,2)
             toray_paths(3,3,2)=   toray_paths(2,3,2)
             toray_paths(3,1,3)=   toray_paths(2,1,3)
             toray_paths(3,2,3)=   toray_paths(2,2,3)
             toray_paths(3,3,3)=   toray_paths(2,3,3)
             root_str(3)       =   root_str(2) 


             !hydra data 4 :
             toray_paths(4,1,1)= '/u/stjohn/toray/hp/v1.60/129_51/toray'
             toray_paths(4,2,1)= '/u/stjohn/toray/hp/v0.97/129_51/toray' 
             toray_paths(4,3,1)= '/u/stjohn/toray/hp/v1.70/129_51/toray' 
             toray_paths(4,1,2)= '/u/stjohn/toray/hp/v1.60/65_51/toray'
             toray_paths(4,2,2)= '/u/stjohn/toray/hp/v0.97/65_51/toray' 
             toray_paths(4,3,2)= '/u/stjohn/toray/hp/v1.70/65_51/toray' 
             toray_paths(4,1,3)= '/u/stjohn/toray/hp/v1.60/129_201/toray'
             toray_paths(4,2,3)= '/u/stjohn/toray/hp/v0.97/129_201/toray' 
             toray_paths(4,3,3)= '/u/stjohn/toray/hp/v1.70/129_201/toray'   
             root_str(4) = '/u/stjohn/toray/hp' 

             !delphi data 5 :
             toray_paths(5,1,1)=   toray_paths(2,1,1)
             toray_paths(5,2,1)=   toray_paths(2,2,1)
             toray_paths(5,3,1)=   toray_paths(2,3,1)
             toray_paths(5,1,2)=   toray_paths(2,1,2)
             toray_paths(5,2,2)=   toray_paths(2,2,2)
             toray_paths(5,3,2)=   toray_paths(2,3,2)
             toray_paths(5,1,3)=   toray_paths(2,1,3)
             toray_paths(5,2,3)=   toray_paths(2,2,3)
             toray_paths(5,3,3)=   toray_paths(2,3,3)
             root_str(5)       =   root_str(2) 

             !cardea data 6 : 
             toray_paths(6,1,1)= '/u/stjohn/toray/osf/v1.60/129_51/toray'
             toray_paths(6,2,1)= '/u/stjohn/toray/osf/v0.97/129_51/toray' 
             toray_paths(6,3,1)= '/u/stjohn/toray/osf/v1.70/129_51/toray'  
             toray_paths(6,1,2)= '/u/stjohn/toray/osf/v1.60/65_51/toray'
             toray_paths(6,2,2)= '/u/stjohn/toray/osf/v0.97/65_51/toray'
             toray_paths(6,3,2)= '/u/stjohn/toray/osf/v1.70/65_51/toray' 
             toray_paths(6,1,3)= '/u/stjohn/toray/osf/v1.60/129_201/toray'
             toray_paths(6,2,3)= '/u/stjohn/toray/osf/v0.97/129_201/toray' 
             toray_paths(6,3,3)= '/u/stjohn/toray/osf/v1.70/129_51/toray'  
             root_str(6) = '/u/stjohn/toray/osf'  

             !katze data 7:
             toray_paths(7,1,1)= toray_paths(6,1,1)
             toray_paths(7,2,1)= toray_paths(6,2,1)
             toray_paths(7,3,1)= toray_paths(6,3,1)
             toray_paths(7,1,2)= toray_paths(6,1,2)
             toray_paths(7,2,2)= toray_paths(6,2,2)
             toray_paths(7,3,2)= toray_paths(6,3,2)
             toray_paths(7,1,3)= toray_paths(6,1,3)
             toray_paths(7,2,3)= toray_paths(6,2,3)
             toray_paths(7,3,3)= toray_paths(6,3,3)
             root_str(7)       = root_str(6)

             !nemsis  data 8:
             toray_paths(8,1,1)=  toray_paths(6,1,1)
             toray_paths(8,2,1)=  toray_paths(6,2,1)
             toray_paths(8,3,1)=  toray_paths(6,3,1)
             toray_paths(8,1,2)=  toray_paths(6,1,2)
             toray_paths(8,2,2)=  toray_paths(6,2,2)
             toray_paths(8,3,2)=  toray_paths(6,3,2)
             toray_paths(8,1,3)=  toray_paths(6,1,3)
             toray_paths(8,2,3)=  toray_paths(6,2,3)
             toray_paths(8,3,3)=  toray_paths(6,3,3)
             root_str(8)       =  root_str(6)

             !ulam data 9 :
             toray_paths(9,1,1)= toray_paths(6,1,1)
             toray_paths(9,2,1)= toray_paths(6,2,1)
             toray_paths(9,3,1)= toray_paths(6,3,1)
             toray_paths(9,1,2)= toray_paths(6,1,2)
             toray_paths(9,2,2)= toray_paths(6,2,2)
             toray_paths(9,3,2)= toray_paths(6,3,2)
             toray_paths(9,1,3)= toray_paths(6,1,3)
             toray_paths(9,2,3)= toray_paths(6,2,3)
             toray_paths(9,3,3)= toray_paths(6,3,3)
             root_str(9)       = root_str(6)

             !irenic data 10:
             toray_paths(10,1,1)= toray_paths(6,1,1)
             toray_paths(10,2,1)= toray_paths(6,2,1)
             toray_paths(10,3,1)= toray_paths(6,3,1)
             toray_paths(10,1,2)= toray_paths(6,1,2)
             toray_paths(10,2,2)= toray_paths(6,2,2)
             toray_paths(10,3,2)= toray_paths(6,3,2)
             toray_paths(10,1,3)= toray_paths(6,1,3)
             toray_paths(10,2,3)= toray_paths(6,2,3)
             toray_paths(10,3,3)= toray_paths(6,3,3)
             root_str(10)       = root_str(6)

             !hera data 11 :
             toray_paths(11,1,1)= toray_paths(6,1,1)
             toray_paths(11,2,1)= toray_paths(6,2,1)
             toray_paths(11,3,1)= toray_paths(6,3,1)
             toray_paths(11,1,2)= toray_paths(6,1,2)
             toray_paths(11,2,2)= toray_paths(6,2,2)
             toray_paths(11,3,2)= toray_paths(6,3,2)
             toray_paths(11,1,3)= toray_paths(6,1,3)
             toray_paths(11,2,3)= toray_paths(6,2,3)
             toray_paths(11,3,3)= toray_paths(6,3,3)
             root_str(11)       = root_str(6)

             !phobos data 12 :
             toray_paths(12,1,1)= toray_paths(6,1,1)
             toray_paths(12,2,1)= toray_paths(6,2,1)
             toray_paths(12,3,1)= toray_paths(6,3,1)
             toray_paths(12,1,2)= toray_paths(6,1,2)
             toray_paths(12,2,2)= toray_paths(6,2,2)
             toray_paths(12,3,2)= toray_paths(6,3,2)
             toray_paths(12,1,3)= toray_paths(6,1,3)
             toray_paths(12,2,3)= toray_paths(6,2,3)
             toray_paths(12,3,3)= toray_paths(6,3,3)
             root_str(12)       = root_str(6)


             !uscws8 data 13 :
             toray_paths(13,1,1)= toray_paths(6,1,1)
             toray_paths(13,2,1)= toray_paths(6,2,1)
             toray_paths(13,3,1)= toray_paths(6,3,1)
             toray_paths(13,1,2)= toray_paths(6,1,2)
             toray_paths(13,2,2)= toray_paths(6,2,2)
             toray_paths(13,3,2)= toray_paths(6,3,2)
             toray_paths(13,1,3)= toray_paths(6,1,3)
             toray_paths(13,2,3)= toray_paths(6,2,3)
             toray_paths(13,3,3)= toray_paths(6,3,3)
             root_str(13)       = root_str(6)

             !hermes data 14:
             toray_paths(14,1,1)= toray_paths(6,1,1)
             toray_paths(14,2,1)= toray_paths(6,2,1)
             toray_paths(14,3,1)= toray_paths(6,3,1)
             toray_paths(14,1,2)= toray_paths(6,1,2)
             toray_paths(14,2,2)= toray_paths(6,2,2)
             toray_paths(14,3,2)= toray_paths(6,3,2)
             toray_paths(14,1,3)= toray_paths(6,1,3)
             toray_paths(14,2,3)= toray_paths(6,2,3)
             toray_paths(14,3,3)= toray_paths(6,3,3)
             root_str(14)       = root_str(6)


             !AETNA data 15 :
             toray_paths(15,1,1)=   toray_paths(2,1,1)
             toray_paths(15,2,1)=   toray_paths(2,2,1)
             toray_paths(15,3,1)=   toray_paths(2,3,1)
             toray_paths(15,1,2)=   toray_paths(2,1,2)
             toray_paths(15,2,2)=   toray_paths(2,2,2)
             toray_paths(15,3,2)=   toray_paths(2,3,2)
             toray_paths(15,1,3)=   toray_paths(2,1,3)
             toray_paths(15,2,3)=   toray_paths(2,2,3)
             toray_paths(15,3,3)=   toray_paths(2,3,3)
             root_str(15)       =   root_str(2) 

             !lohan3 data 16 :
             toray_paths(16,1,1)=   '/c/idl/source/autoonetwo/server/bin/toray1.8_noclob_pgf'
             toray_paths(16,2,1)=   toray_paths(16,1,1) ! no change in grid ?
             toray_paths(16,3,1)=   toray_paths(16,1,1)
             toray_paths(16,1,2)=   toray_paths(2,1,2)
             toray_paths(16,2,2)=   toray_paths(2,2,2)
             toray_paths(16,3,2)=   toray_paths(2,3,2)
             toray_paths(16,1,3)=   toray_paths(2,1,3)
             toray_paths(16,2,3)=   toray_paths(2,2,3)
             toray_paths(16,3,3)=   toray_paths(2,3,3)
             root_str(16)       =   root_str(2) 


             !ranier data 17 :
             toray_paths(17,1,1)=   toray_paths(2,1,1)
             toray_paths(17,2,1)=   toray_paths(2,2,1)
             toray_paths(17,3,1)=   toray_paths(2,3,1)
             toray_paths(17,1,2)=   toray_paths(2,1,2)
             toray_paths(17,2,2)=   toray_paths(2,2,2)
             toray_paths(17,3,2)=   toray_paths(2,3,2)
             toray_paths(17,1,3)=   toray_paths(2,1,3)
             toray_paths(17,2,3)=   toray_paths(2,2,3)
             toray_paths(17,3,3)=   toray_paths(2,3,3)
             root_str(17)       =   root_str(2) 


             !zeus  data 18 :
             toray_paths(18,1,1)=   toray_paths(2,1,1)
             toray_paths(18,2,1)=   toray_paths(2,2,1)
             toray_paths(18,3,1)=   toray_paths(2,3,1)
             toray_paths(18,1,2)=   toray_paths(2,1,2)
             toray_paths(18,2,2)=   toray_paths(2,2,2)
             toray_paths(18,3,2)=   toray_paths(2,3,2)
             toray_paths(18,1,3)=   toray_paths(2,1,3)
             toray_paths(18,2,3)=   toray_paths(2,2,3)
             toray_paths(18,3,3)=   toray_paths(2,3,3)
             root_str(18)       =   root_str(2) 


             !uscwsd data 19 :
             toray_paths(19,1,1)=   toray_paths(4,1,1)
             toray_paths(19,2,1)=   toray_paths(4,2,1)
             toray_paths(19,3,1)=   toray_paths(4,3,1)
             toray_paths(19,1,2)=   toray_paths(4,1,2)
             toray_paths(19,2,2)=   toray_paths(4,2,2)
             toray_paths(19,3,2)=   toray_paths(4,3,2)
             toray_paths(19,1,3)=   toray_paths(4,1,3)
             toray_paths(19,2,3)=   toray_paths(4,2,3)
             toray_paths(19,3,3)=   toray_paths(4,3,3)
             root_str(19)       =   root_str(4) 


             !EOS data 20 :
             toray_paths(20,1,1)=   toray_paths(2,1,1)
             toray_paths(20,2,1)=   toray_paths(2,2,1)
             toray_paths(20,3,1)=   toray_paths(2,3,1)
             toray_paths(20,1,2)=   toray_paths(2,1,2)
             toray_paths(20,2,2)=   toray_paths(2,2,2)
             toray_paths(20,3,2)=   toray_paths(2,3,2)
             toray_paths(20,1,3)=   toray_paths(2,1,3)
             toray_paths(20,2,3)=   toray_paths(2,2,3)
             toray_paths(20,3,3)=   toray_paths(2,3,3)
             root_str(20)       =   root_str(2)
  
             !PHOEBE data 21 :
             toray_paths(21,1,1)=   toray_paths(2,1,1)
             toray_paths(21,2,1)=   toray_paths(2,2,1)
             toray_paths(21,3,1)=   toray_paths(2,3,1)
             toray_paths(21,1,2)=   toray_paths(2,1,2)
             toray_paths(21,2,2)=   toray_paths(2,2,2)
             toray_paths(21,3,2)=   toray_paths(2,3,2)
             toray_paths(21,1,3)=   toray_paths(2,1,3)
             toray_paths(21,2,3)=   toray_paths(2,2,3)
             toray_paths(21,3,3)=   toray_paths(2,3,3)
             root_str(21)       =   root_str(2) 

             !HESTIA data 22 :
             toray_paths(22,1,1)=   toray_paths(2,1,1)
             toray_paths(22,2,1)=   toray_paths(2,2,1)
             toray_paths(22,3,1)=   toray_paths(2,3,1)
             toray_paths(22,1,2)=   toray_paths(2,1,2)
             toray_paths(22,2,2)=   toray_paths(2,2,2)
             toray_paths(22,3,2)=   toray_paths(2,3,2)
             toray_paths(22,1,3)=   toray_paths(2,1,3)
             toray_paths(22,2,3)=   toray_paths(2,2,3)
             toray_paths(22,3,3)=   toray_paths(2,3,3)
             root_str(22)       =   root_str(2) 

             !ELCAP data 23:
             toray_paths(23,1,1)=   toray_paths(2,1,1)
             toray_paths(23,2,1)=   toray_paths(2,2,1)
             toray_paths(23,3,1)=   toray_paths(2,3,1)
             toray_paths(23,1,2)=   toray_paths(2,1,2)
             toray_paths(23,2,2)=   toray_paths(2,2,2)
             toray_paths(23,3,2)=   toray_paths(2,3,2)
             toray_paths(23,1,3)=   toray_paths(2,1,3)
             toray_paths(23,2,3)=   toray_paths(2,2,3)
             toray_paths(23,3,3)=   toray_paths(2,3,3)
             root_str(23)       =   root_str(2) 


             !AJAX data 24 :
             toray_paths(24,1,1)=   toray_paths(2,1,1)
             toray_paths(24,2,1)=   toray_paths(2,2,1)
             toray_paths(24,3,1)=   toray_paths(2,3,1)
             toray_paths(24,1,2)=   toray_paths(2,1,2)
             toray_paths(24,2,2)=   toray_paths(2,2,2)
             toray_paths(24,3,2)=   toray_paths(2,3,2)
             toray_paths(24,1,3)=   toray_paths(2,1,3)
             toray_paths(24,2,3)=   toray_paths(2,2,3)
             toray_paths(24,3,3)=   toray_paths(2,3,3)
             root_str(24)       =   root_str(2) 


             !BOB5 data 25 :

             toray_paths(25,1,1)= 'not yet defined'
             toray_paths(25,2,1)= 'not yet defined'
             toray_paths(25,3,1)= 'not yet defined'
             toray_paths(25,1,2)= 'not yet defined'
             toray_paths(25,2,2)= 'not yet defined'
             toray_paths(25,3,2)= 'not yet defined'
             toray_paths(25,1,3)= 'not yet defined'
             toray_paths(25,2,3)= 'not yet defined'
             toray_paths(25,3,3)= 'not yet defined' 
             root_str(25) = '/p/linux/toray' 


             !PLUTO data 26 :

             toray_paths(26,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(26,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(26,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(26,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(26,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(26,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(26,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(26,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(26,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(26) = '/p/linux/toray' 


             !STAR12 data: LOHAN4 27:

             toray_paths(27,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(27,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(27,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(27,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(27,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(27,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(27,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(27,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(27,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(27) = '/p/linux/toray' 

             !lohan4 data: (32 bit toray) NFRC 28 :

             toray_paths(28,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(28,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(28,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(28,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(28,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(28,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(28,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(28,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(28,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(28) = '/p/linux/toray' 


             !STAR1 data:

             toray_paths(35,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(35,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(35,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(35,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(35,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(35,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(35,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(35,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(35,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(35) = '/p/linux/toray' 

             !STAR2 data:

             toray_paths(36,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(36,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(36,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(36,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(36,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(36,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(36,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(36,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(36,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(36) = '/p/linux/toray'

 
             !STAR3 data:

             toray_paths(37,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(37,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(37,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(37,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(37,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(37,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(37,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(37,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(37,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(37) = '/p/linux/toray' 

             !STAR4 data:

             toray_paths(38,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(38,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(38,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(38,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(38,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(38,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(38,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(38,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(38,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(38) = '/p/linux/toray' 

             !STAR5 data:

             toray_paths(39,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(39,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(39,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(39,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(39,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(39,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(39,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(39,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(39,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(39) = '/p/linux/toray' 

             !STAR6 data:

             toray_paths(40,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(40,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(40,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(40,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(40,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(40,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(40,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(40,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(40,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(40) = '/p/linux/toray' 


             !BENTEN  data:

             toray_paths(41,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(41,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(41,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(41,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(41,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(41,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(41,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(41,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(41,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(41) = '/p/linux/toray' 

             !LOHAN5  data:

             toray_paths(42,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(42,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(42,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(42,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(42,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(42,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(42,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(42,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(42,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(42) = '/p/linux/toray' 


             !LOHAN6  data:

             toray_paths(43,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(43,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(43,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(43,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(43,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(43,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(43,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(43,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(43,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(43) = '/p/linux/toray' 


             !LOHAN7  data:

             toray_paths(44,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(44,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(44,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(44,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(44,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(44,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(44,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(44,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(44,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(44) = '/p/linux/toray' 

             !LOAHN8  data:

             toray_paths(45,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(45,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(45,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(45,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(45,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(45,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(45,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(45,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(45,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(45) = '/p/linux/toray' 

             !LOAHN9  data:

             toray_paths(46,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(46,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(46,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(46,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(46,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(46,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(46,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(46,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(46,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(46) = '/p/linux/toray' 

            !LOAHN10  data:

             toray_paths(47,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(47,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(47,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(47,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(47,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(47,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(47,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(47,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(47,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(47) = '/p/linux/toray' 


            !LOAHN11  data:

             toray_paths(48,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(48,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(48,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(48,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(48,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(48,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(48,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(48,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(48,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(48) = '/p/linux/toray' 

            !LOAHN12  data:

             toray_paths(49,1,1)= '/p/linux/toray/v1.60/129_51/toray'
             toray_paths(49,2,1)= '/p/linux/toray/v0.97/129_51/toray'
             toray_paths(49,3,1)= '/p/linux/toray/v1.70/129_51/toray' 
             toray_paths(49,1,2)= '/p/linux/toray/v1.60/65_51/toray'
             toray_paths(49,2,2)= '/p/linux/toray/v0.97/65_51/toray'
             toray_paths(49,3,2)= '/p/linux/toray/v1.70/65_51/toray' 
             toray_paths(49,1,3)= '/p/linux/toray/v1.60/129_201/toray'
             toray_paths(49,2,3)= '/p/linux/toray/v0.97/129_201/toray' 
             toray_paths(49,3,3)= '/p/linux/toray/v1.70/129_201/toray'   
             root_str(49) = '/p/linux/toray' 



             !----------------------------------HSJ 2-14-2013--------------------------
             !toray_paths(ncpu_arch,n_versions,grid_types) set all to same executable:
             !----------------------------------HSJ 2-14-2013--------------------------
             toray_paths(:,:,:) = '/task/imd/toray/xtoray'
             root_str(:) = '/task/imd/toray'


         END SUBROUTINE set_rf_data



        SUBROUTINE get_preplt(kj,ncrt,nout,preplt_path_out,len_str, &
                              preplt_py)
!-------------------------------------------------------------------------
!          subroutine returns a fully qualified executable name for
!          preplt in preplt_path_out
!  INPUT
!  kj            radial grid size (these are used to encode path)
!  ncrt,nout     unit nos. for screen and outone files
!       NOTE nout is disconnected when this routine is called


!  OUTPUT
!  preplt_path_out a fully qualified executable name for preplt
!                 with no leading blanks
!  preplt_py      a fully qualified executable name for  python
!                 program "file_attrib.py"
!  len_str        lengh of siginificant part of preplt_path_out
!
!-------------------------------------------------------HSJ--1/15/03------
           USE io,                   ONLY : versid

           IMPLICIT NONE
           INTEGER, INTENT(in) :: kj,ncrt,nout
           INTEGER, INTENT(out) :: len_str
           CHARACTER(len=*),INTENT(out) :: preplt_path_out
           CHARACTER(len=*),INTENT(out) :: preplt_py
           INTEGER i1,i2,j,tgrid
           CHARACTER charkj*3,vid*5
!           character *12 storesp
           LOGICAL BACK


             PRINT *,'get_preplt,plotcoedid = ',plotcodeid 

!              taurus data:
!            for easy reading/changing do it this way rather than compound
!            root_str to get preplt_paths
             root_str_prplt(1) = '/usr/local/bin/preplt/'
             preplt_paths(1,1) = '/usr/local/bin/preplt/51/preplt' 
             preplt_paths(1,2) = '/usr/local/bin/preplt/201/preplt' 

             !venus
             preplt_paths(1,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(1,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(1,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(1)  =   '/p/linux/preplt/' 

             !delphi2 data:
             preplt_paths(2,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(2,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(2,3)  =   '/p/linux/preplt/101/preplt' 
             root_str_prplt(2)  =   '/p/linux/preplt/' 

 
             !lohan2 data:
             preplt_paths(3,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(3,2)  =   '/p/linux/preplt/201/preplt' 
             root_str_prplt(3)  =   '/p/linux/preplt/' 

             !hydra data:
             preplt_paths(4,1) = '/u/stjohn/preplt/hp/51/preplt'
             preplt_paths(4,2) = '/u/stjohn/preplt/hp/201/preplt' 
             root_str_prplt(4) = '/u/stjohn/preplt/hp/' 

             !delphi data:
             preplt_paths(5,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(5,2)  =   '/p/linux/preplt/201/preplt' 
             preplt_paths(5,3)  =   '/p/linux/preplt/101/preplt' 
             root_str_prplt(5)  =   '/p/linux/preplt/' 


             !cardea data:
             preplt_paths(6,1) = '/u/stjohn/preplt/osf/51/preplt'
             preplt_paths(6,2) = '/u/stjohn/preplt/osf/201/preplt'  
             root_str_prplt(6) = '/u/stjohn/preplt/osf/'  

             !katze data:
             preplt_paths(7,1) = '/u/stjohn/preplt/osf/51/preplt'
             preplt_paths(7,2) = '/u/stjohn/preplt/osf/201/preplt'  
             root_str_prplt(7) = '/u/stjohn/preplt/osf/'  

             !nemsis data:
             preplt_paths(8,1) = '/u/stjohn/preplt/osf/51/preplt'
             preplt_paths(8,2) = '/u/stjohn/preplt/osf/201/preplt'  
             root_str_prplt(8) = '/u/stjohn/preplt/osf/'
  
             !ulam data:
             preplt_paths(9,1) = '/u/stjohn/preplt/osf/51/preplt'
             preplt_paths(9,2) = '/u/stjohn/preplt/osf/201/preplt'  
             root_str_prplt(9) = '/u/stjohn/preplt/osf/'  

             !irenic data:
             preplt_paths(10,1) = '/u/stjohn/preplt/osf/51/preplt'
             preplt_paths(10,2) = '/u/stjohn/preplt/osf/201/preplt'  
             root_str_prplt(10) = '/u/stjohn/preplt/osf/'  

             !hera  data:
             preplt_paths(11,1) = '/u/stjohn/preplt/osf/51/preplt'
             preplt_paths(11,2) = '/u/stjohn/preplt/osf/201/preplt'  
             root_str_prplt(11) = '/u/stjohn/preplt/osf/'  

             !phobos data:
             preplt_paths(12,1) = '/u/stjohn/preplt/osf/51/preplt'
             preplt_paths(12,2) = '/u/stjohn/preplt/osf/201/preplt'  
             root_str_prplt(12) = '/u/stjohn/preplt/osf/'  

             !uscws8 data:
             preplt_paths(13,1) = '/u/stjohn/preplt/osf/51/preplt'
             preplt_paths(13,2) = '/u/stjohn/preplt/osf/201/preplt'  
             root_str_prplt(13) = '/u/stjohn/preplt/osf/'  

             !hermes data:
             preplt_paths(14,1) = '/u/stjohn/preplt/osf/51/preplt'
             preplt_paths(14,2) = '/u/stjohn/preplt/osf/201/preplt'  
             preplt_paths(14,3) = '/u/stjohn/preplt/osf/101/preplt'  
             root_str_prplt(14) = '/u/stjohn/preplt/osf/'  
 
             !AETNA data:
             preplt_paths(15,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(15,2)  =   '/p/linux/preplt/201/preplt' 
             root_str_prplt(15)  =   '/p/linux/preplt/' 

             !lohan3 data:
             preplt_paths(16,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(16,2)  =   '/p/linux/preplt/201/preplt' 
             preplt_paths(16,3)  =   '/p/linux/preplt/101/preplt' 
             root_str_prplt(16)  =   '/p/linux/preplt/' 

 
             !RANIER data:
             preplt_paths(17,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(17,2)  =   '/p/linux/preplt/201/preplt' 
             root_str_prplt(17)  =   '/p/linux/preplt/' 
 
             !zeus data:
             preplt_paths(18,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(18,2)  =   '/p/linux/preplt/201/preplt' 
             preplt_paths(18,3)  =   '/p/linux/preplt/101/preplt' 
             root_str_prplt(18)  =   '/p/linux/preplt/' 


             !uscwsd data:
             preplt_paths(19,1) = '/u/stjohn/preplt/hp/51/preplt'
             preplt_paths(19,2) = '/u/stjohn/preplt/hp/201/preplt' 
             root_str_prplt(19) = '/u/stjohn/preplt/hp/' 


             !EOS data:
             preplt_paths(20,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(20,2)  =   '/p/linux/preplt/201/preplt' 
             preplt_paths(20,3)  =   '/p/linux/preplt/101/preplt' 
             root_str_prplt(20)  =   '/p/linux/preplt/' 


             !PHOEBE data:
             preplt_paths(21,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(21,2)  =   '/p/linux/preplt/201/preplt' 
             preplt_paths(21,3)  =   '/p/linux/preplt/101/preplt' 
             root_str_prplt(21)  =   '/p/linux/preplt/' 

             !HESTIA data:
             preplt_paths(22,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(22,2)  =   '/p/linux/preplt/201/preplt' 
             preplt_paths(22,3)  =   '/p/linux/preplt/101/preplt' 
             root_str_prplt(22)  =   '/p/linux/preplt/' 


             !ELCAP data:
             preplt_paths(23,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(23,2)  =   '/p/linux/preplt/201/preplt' 
             root_str_prplt(23)  =   '/p/linux/preplt/' 


             !AJAX data:
             preplt_paths(24,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(24,2)  =   '/p/linux/preplt/201/preplt' 
             root_str_prplt(24)  =   '/p/linux/preplt/' 

             !BOB5 data:
             preplt_paths(25,1)  =   'not yet defined' 
             preplt_paths(25,2)  =   'not yet defined'    
             root_str_prplt(25)  =   'not yet defined' 


             !pluto  data:
             preplt_paths(26,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(26,2)  =   '/p/linux/preplt/201/preplt' 
             root_str_prplt(26)  =   '/p/linux/preplt/' 




             !lohan4:
             preplt_paths(27,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(27,2)  =   '/p/linux/preplt/201/preplt' 
             preplt_paths(27,3)  =   '/p/linux/preplt/101/preplt' 
             root_str_prplt(27)  =   '/p/linux/preplt/' 



             !NFRC data: 
             preplt_paths(28,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(28,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(28,3)  =   '/p/linux/preplt/101/preplt' 
             root_str_prplt(28)  =   '/p/linux/preplt/' 

             !JAGUAR data:
             preplt_paths(29,1)  =   '/p/linux/preplt/51/prepltv4.0_51'  !jmp
             preplt_paths(29,2)  =   '/p/linux/preplt/201/preplt' !jmp 
             root_str_prplt(29)  =   '/p/linux/preplt/'           !jmp

             !RFPLASMA data:
             preplt_paths(30,1)  =   '/p/linux/preplt/51/prepltv4.0_51'  !jmp
             preplt_paths(30,2)  =   '/p/linux/preplt/201/preplt' !jmp 
             root_str_prplt(30)  =   '/p/linux/preplt/'           !jmp

             !SNUIBM data:
             preplt_paths(31,1)  =   '/p/linux/preplt/51/prepltv4.0_51'  !jmp
             preplt_paths(31,2)  =   '/p/linux/preplt/201/preplt' !jmp 
             root_str_prplt(31)  =   '/p/linux/preplt/'           !jmp

             !ISIS1:
             preplt_paths(32,1)  =   '/p/linux/preplt/51/prepltv4.0_51'  !jmp
             preplt_paths(32,2)  =   '/p/linux/preplt/201/preplt' !jmp 
             root_str_prplt(32)  =   '/p/linux/preplt/'           !jmp

             !isis2 data: (32 bit preplt) jmp.ornl
             preplt_paths(33,1)  =   '/data1/onetwo/preplt/51/prepltv4.0_51'  !jmp
             preplt_paths(33,2)  =   '/data1/onetwo/preplt/201/preplt' !jmp 
             root_str_prplt(33)  =   '/data1/onetwo/preplt/'           !jmp

             !BENTEN data: (32 bit preplt) jmp.ornl
             preplt_paths(34,1)  =   '/data1/onetwo/preplt/51/prepltv4.0_51'  !jmp 
             preplt_paths(34,2)  =   '/data1/onetwo/preplt/201/preplt' !jmp 
             root_str_prplt(34)  =   '/data1/onetwo/preplt/'           !jmp


             !LOHAN5
             preplt_paths(35,1)  =    '/p/linux/preplt/51/preplt'           
             preplt_paths(35,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(35,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(35)  =   '/p/linux/preplt/' 

             !HEAD 
             preplt_paths(36,1)  =    '/p/linux/preplt/51/preplt'           
             preplt_paths(36,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(36,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(36)  =   '/p/linux/preplt/' 
             !NODE01 
             preplt_paths(37,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(37,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(37,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(37)  =   '/p/linux/preplt/' 

             !NODE02 
             preplt_paths(38,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(38,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(38,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(38)  =   '/p/linux/preplt/' 

             !!NODE03  
             preplt_paths(39,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(39,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(39,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(39)  =   '/p/linux/preplt/' 

             !NODE04 
             preplt_paths(40,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(40,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(40,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(40)  =   '/p/linux/preplt/' 

             !NODE05 
             preplt_paths(41,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(41,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(41,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(41)  =   '/p/linux/preplt/' 

             !NODE06
             preplt_paths(42,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(42,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(42,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(42)  =   '/p/linux/preplt/'
 
             !LOHAN1
             preplt_paths(43,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(43,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(43,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(43)  =   '/p/linux/preplt/' 

             !STAR
             preplt_paths(44,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(44,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(44,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(44)  =   '/p/linux/preplt/' 
             !STAR1
             preplt_paths(45,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(45,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(45,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(45)  =   '/p/linux/preplt/' 
             !STAR2
             preplt_paths(46,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(46,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(46,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(46)  =   '/p/linux/preplt/' 
             !STAR3
             preplt_paths(47,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(47,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(47,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(47)  =   '/p/linux/preplt/' 
             !STAR4
             preplt_paths(48,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(48,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(48,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(48)  =   '/p/linux/preplt/' 
             !STAR5
             preplt_paths(49,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(49,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(49,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(49)  =   '/p/linux/preplt/' 
             !STAR6
             preplt_paths(50,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(50,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(50,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(50)  =   '/p/linux/preplt/' 
             !STAR7
             preplt_paths(51,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(51,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(51,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(51)  =   '/p/linux/preplt/' 
             !STAR8
             preplt_paths(52,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(52,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(52,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(52)  =   '/p/linux/preplt/' 
             !STAR9
             preplt_paths(53,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(53,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(53,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(53)  =   '/p/linux/preplt/' 
             !STAR10
             preplt_paths(54,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(54,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(54,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(54)  =   '/p/linux/preplt/' 
             !STAR11
             preplt_paths(55,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(55,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(55,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(55)  =   '/p/linux/preplt/' 
             !STAR12
             preplt_paths(56,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(56,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(56,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(56)  =   '/p/linux/preplt/' 

             !VIDAR
             preplt_paths(1,1)  =   '/p/linux/preplt/51/preplt' 
             preplt_paths(1,2)  =   '/p/linux/preplt/201/preplt'
             preplt_paths(1,3)  =   '/p/linux/preplt/101/preplt'
             root_str_prplt(1)  =   '/p/linux/preplt/' 

             CALL get_hostname(ncrt,nout)
        

!          decode version id:
           vid ='0.00'
           i1 =INDEX(versid,'V')
           IF(i1 == 0)THEN
               i1 =INDEX(versid,'v')
           ENDIF
           IF(i1 == 0)THEN
              WRITE(ncrt,611)
611           FORMAT(2x,'error in decoding version number in sub get_prplt',/, &
                   2x,'you will have to run preplt manually')
           ENDIF
           IF(i1 .ne. 0)THEN
              i2  = LEN_TRIM(versid)
              vid = versid(i1:i2)
           ENDIF

!          only other piece of informatin required for preplt is
!          kj, the transport compiled grid size
           DO j=1,grid_types_prplt
              i1=LEN_TRIM(ADJUSTL(root_str_prplt(host_index)))
              i2 =INDEX(preplt_paths(host_index ,j),'preplt',BACK =.TRUE.)
              i1=i1+1
              i2 =i2-2
              IF(i2 .GT. i1)THEN
                 charkj =preplt_paths(host_index ,j)(i1:i2)
!                 print *,'charkj, get_preplt =',charkj
                 WRITE(storesp,FMT='(a)')charkj(1:LEN_TRIM(charkj))
                 READ(storesp,FMT='(i3)')tgrid
!                 print *,'tgrid =',tgrid
                 IF(tgrid .EQ. kj)THEN
                    len_str = LEN_TRIM(preplt_paths(host_index ,j))
                    IF(is_32_bit)THEN
                       i1=i1-1 ! now i1 is  at the / after preplt
                       preplt_path_out=preplt_paths(host_index ,j)(1:i1)//'32bit/'//preplt_paths(host_index ,j)(i1+1:len_str)//vid 
                    ELSE
                       preplt_path_out =  preplt_paths(host_index ,j)(1:len_str)//vid 
                    ENDIF
                    go to 400
                 ENDIF
              ENDIF 
           ENDDO
           WRITE(ncrt,610)
610        FORMAT(2x,'error in decoding charkj in sub get_preplt',/, &
           2x,'you will have to run preplt manually')
           WRITE(ncrt,FMT='(" tgrid,kj = ",2x,i5,2x,i5)')tgrid,kj
           RETURN

 400       len_str =LEN_TRIM(preplt_path_out)


 
           WRITE(ncrt,600)preplt_path_out(1:len_str)
!           write(nout,600)preplt_path_out(1:len_str)
 600       FORMAT(2x,'sub get_preplt has set the preplt  executeable ',/, &
                 2x,'to :  ',a)
!          get python file path:
!           preplt_py = preplt_path_out(1:len_str-6)//'file_attrib.py' 
!           above line changed HSJ 6/14/07 to accomodate more general
!          preplt file names: 
           preplt_py = preplt_path_out(1:SCAN(preplt_path_out,'/', &
                BACK = .TRUE.))//'file_attrib.py'
        END SUBROUTINE get_preplt





        SUBROUTINE get_pedestal(ncrt,nout,pedestal_path_out,len_str)
!-------------------------------------------------------------------------
!          subroutine returns a fully qualified executable name for
!          PEDESTAL  in pedestal_path_out
!  INPUT
!  ncrt,nout     unit nos. for screen and outone files
!       NOTE nout is disconnected when this routine is called
!  pedestal_path  user supplied path (from inone )

!  OUTPUT
!  pedestal_path_out a fully qualified executable name for pedestal
!                 with no leading blanks
!  len_str        lengh of siginificant part of pedestal_path_out
!
!-------------------------------------------------------HSJ--6/25/04------

           IMPLICIT NONE
           INTEGER, INTENT(in) :: ncrt,nout
           INTEGER, INTENT(out) :: len_str
           CHARACTER(len=*),INTENT(out) :: pedestal_path_out
           INTEGER i1,i2,j,tgrid
           CHARACTER charkj*3
!           character *12 storesp


             !venus data :  
             pedestal_paths(1)    =  '/p/linux/pedestal/pedestal' 

             !delphi2 data:
             pedestal_paths(2)    =  '/p/linux/pedestal/pedestal' 

             !lohan2 data:
             pedestal_paths(3)   =    pedestal_paths(2) 

             !hydra data:
             pedestal_paths(4)   =  'none'

             !delphi data:
             pedestal_paths(5)   =   pedestal_paths(2) 

             !cardea data:
             pedestal_paths(6)   = 'none'

             !katze data:
             pedestal_paths(7)   = 'none'

             !nemsis data:
             pedestal_paths(8)   = 'none'

             !ulam data:
             pedestal_paths(9)   = 'none'

             !irenic data:
             pedestal_paths(10)   = 'none'

             !hera  data:
             pedestal_paths(11)   = 'none'

             !phobos data:
             pedestal_paths(12)   = 'none'

             !uscws8 data:
             pedestal_paths(13)   = 'none'
 
             !hermes data:
             pedestal_paths(14)   = 'none'

             !AETNA data:
             pedestal_paths(15)   =    pedestal_paths(2) 

             !lohan3 data:
             pedestal_paths(16)   =    pedestal_paths(2) 

             !RANIER data:
             pedestal_paths(17)   =    pedestal_paths(2)

             !zeus data:
             pedestal_paths(18)   =    pedestal_paths(2)

             !USCWSD  data:
             pedestal_paths(19)   =    pedestal_paths(2)

             !EOS data:
             pedestal_paths(20)   =    pedestal_paths(2) 

             !PHOEBE data:
             pedestal_paths(21)   =    pedestal_paths(2) 

             !HESTIA data:
             pedestal_paths(22)   =    pedestal_paths(2) 


             !ELCAP data:
             pedestal_paths(23)   =    pedestal_paths(2) 

             !AJAX  data:
             pedestal_paths(24)   =    pedestal_paths(2) 

             !BOB5  data:
             pedestal_paths(25)    =  'not yet defined' 

             !PLUTO data:
             pedestal_paths(26)    =  '/p/linux/pedestal/pedestal' 

             !star12 data:
             pedestal_paths(27)    =  '/p/linux/pedestal/pedestal' 
      
             !lohan4 data:
             pedestal_paths(28)    =  pedestal_paths(2)

             !nfrc data:
             pedestal_paths(29)    =  pedestal_paths(2)

             !rfplasma data: 
             pedestal_paths(30)    =  pedestal_paths(2)

             !JAGUAR
             pedestal_paths(31)    =  pedestal_paths(2)

             !SNUIBM
             pedestal_paths(32)    =  pedestal_paths(2)

             !ISIS1-2
             pedestal_paths(33)    =  pedestal_paths(2)
             pedestal_paths(34)    =  pedestal_paths(2)

             ! STAR!-6         
             pedestal_paths(35)    =  pedestal_paths(2)
             pedestal_paths(36)    =  pedestal_paths(2)         
             pedestal_paths(37)    =  pedestal_paths(2)         
             pedestal_paths(38)    =  pedestal_paths(2)         
             pedestal_paths(39)    =  pedestal_paths(2)         
             pedestal_paths(40)    =  pedestal_paths(2)  
 
             !BENTEN      
             pedestal_paths(41)    =  pedestal_paths(2)    
 
             !lohan5      
             pedestal_paths(42)    =  pedestal_paths(2)    
 
             !loahn6     
             pedestal_paths(43)    =  pedestal_paths(2)    
 
             !loahn7     
             pedestal_paths(44)    =  pedestal_paths(2)    

             !loahn8     
             pedestal_paths(45)    =  pedestal_paths(2)    

             !loahn9     
             pedestal_paths(46)    =  pedestal_paths(2) 
   
             !loahn10     
             pedestal_paths(47)    =  pedestal_paths(2) 
   
             !loahn11     
             pedestal_paths(48)    =  pedestal_paths(2)
    
             !loahn12     
             pedestal_paths(49)    =  pedestal_paths(2)    


             ! HSJ 2/14/2013
             pedestal_paths(:) = 'none'


             IF(LEN_TRIM(pedestal_path) == 0)then
                CALL get_hostname(ncrt,nout)

                pedestal_path_out = pedestal_paths(host_index)
             ELSE
                pedestal_path_out = pedestal_path
             ENDIF

             len_str =LEN_TRIM(pedestal_path_out)


        END SUBROUTINE get_pedestal






        SUBROUTINE get_curray(ncrt,nout,curray_path_out,len_str)
! --------------------------------------------------------------------------------
!       get fully qualified name of curray executable on this host
! --------------------------------------------------------------------------------
      USE mhdpar ,ONLY :  nw,nh
      IMPLICIT NONE
           INTEGER, INTENT(out) :: len_str
           INTEGER, INTENT(in) :: ncrt,nout
           CHARACTER(len=*),INTENT(out) :: curray_path_out
           INTEGER j,i2,grid_index


           len_str = LEN_TRIM(curray_path)
           IF(len_str .EQ. 0)THEN  !no user specified path, get default path on this machine:
              ! find what machine Onetwo is running on:
              nchr = ihost_name (host_name)
              IF (nchr .LE. 0) THEN
                host_name = 'indeterminate'
                nchr    =  13
              END IF
              WRITE (ncrt, '(/ a)') ' host name: ' // host_name(1:nchr)
!             write (nout, '(/ a)') ' host name: ' // host_name(1:nchr)
              IF(host_name == 'indeterminate')THEN
                 WRITE (ncrt, 4)
!                write (nout, 4)
 4               FORMAT( ' Error: invalid host found.  Cant lookup ',/, &
                 ' location of preplt  executeable. onetwo terminated')
                 CALL STOP('Subroutine get_curray',0)
              ENDIF
              !ihost gave a valid  host name,get the required associations:
              DO j = 1,ncpu_arch
                 i2 = LEN_TRIM(host_names(j))
                 IF(j == 1)THEN ! venus names are special
                    IF(ADJUSTL(host_name(1:nchr-1)) == host_names(1)(1:5))THEN
                       host_index = 1
                       nchr = 5
                       go to 300
                    ENDIF
                 ELSE
                    IF(host_name(1:nchr) == host_names(j)(1:i2))THEN
                       host_index = j
                       nchr = i2
                       go to 300
                    ENDIF
                 ENDIF
              ENDDO 
              WRITE(ncrt,5)host_name(1:nchr)
!             write(nout,5)host_name(1:nchr)
  5           FORMAT(' The host that Onetwo is running on ',a,/, &
                 ' is not a member  of the Onetwo reckognized '/, &
                 ' hosts for spawning Curray. Onetwo must stop')
              CALL STOP('Subroutine get_curray',0)



 300          CONTINUE  !host_index,host_name(1:nchr) are now known
              PRINT *,'get_curray,host_index =',host_index 

               grid_index = -1
               IF(nw ==  129 .AND. nh == 129)grid_index =1
               IF(nw ==  65 .AND. nh == 65)  grid_index =2
               IF(grid_index <= 0)THEN
                  PRINT *,'Curray not available for nw,nh =',nw,nh
                  CALL STOP("get_curray",1)
               ENDIF
              curray_path_out = curray_paths(host_index,grid_index)

          ELSE ! user provided path :
              curray_path_out = curray_path
          ENDIF
          len_str = LEN_TRIM(curray_path_out)

        END SUBROUTINE get_curray




        SUBROUTINE get_genray(ncrt,nout,genray_path_out,len_str)
! --------------------------------------------------------------------------------
!       get fully qualified name of genray executable on this host
!       **** PRESENTLY, ONLY USER SUPPLIED ****
! --------------------------------------------------------------------------------

      IMPLICIT NONE
           INTEGER, INTENT(out) :: len_str
           INTEGER, INTENT(in) :: ncrt,nout
           CHARACTER(len=*),INTENT(out) :: genray_path_out
           INTEGER j,i2


           write(*,*)'get_genray: write of genray_path',genray_path
           len_str = LEN_TRIM(genray_path)
           IF(len_str .EQ. 0)THEN  !no user specified path, 
                                   !get default path on this machine:
              WRITE(ncrt,6)
!             write(nout,6)
  6           FORMAT(' Presently need to specify genray_ path',/, &
                 '  Onetwo must stop')
              CALL STOP('Subroutine get_genray',0)

              ! find what machine Onetwo is running on:
              nchr = ihost_name (host_name)
              IF (nchr .LE. 0) THEN
                host_name = 'indeterminate'
                nchr    =  13
              END IF
              WRITE (ncrt, '(/ a)') ' host name: ' // host_name(1:nchr)
!             write (nout, '(/ a)') ' host name: ' // host_name(1:nchr)
              IF(host_name == 'indeterminate')THEN
                 WRITE (ncrt, 4)
!                write (nout, 4)
  4              FORMAT( ' Error: invalid host found.  Cant lookup ',/, &
                 ' location of preplt  executable. onetwo terminated')
                 CALL STOP('Subroutine get_genray',0)
              ENDIF
              !ihost gave a valid host name, get the required associations:

              DO j = 1,ncpu_arch
                 i2 = LEN_TRIM(host_names(j))
                 IF(j == 1)THEN ! venus names are special
                    IF(ADJUSTL(host_name(1:nchr-1)) == host_names(1)(1:5))THEN
                       host_index = 1
                       nchr = 5
                       go to 300
                    ENDIF
                 ELSE
                    IF(host_name(1:nchr) == host_names(j)(1:i2))THEN
                       host_index = j
                       nchr = i2
                       go to 300
                    ENDIF
                 ENDIF
              ENDDO 
              WRITE(ncrt,5)host_name(1:nchr)
!             write(nout,5)host_name(1:nchr)
  5           FORMAT(' The host that Onetwo is running on ',a,/, &
                 ' is not a member  of the Onetwo reckognized '/, &
                 ' hosts for spawning GENRAY. Onetwo must stop')
              CALL STOP('Subroutine get_genray',0)



 300          CONTINUE  !host_index,host_name(1:nchr) are now known
              PRINT *,'get_genray,host_index =',host_index 

!              genray_path_out = genray_paths(host_index,grid_index)

          ELSE ! user provided path :
              genray_path_out = genray_path
          ENDIF
          len_str = LEN_TRIM(genray_path_out)


!****************************************************************
!         Temporary section (BobH, Aug 30, 2007) to determine if
!         running on BOB5/DELPHI/LOHAN4/LOHAN%
!         (the only machines checked out).
              ! find what machine Onetwo is running on:
              nchr = ihost_name (host_name)
              IF (nchr .LE. 0) THEN
                host_name = 'indeterminate'
                nchr    =  13
              END IF
              WRITE (ncrt, '(/ a)') ' host name: ' // host_name(1:nchr)
!             write (nout, '(/ a)') ' host name: ' // host_name(1:nchr)
              IF(host_name == 'indeterminate')THEN
                 WRITE (ncrt, 7)
!                write (nout, 7)
  7              FORMAT( ' Error: invalid host found.  Cant lookup ',/, &
                 ' location of preplt  executable. onetwo terminated')
                 CALL STOP('Subroutine get_genray',0)
              ENDIF
              !ihost gave a valid  host name, get the required associations:

              DO j = 1,ncpu_arch
                 i2 = LEN_TRIM(host_names(j))
                 IF(j == 1)THEN ! venus names are special
                    IF(ADJUSTL(host_name(1:nchr-1)) == host_names(1)(1:5))THEN
                       host_index = 1
                       !nchr = 5
                       go to 400
                    ENDIF
                 ELSE
                    IF(host_name(1:nchr) == host_names(j)(1:i2))THEN
                       host_index = j
                       nchr = i2
                       go to 400
                    ENDIF
                 ENDIF
              ENDDO 
              WRITE(ncrt,8)host_name(1:nchr)
!             write(nout,8)host_name(1:nchr)
  8           FORMAT(' The host that Onetwo is running on ',a,/, &
                 ' is not a member  of the Onetwo hosts. '/, &
                 ' Onetwo must stop')
              CALL STOP('Subroutine get_genray',0)


  400         IF(host_name(1:nchr) .ne. 'DELPHI'   .and.      &
                 host_name(1:nchr) .ne. 'LOHAN4'   .and.      &
                 host_name(1:nchr) .ne. 'LOHAN5'   .and.      &
                 host_name(1:nchr) .ne. 'BOB5'     .and.      &
                 host_name(1:nchr) .ne. 'HEAD'     .and.      &
                 host_name(1:nchr-1)  .ne. 'VENUS'   ) THEN
                 WRITE(ncrt,9)host_name(1:nchr)
!                write(nout,9)host_name(1:nchr)
  9              FORMAT('The host that Onetwo is running on ',a,/, &
                 'is not DELPHI/LOHAN4/LOHAN5/BOB5/VENUS, presently only '/, &
                 'hosts for spawning GENRAY. Onetwo must stop')
                 CALL STOP('Subroutine get_genray',0)
              ENDIF
                 
!****************************************************************



        END SUBROUTINE get_genray




        SUBROUTINE get_toq(ncrt,nout,toq_path_out,len_str)
! --------------------------------------------------------------------------------
!       get fully qualified name of toq executable on this host
! --------------------------------------------------------------------------------
      USE mhdpar ,ONLY :  nw,nh
      IMPLICIT NONE
           INTEGER, INTENT(out) :: len_str
           INTEGER, INTENT(in) :: ncrt,nout
           CHARACTER(len=*),INTENT(out) :: toq_path_out
           INTEGER j,i2,grid_index


           len_str = LEN_TRIM(toq_path)
           IF(len_str .EQ. 0)THEN  !no user specified path, get default path on this machine:
              ! find what machine Onetwo is running on:
              nchr = ihost_name (host_name)
              IF (nchr .LE. 0) THEN
                host_name = 'indeterminate'
                nchr    =  13
              END IF
              WRITE (ncrt, '(/ a)') ' host name: ' // host_name(1:nchr)
!             write (nout, '(/ a)') ' host name: ' // host_name(1:nchr)
              IF(host_name == 'indeterminate')THEN
                 WRITE (ncrt, 4)
!                write (nout, 4)
 4               FORMAT( ' Error: invalid host found.  Cant lookup ',/, &
                 ' location of preplt  executeable. onetwo terminated')
                 CALL STOP('Subroutine get_toq',0)
              ENDIF
              !ihost gave a valid  host name,get the required associations:
              DO j = 1,ncpu_arch
                 i2 = LEN_TRIM(host_names(j))
                 IF(j == 1)THEN ! venus names are special
                    IF(ADJUSTL(host_name(1:nchr-1)) == host_names(1)(1:5))THEN
                       host_index = 1
                       nchr = 5
                       go to 300
                    ENDIF
                 ELSE
                    IF(host_name(1:nchr) == host_names(j)(1:i2))THEN
                       host_index = j
                       nchr = i2
                       go to 300
                    ENDIF
                 ENDIF
              ENDDO 

              WRITE(ncrt,5)host_name(1:nchr)
!             write(nout,5)host_name(1:nchr)
  5           FORMAT(' The host that Onetwo is running on ',a,/, &
                 ' is not a member  of the Onetwo reckognized '/, &
                 ' hosts for spawning Toq. Onetwo must stop')
              CALL STOP('Subroutine get_toq',0)



 300          CONTINUE  !host_index,host_name(1:nchr) are now known

               grid_index = -1
               IF(nw ==  129 .AND. nh == 129)grid_index =1
               IF(nw ==  65 .AND. nh == 65)  grid_index =2
               IF(grid_index <= 0)THEN
                  PRINT *,'Toq not available for nw,nh =',nw,nh
                  CALL STOP("get_toq",1)
               ENDIF
              toq_path_out = toq_paths(host_index,grid_index)

          ELSE ! user profided path :
              toq_path_out = toq_path
          ENDIF
          len_str = LEN_TRIM(toq_path_out)
        END SUBROUTINE get_toq



        SUBROUTINE get_toq_base(ncrt,nout,toq_base_out,len_str)
! --------------------------------------------------------------------------------
!       get toq base name on this host to identify locsation of fixb129xy.sc
! --------------------------------------------------------------------------------
      USE mhdpar ,ONLY :  nw,nh
      IMPLICIT NONE
           INTEGER, INTENT(out) :: len_str
           INTEGER, INTENT(in) :: ncrt,nout
           CHARACTER(len=*),INTENT(out) :: toq_base_out
           INTEGER j,i2,grid_index

              ! find what machine Onetwo is running on:
              nchr = ihost_name (host_name)
              IF (nchr .LE. 0) THEN
                host_name = 'indeterminate'
                nchr    =  13
              END IF
              WRITE (ncrt, '(/ a)') ' host name: ' // host_name(1:nchr)
!             write (nout, '(/ a)') ' host name: ' // host_name(1:nchr)
              IF(host_name == 'indeterminate')THEN
                 WRITE (ncrt, 4)
!                write (nout, 4)
 4               FORMAT( ' Error: invalid host found.  Cant lookup ',/, &
                 ' location of preplt  executeable. onetwo terminated')
                 CALL STOP('Subroutine get_toq',0)
              ENDIF
              !ihost gave a valid  host name,get the required associations:

              DO j = 1,ncpu_arch
                 i2 = LEN_TRIM(host_names(j))
                 IF(j == 1)THEN ! venus names are special
                    IF(ADJUSTL(host_name(1:nchr-1)) == host_names(1)(1:5))THEN
                       host_index = 1
                       nchr = 5
                       go to 300
                    ENDIF
                 ELSE
                    IF(host_name(1:nchr) == host_names(j)(1:i2))THEN
                       host_index = j
                       nchr = i2
                       go to 300
                    ENDIF
                 ENDIF
              ENDDO 
              WRITE(ncrt,5)host_name(1:nchr)
!             write(nout,5)host_name(1:nchr)
  5           FORMAT(' The host that Onetwo is running on ',a,/, &
                 ' is not a member  of the Onetwo reckognized '/, &
                 ' hosts for spawning Toq. Onetwo must stop')
              CALL STOP('Subroutine get_toq',0)



 300          CONTINUE  !host_index,host_name(1:nchr) are now known


              toq_base_out = toq_base(host_index)
              len_str = LEN_TRIM(toq_base_out)

        END SUBROUTINE get_toq_base



        SUBROUTINE get_nubeam(ncrt,nout,nubeam_path_out, &
                              nubeam_setup_out,len_str)
! --------------------------------------------------------------------------------
!       get fully qualified name of nubeam  executable on this host
! --------------------------------------------------------------------------------
      IMPLICIT NONE
           INTEGER, INTENT(out) :: len_str
           INTEGER, INTENT(in) :: ncrt,nout
           CHARACTER(len=*),INTENT(out) :: nubeam_path_out,nubeam_setup_out
           INTEGER j,i2



           IF(len_str .EQ. 0)THEN  !no user specified path, get default path on this machine:
              ! find what machine Onetwo is running on:
              nchr = ihost_name (host_name)
              IF (nchr .LE. 0) THEN
                host_name = 'indeterminate'
                nchr    =  13
              END IF
              WRITE (ncrt, '(/ a)') ' host name: ' // host_name(1:nchr)
!             write (nout, '(/ a)') ' host name: ' // host_name(1:nchr)
              IF(host_name == 'indeterminate')THEN
                 WRITE (ncrt, 4)
!                write (nout, 4)
 4               FORMAT( ' Error: invalid host found.  Cant lookup ',/, &
                 ' location of preplt  executeable. onetwo terminated')
                 CALL STOP('Subroutine get_nubeam',0)
              ENDIF
              !ihost gave a valid  host name,get the required associations:

              DO j = 1,ncpu_arch
                 i2 = LEN_TRIM(host_names(j))
                 IF(j == 1)THEN ! venus names are special
                    IF(ADJUSTL(host_name(1:nchr-1)) == host_names(1)(1:5))THEN
                       host_index = 1
                       nchr = 5
                       go to 300
                    ENDIF
                 ELSE
                    IF(host_name(1:nchr) == host_names(j)(1:i2))THEN
                       host_index = j
                       nchr = i2
                       go to 300
                    ENDIF
                 ENDIF
              ENDDO 
              WRITE(ncrt,5)host_name(1:nchr)
!             write(nout,5)host_name(1:nchr)
  5           FORMAT(' The host that Onetwo is running on ',a,/, &
                 ' is not a member  of the Onetwo reckognized '/, &
                 ' hosts for spawning Nubeam. Onetwo must stop')
              CALL STOP('Subroutine get_nubeam',0)



 300          CONTINUE  !host_index,host_name(1:nchr) are now known

              nubeam_path_out = nubeam_paths(host_index)
              nubeam_setup_out =  nubeam_setup(host_index)
          ELSE ! user provided path :
              nubeam_path_out = nubeam_path
              nubeam_setup_out  = nubeam_setup_ext
          ENDIF
          len_str = LEN_TRIM(nubeam_path_out)
        END SUBROUTINE get_nubeam





        SUBROUTINE get_fastcd(ncrt,nout,fcd_path)
!  ----------------------------------------------------------------------------------------
!  ---------- return path to fast cd executable:
!  ----------------------------------------------------------------------------HSJ---------
         IMPLICIT NONE
         INTEGER  strleng,len_str
         INTEGER, INTENT(in)    :: ncrt,nout
         LOGICAL, INTENT(inout) :: fcd_path
         LOGICAL  file_exists
         IF(fcd_path)RETURN
                strleng = LEN_TRIM(fastcd_path)
                IF(strleng .LE. 0) THEN
                   !user did not specify path for fastcd
                   !we will set the path according to what machine
                   !we are running on:

                   CALL get_hostname(ncrt,nout)

                   IF(host_index > 0)THEN
                        fastcd_path = fastcd_paths(host_index,1)
                        fcd_path = .TRUE.
                   ELSE
                     WRITE(nout,1)host_name(1:nchars_host)
                     WRITE(ncrt,1)host_name(1:nchars_host)
                     CALL STOP('cray102.f, no fastcd executeable found',1)
 1                   FORMAT(2x,'ERROR, could not determine path',/,    &
                            '  to fastcd executeable  on this host:', a)
                   ENDIF
                ELSE                !check user path for valid executeable
                   fcd_path = .TRUE.
                   PRINT *,' Fastcd  path specified option'
                   len_str = LEN_TRIM(fastcd_path)
                   file_exists = .FALSE.
                   INQUIRE(FILE = fastcd_path(1:len_str),EXIST =  file_exists)
                  IF(.NOT. file_exists)THEN
                      WRITE(ncrt,2)fastcd_path(1:len_str)
                      WRITE(nout,2)fastcd_path(1:len_str)
  2                   FORMAT(2x,'ERROR specified fastcd not ',/, 2x,'accessible:',/,2x,a)
                      CALL STOP('fastcd_path problem',0)
                  ENDIF
                ENDIF
          print *,'fastcd_path on output ',fastcd_path
          RETURN
        END SUBROUTINE get_fastcd




         SUBROUTINE get_toray(nw,nh,kj,ncrt,nout,    &
                             toray_path_out,len_str)
!-------------------------------------------------------------------------
!          subroutine returns a fully qualified executable name for
!          toray in toray_path_out
!  INPUT
!  nw,nh         mhd grid size 
!  kj            radial grid size (these are used to encode path)
!  ncrt,nout     unit nos. for screen and outone files
!  toray_path    user supplied fully qualified executeable name.
!                this routine checks to see if it is valid. If it is
!                then it will be copied to toray_path_out on output.
!                Must be set to zero length string if it is not used.

!  toray_version 3 digit version number ( eq 0.97, 1.41, etc)  of toray
!                to use, checked only if toray_path  is not set.
!                if toray_version is set then toray path out will 
!                be the executable
!                version that was requested, if such a version exists
!                on the architecture Onetwo is running on. Otherwise
!                a fatal error is signalled. toray_version must be set to 0
!                if a particular version is not requested.

!  OUTPUT
!  toray_path_out a fully qualified executable name for toray
!                 with no leading blanks
!  len_str        lengh of siginificant part of toray_path_out
!  toray_version  version of toray that will be returned in 
!                 (in ext_prog_info)
!
!
! -------------------------------------------------------HSJ--1/15/03------


!           USE ext_prog_info
           IMPLICIT NONE
           INTEGER nw,nh,kj,ncrt,nout
           INTEGER i1,i2,i3,j,i,jsave
           REAL vn 
           INTEGER len_str,len_str2
           CHARACTER *(*) toray_path_out
           CHARACTER char_grid*12, charkj*3,charnw*3,charnh*3,    &
                     charv*8
!           character *12 storesp
           CHARACTER buildstr*256
           CHARACTER *128 buildstr1
           LOGICAL ex
           ex =.FALSE.




           ! first find what machine Onetwo is running on:
           nchr = ihost_name (host_name)
!           print *,'nchr =',nchr
           IF (nchr .LE. 0) THEN
             host_name = 'indeterminate'
             nchr    =  13
           END IF
           WRITE (ncrt, '(/ a)') ' host name: ' // host_name(1:nchr)
!           write (ncrt, '(/ a)') ' host name: ' // host_name
           WRITE (nout, '(/ a)') ' host name: ' // host_name(1:nchr)
           IF(host_name == 'indeterminate')THEN
              WRITE (ncrt, 4)
              WRITE (nout, 4)
 4            FORMAT( ' Error: invalid host found.  Cant lookup ',/, &
             ' location of toray executeable. onetwo terminated')
              CALL STOP('Subroutine get_toray',0)
           ENDIF
           !ihost gave a valid  host name,get the required associations:
           DO j = 1,ncpu_arch
              i2 = LEN_TRIM(host_names(j))
              IF(j == 1)THEN ! venus names are special
                    IF(ADJUSTL(host_name(1:nchr-1)) == host_names(1)(1:5))THEN
                       host_index = 1
                       nchr = 5
                       go to 300
                    ENDIF
              ELSE
                   IF(host_name(1:nchr) == host_names(j)(1:i2))THEN
                      host_index = j
                      go to 300
                   ENDIF
              ENDIF
           ENDDO 
           WRITE(ncrt,5)host_name(1:nchr)
           WRITE(nout,5)host_name(1:nchr)
  5        FORMAT(' The host that Onetwo is running on ',a,/,            &
                 ' is not a member  of the Onetwo reckognized '/,        &
                 ' hosts for spawning Toray. Onetwo must stop')
           CALL STOP('Subroutine get_toray',0)



 300       CONTINUE  !host_index,host_name(1:nchr) are now known

           ! older (pre 5.45 versions of onetwo force specification of input path
           ! v5.45 allows default version to be used:
            IF(LEN_TRIM(toray_path) == 0) THEN
                  toray_path = toray_paths(host_index,1,1)
                  TORAY_VERSION = 1.8
                  WRITE(ncrt,FMT='("Default toray selected by Onetwo : ",a)')LEN_TRIM(toray_path)
                  WRITE(ncrt,FMT='("Default toray version selected by Onetwo : ",f10.3)')toray_version
            ENDIF
    tpath:        IF(LEN_TRIM(toray_path) == 0)THEN



                     write(ncrt,301) ! added 4/07/09 HSJ 
 301                  FORMAT(2x,'ERROR: fully qualified Toray executable name, ',/, &
                       2x,'(eg TORAY_PATH = ),must be given in Onetwo version 5.2 and beyond')
                    CALL EXIT 
               WRITE(storesp,'(i3)')kj
               READ(storesp,'(a)')charkj
               charkj= ADJUSTL(charkj)     !move leading blanks to end
               WRITE(storesp,FMT='(i3)')nw
               READ(storesp,FMT='(a)')charnw
               charnw= ADJUSTL(charnw)     !move leading blanks to end
               IF(nw .NE. nh)THEN
                    WRITE(storesp,'(i2)')nh
                    READ(storesp,FMT='(a)')charnh
!                    print *,'charnw =',charnh
                    charnh = ADJUSTL(charnh)
                    buildstr1 = TRIM(charnw)//'_'//TRIM(charnh)//'_'//                  &
                                       TRIM(charkj)//'/toray'
               ELSE
                   buildstr1 = TRIM(charnw)//'_'//TRIM(charkj)//'/toray'
               ENDIF


              !is there a requested version ?:
              IF( toray_version == 0.0 )THEN
!                 print *,'default case'
                 !the default case, Onetwo will decide what version of
                 !Toray to run. We will search from the most recent
                 !to the oldest version and use the first one that we
                 !find.

                 DO j = n_versions,1,-1
                    buildstr =                                       &
                     ADJUSTL(buildstr1(1:LEN_TRIM(buildstr1)))
                    WRITE(storesp,FMT='(f4.2)')versions(j)
                    READ(storesp,FMT='(a)')charv
                    charv = ADJUSTL(charv)
                    buildstr = 'v'//TRIM(charv)//'/'//TRIM(buildstr)
                    len_str = LEN_TRIM(buildstr)
                    len_str2 = LEN_TRIM(root_str(host_index))
                    buildstr = root_str(host_index)(1:len_str2)//    &
                              '/'//ADJUSTL(buildstr(1:len_str))
                    len_str = LEN_TRIM(buildstr)
!                    print *,'final str =',buildstr(1:len_str)
                    INQUIRE(FILE = buildstr(1:len_str),EXIST  = ex)
                    IF(ex)THEN
                      toray_path_out = buildstr(1:len_str)
                      len_str = LEN_TRIM(toray_path_out)
                      toray_version = versions(j)  !needed in sub ech
                      go to 100
                    ENDIF
                 ENDDO
                    !no acceptable toray found ?? programmer goof!!
                    WRITE(ncrt,7)buildstr(1:len_str),host_name
                    WRITE(nout,7)buildstr(1:len_str),host_name
 7                  FORMAT(2x,'ERROR , could not find a valid Toray ',/,           &
                           2x,'last path tried is :',a,/,                          &
                           2x,'on host :',a)
                    CALL STOP('Subroutine get_toray',0)
 100                CONTINUE
              ELSE             ! user has supplied a preferred version number
                               ! but has not indicated where to find 
                               ! that version of toray. Check
                               ! the default  paths to see if we have an
                               ! appropriate version of toray available.
                 DO j=1,n_versions
                    IF(ABS(toray_version- versions(j)) .LT. vertol)go to 200
                 ENDDO
                 !requested version not available
                 WRITE(ncrt,8)toray_version,host_name(1:nchr)
                 WRITE(nout,8)toray_version,host_name(1:nchr)
 8               FORMAT(2x,'ERROR , requested version of toray ',f6.2, &
                 /,'  on host ',a,' is not available')
                 WRITE(ncrt,9)(versions(j),j=1,n_versions)
                 WRITE(nout,9)(versions(j),j=1,n_versions)
 9               FORMAT(2x,'available versions are ',/,                &
                                              (4(2x,f6.2)))
                 CALL STOP('Subroutine get_toray',0)
  200            CONTINUE
                 !we have an acceptable version number find the path
                 !on the given host by searching for matching version 
                 !and grids in use:
                 DO i = 1,grid_types
                    DO j=1,n_versions
                      i3 =LEN_TRIM(toray_paths(host_index,j,i))
                      i2 = 1
                      !search for '/v'. then next 4 characters form version #
 210                  i1 = INDEX(toray_paths(host_index,j,i)(i2:i3),'/v')
                      IF(i1 .GT. 0)THEN
                        i2=i1+2
                        charv = toray_paths(host_index,j,i)(i2:i2+3)
                        !does charv represent a version #?
                        READ(charv,'(f5.2)',ERR = 220)vn
                        IF(ABS(vn - toray_version) .LT. vertol)THEN
                            !we have a match
                            toray_path_out = ADJUSTL(TRIM(toray_paths(host_index,j,i)))
                            len_str = LEN_TRIM(toray_path_out)
                            jsave = j
                            go to 240
                        ENDIF
                      ELSE               !i1 =0 ,'/v' not found
                         go to 221       !go to next version
                      ENDIF
 220                  IF(i2 .LT. i3)go to 210  !no valid match,search for 
                                               !further occurances of \v
                                               !in the same string
 221                  CONTINUE
                    ENDDO                      !end loop over versions
                    !no version matches (should not happen)
                    WRITE(ncrt,7)toray_version,host_name
                    WRITE(nout,7)toray_version,host_name
                    CALL STOP('sub get_toray program error 1',0)
 240                CONTINUE                      !match found for version
                                                  !what about grid ?

!                    print *,'searching for grid'
                    !search for '/toray from the right towards the left:
                    i3 = LEN_TRIM(toray_paths(host_index,jsave,i))
                    i1 = INDEX(toray_paths(host_index,jsave,i), '/toray',BACK =.TRUE.)
                    i1=i1-1
                    i2 = INDEX(toray_paths(host_index,jsave,i)(1:i1),'/'  ,BACK =.TRUE.)
                    i2 =i2+1
                    !grid is  in positions i2 to i1
                    char_grid = toray_paths(host_index,jsave,i)(i2:i1)
!                    print *,'char_grid =',char_grid
                    !decode char_grid to get nw,nh,kj
                    !search for number of '_' present. There will be 2
                    !if nw .ne. nh and one if nw .eq.nh:
                    i1 = INDEX(char_grid,'_')
                    i2 = INDEX(char_grid,'_',BACK = .TRUE.)
                    IF(i1 .EQ. i2)THEN          !nw = nh
                       !get nw and kj
                       READ(char_grid,'(i3,x,i3)')nw_rf,kj_rf
!                       print *,'nwl,kj1 =',nw_rf,kj_rf
                       IF(nw_rf .EQ. nw .AND. kj_rf .EQ. kj)THEN
                          !success
                          toray_path_out = toray_paths(host_index,  jsave,i)
                          len_str = LEN_TRIM(toray_path_out)
                          go to 260
                       ELSE
                          !grid didnt match check next case
                          go to 290
                       ENDIF
                     ELSE                        !nw .ne. nh
                       READ(char_grid,'(i3,x,i3,x,i3)')nw_rf,nh_rf,kj_rf
!                       print *,'nw1,nh1,kj1 =',nw_rf,nh_rf,kj_rf
                       IF(nw_rf .EQ. nw .AND. kj_rf .EQ. kj   .AND. nh_rf .EQ.  nh)THEN
                          !success
                          toray_path_out = toray_paths(host_index,jsave,i)
                          len_str = LEN_TRIM(toray_path_out)
                          go to 260
                       ELSE
                          !grid didnt match check next case
                          go to 290
                       ENDIF
                     ENDIF
 290                 CONTINUE   
                 ENDDO    !end loop over grid types
                 !checked all versions,grids on machine pointed to by
                 !host_index and did not find a match.
                  WRITE(ncrt,11)toray_version,host_name
                  WRITE(nout,11)toray_version,host_name
 11               FORMAT(2x,'ERROR did not find toray version',f6.3,/,   &
                        2x,'on this machine (',a,')')
                  CALL STOP('Subroutine get_toray',0)
 260             CONTINUE
              ENDIF  !toray version set or not
           ELSE  tpath         ! user has supplied afull path to  executeable
                               ! toray. Check to see if it will work:
              PRINT *,'Toray path specified option'
              len_str = LEN_TRIM(toray_path)
              toray_path_out = toray_path(1:len_str)
              !toray_path_out has to be the fully qualified executable 
              !name of toray:
              INQUIRE(FILE = toray_path_out(1:len_str),EXIST = ex)
              IF(.NOT. ex)THEN
                  WRITE(ncrt,6)toray_path_out(1:len_str)
                  WRITE(nout,6)toray_path_out(1:len_str)
  6               FORMAT(2x,'ERROR specified toray not accessible:',/,  2x,a)
                  CALL STOP('Subroutine get_toray',0)
              ENDIF
              !user has supllied a toray version. We will assume that
              !this version matches the onetwo that is being run:
              nw_rf = nw
              nh_rf = nh
              kj_rf = kj
              !version number must be given for this otpion to work
              IF(toray_version .LE. 0.0)THEN
                  WRITE(ncrt,601)toray_path,toray_version
                  WRITE(nout,601)toray_path,toray_version
 601             FORMAT(2x, 'Error, toray_version must be given',/,   &
                       2x,'if toray_path option is used',/,           &
                       2x,'toray_path =',a,/,                         &
                       2x,'toray_version =',f10.3)
                 CALL STOP('get_toray',0)
              ENDIF
           ENDIF  tpath  !toray path set or not
           WRITE(ncrt,600)toray_path_out(1:len_str)
           WRITE(nout,600)toray_path_out(1:len_str)
 600       FORMAT(2x,'sub get_toray has set the toray executeable ',/, &
                 2x,'to :  ',a)
           
           RETURN
           END SUBROUTINE  get_toray
 

       SUBROUTINE get_xsect_path(ncrt,nout)
! ----------------------------------------------------------------------       
! --- get  location of cross section data files
! ---------------------------------------------------------------HSJ----
!
       IMPLICIT NONE
       INTEGER NCRT,NOUT
       INTEGER, EXTERNAL :: LENGTH

        CALL get_hostname(ncrt,nout)
        SELECT CASE (host_index)
           CASE (1,2,3,5,17,18,20,21,22,23,24,25,26,27,28,35,36,37,38,39,40,41,42,43,44,&
                 45,46,47,48,49,50,51,52,53,56,57)
               onetwo_xsct = "/p/linux/onetwo/"
               print *,'onetwo_xsct=',onetwo_xsct !jmp
           CASE (4)
               onetwo_xsct = "/d/hp"
           CASE (6,7,8,9,10,11,12,13,14,15,16)
               onetwo_xsct = "/d/osf"
           CASE (29) !jmp
               onetwo_xsct = "/home/myung/ppp/onetwo" !jmp
               print *,'onetwo_xsct=',onetwo_xsct !jmp
           CASE (30) !jmp
               onetwo_xsct = "/autofs/spin/home/jbp/ppp/onetwo" !jmp
           CASE (31) !jmp.ibm
               onetwo_xsct = "/home/jbp/onetwo" !jmp.ibm
           CASE (32) !jmp.ibm
               onetwo_xsct = "/home/fusma/user/haha/snu" !jmp.ibm
           CASE (33) !jmp.ornl
               onetwo_xsct = "/data1/onetwo/" !jmp.ornl
           CASE (34) !jmp.ornl
               onetwo_xsct = "/data1/onetwo/" !jmp.ornl
           CASE (54,55)
              WRITE(nout,2)host_name(1:nchars_host)
              WRITE(ncrt,2)host_name(1:nchars_host)
2             FORMAT(2x," ERROR, cross section path /p/linux/onetwo ",/,&
                   "       is not mounted on STAR10,STAR11")
              CALL STOP('get_xsect_path, invalid host name',1)
           CASE default
               WRITE(nout,1)host_name(1:nchars_host)
               WRITE(ncrt,1)host_name(1:nchars_host)
               CALL STOP('get_xsect_path, invalid host name',1)
 1            FORMAT(2x,'ERROR, could not determine path to cross ',/, &
                2x,'section files on this host:', a)
        END SELECT

        nchars_12    =   LENGTH  (onetwo_xsct)

      END SUBROUTINE get_xsect_path







          SUBROUTINE get_hostname(ncrt,nout)
!  -----------------------------------------------------------------------------------------
!  ------ return the host name code is running on ------------------------------------------
!  ------ must match one of the entries in host_names ---------------------------HSJ --------
!-------------------------------------------------------------------------------------------
            IMPLICIT NONE
             
            INTEGER ncrt,nout,j,i2
              nchr = ihost_name (host_name)
              IF (nchr .LE. 0) THEN
                host_name = 'indeterminate'
                nchr    =  13
              END IF
              WRITE (ncrt, '(/ a)') ' host name: ' // host_name(1:nchr)
!             write (nout, '(/ a)') ' host name: ' // host_name(1:nchr)
              IF(host_name == 'indeterminate')THEN
                 WRITE (ncrt, 4)
!                write (nout, 4)
 4               FORMAT( ' Error: invalid host found.  Cant lookup ',/, &
                 ' location of preplt  executeable. onetwo terminated')
                 CALL STOP('Subroutine get_hostname',0)
              ENDIF
              !ihost gave a valid  host name,get the required associations:
              DO j = 1,ncpu_arch
                 i2 = LEN_TRIM(host_names(j))
                 IF(j == 1)THEN ! venus names are special
                    IF(ADJUSTL(host_name(1:nchr-1)) == host_names(1)(1:5))THEN
                       host_index = 1
                       nchr = 5
                       go to 300
                    ENDIF
                 ELSE
! print *,'j,nchar,i2 =',j,nchr,i2
! print *,'host_name(1:nchr) host_names(j)(1:i2) =',host_name(1:nchr),host_names(j)(1:i2) 
                    IF(host_name(1:nchr) == host_names(j)(1:i2))THEN
                       host_index = j
                       nchr = i2
 !print *,'going to 300'
                       go to 300
                    ENDIF
                 ENDIF
! print *,'j,ncpu_arch= ',j,ncpu_arch
              ENDDO 
              WRITE(ncrt,5)host_name(1:nchr)
!             write(nout,5)host_name(1:nchr)
  5           FORMAT(' The host that Onetwo is running on ',a,/, &
                 ' is not reckognized')
              CALL STOP('Subroutine get_hostname',0)



 300          CONTINUE  !host_index,host_name(1:nchr) are now known
              host_name = host_names(host_index)(1:nchr)
              nchars_host = nchr
              DO j=1,ncpu_arch_32
                 if(host_name == host_names_32(j))is_32_bit = .TRUE.
              ENDDO
              RETURN

          END SUBROUTINE get_hostname

          FUNCTION ihost_name(host_name) !jmp.ibm.par start

             IMPLICIT NONE
             INTEGER,EXTERNAL :: IHOST
             INTEGER ihost_name
             CHARACTER *(*) host_name

             IF(host_name_ext .ne. 'undefined') THEN

                ihost_name = LEN_TRIM(host_name_ext)
                host_name = host_name_ext(1:ihost_name)

             ELSE
             
                ihost_name = IHOST(host_name)

             ENDIF

             RETURN

          END FUNCTION !jmp.ibm.par start

       END MODULE ext_prog_info
