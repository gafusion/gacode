   SUBROUTINE write_statefile(rtyp)  
!------------------------------------------------------------------------------
! -- translate to statefile outputs and call shared statefile writing routines
!     rtyp                not used
!     irwflag             read/write flag, 
!                         irwflag = 0  ==>  WRITE the data
!                         irwflag = 1  ==>  READ  the data
!
!     iterdb_file_name    on input holds the name of the state file read in
!                         This subroutine changes the name to the
!                         predefined value in the P_nfreya namelist input
!                         We do it this way so, for example, Onetwo knows
!                         what file to read because onewo specified the file name in
!                         iwhen it creatd the P_nfreya nmaelist file.
!------------------------------------------------------------------------------
       USE nrtype,                                          ONLY : I4B,DP
       
       USE Plasma_properties ,                              ONLY : neut_beam

       USE zonal_data,                                      ONLY :  mf,mfm1

       USE iterdbmd_gcnmp,                                  ONLY : iterdb_file_name,         &
                                                                   irwflag,idb_sufix
                                                                   
       USE io_gcnmp,                                        ONLY : switch_statefile_output,   &
                                                                   statefile_output_name

       USE common_constants,                                ONLY : izero



       IMPLICIT NONE

       INTEGER(I4b) i,out_namel,rtyp
       CHARACTER *4 tail
       LOGICAL txt_set,nc_set


       ! determine if  statefile_output_name has an appropriate suffix:
       txt_set = .FALSE.
       nc_set  = .FALSE.
       out_namel  = LEN_TRIM(statefile_output_name)
       tail = statefile_output_name(out_namel-3:out_namel+1)
       CALL to_upper_case(tail)
       IF( tail == '.TXT')THEN
          txt_set = .TRUE.
       ELSEIF(tail(2:4) == '.NC')THEN
          nc_set = .TRUE.
       ENDIF
           





       iterdb_file_name = statefile_output_name   ! name of output file to create
       irwflag = izero                            ! we will write the output statefile

       IF ( idb_sufix == '.NC' .AND. .NOT. switch_statefile_output )THEN
          ! input and output in netcdf form
          IF(txt_set)THEN        ! strip .txt and add .nc
             iterdb_file_name =    statefile_output_name(1:out_namel-4)//'.nc'
          ELSEIF(nc_set)THEN     ! consistent
             iterdb_file_name =    statefile_output_name
          ELSE                   ! suffix not set so set it:
             iterdb_file_name =    ADJUSTL(statefile_output_name(1:out_namel) // '.nc')
          ENDIF
          CALL iter_dbase_nc
       ELSE IF ( idb_sufix == '.NC' .AND. switch_statefile_output )THEN
          ! input in netcdf form, output in text form
          IF(txt_set)THEN        ! consistent
             iterdb_file_name =    statefile_output_name(1:out_namel)
          ELSEIF(nc_set)THEN     ! strip .nc and add  .txt
             iterdb_file_name =    ADJUSTL( statefile_output_name(1:out_namel-3) // '.txt')
          ELSE                   ! suffix not set so set it:
             iterdb_file_name =    ADJUSTL(statefile_output_name(1:out_namel) // '.txt')
          ENDIF
          CALL iter_dbase_txt
       ELSE IF ( idb_sufix .NE. '.NC' .AND. switch_statefile_output )THEN
          ! input in text  form, output in netcdf  form
          IF(txt_set)THEN        ! strip .txt and and .nc
             iterdb_file_name =    statefile_output_name(1:out_namel-4)//'.nc'
          ELSEIF(nc_set)THEN     ! consistent
             iterdb_file_name =    statefile_output_name
          ELSE                   ! suffix not set so set it:
             iterdb_file_name =    ADJUSTL(statefile_output_name(1:out_namel) // '.nc')
          ENDIF
          CALL iter_dbase_nc
       ELSE
          ! input in text ,output in text form
          IF(nc_set)THEN          ! strip .txt and and .nc
             iterdb_file_name =    statefile_output_name(1:out_namel-3)//'.txt'
          ELSEIF(txt_set)THEN     ! consistent
             iterdb_file_name =    statefile_output_name
          ELSE                    ! suffix not set so set it:
             iterdb_file_name =    ADJUSTL(statefile_output_name(1:out_namel) // '.txt')
          ENDIF
          CALL iter_dbase_txt
       ENDIF


       RETURN


   END SUBROUTINE write_statefile
