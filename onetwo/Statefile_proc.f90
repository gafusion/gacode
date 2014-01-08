   SUBROUTINE Statefile_proc
!----------------------------------------------------------------------
! -- Process statefile for Onetwo:
!--------------------------------------------------------HSJ-----------
     USE nrtype,                                      ONLY : I4B,DP
     USE iterdbmd_gcnmp,                              ONLY : iterdb_file_name
     USE gcnmp_input,                                 ONLY : write_iterdb_txt,           &
                                                             switch_iterdb_output,       &
                                                             gcnmp_iterdb_filename
                                                             
     USE soln2d,                                      ONLY : ifixshap
     USE mhdcom,                                      ONLY : mhdmethd
     USE iterdbmd,                                    ONLY : iterdb,create_GCNMP_input,  &
                                                             create_XPTOR_input,irwflag, &
                                                             statefile_type
     USE iterdbmd_gcnmp,                              ONLY : irwflag_gcnmp => irwflag
     USE  set_12_gcnmp_vars,                          ONLY : set_gcnmp_vars
     USE file_proc,                                   ONLY : append_time_to_filename,    &
                                                             rename_disk_file
     USE solcon,                                      ONLY : time



     IMPLICIT NONE 

   

     irwflag_gcnmp = 0   ! write new statefile

      IF (iterdb .EQ.  1 .AND. create_XPTOR_input )THEN
          irwflag =0  ! write Xptor  text iterdb file
          CALL iter_dbase129
      ENDIF

!      write(6,FMT='("statefile_type,switch_iterdb_output ,mhdmethd =", &
!            x,i3,x,i3,x,a)')statefile_type,switch_iterdb_output ,mhdmethd ! 8888999cat fort

!      write(6,FMT='("iterdb,ifixshap,switch_iterdb_output,create_gcnmp_imnput=", &
!           x,i3,x,i3,x,i3,x,l5)')iterdb,ifixshap,switch_iterdb_output,create_GCNMP_input 
!      write(6,FMT='("statefile_type =",x,i3)')statefile_type
      irwflag_gcnmp = 0   ! write new statefile
      IF (ifixshap .EQ. 1 .AND. iterdb .EQ.  2 .AND.                                     &
                   .NOT. create_GCNMP_input )  CALL iter_dbase !returns without writing

 


       IF (iterdb .EQ.  1 .AND. create_GCNMP_input .AND. mhdmethd .NE.'tdem' )THEN 
           iterdb_file_name = ADJUSTL(gcnmp_iterdb_filename(1:LEN_TRIM(gcnmp_iterdb_filename)))

          CALL set_gcnmp_vars !set_12_gcnmp_vars.f90
 
          IF(statefile_type == 0 .AND. switch_iterdb_output ==1)THEN
           !here the gcnmp input file was in text mode and the gcnmp
           !output file is in netcdf form
              iterdb_file_name = ADJUSTL(gcnmp_iterdb_filename(1:LEN_TRIM(gcnmp_iterdb_filename)))//'.nc'

              CALL append_time_to_filename (time,iterdb_file_name)

              CALL iter_dbase_nc

              CALL rename_disk_file(iterdb_file_name)
 
          ELSEIF(statefile_type == 0 .AND. switch_iterdb_output ==0)THEN
           !here the gcnmp input file was in text mode and the gcnmp
           !output file is in text  form
              iterdb_file_name = ADJUSTL(gcnmp_iterdb_filename(1:LEN_TRIM(gcnmp_iterdb_filename)))//'.txt'
              CALL append_time_to_filename (time,iterdb_file_name)

              CALL iter_dbase_txt
              CALL rename_disk_file(iterdb_file_name)
           ELSEIF(statefile_type == 1  .AND. switch_iterdb_output ==1)THEN

           !here the gcnmp input file was in netcdf mode and the gcnmp
           !output file is in text form
              iterdb_file_name = ADJUSTL(gcnmp_iterdb_filename(1:LEN_TRIM(gcnmp_iterdb_filename)))//'.txt'
              CALL append_time_to_filename (time,iterdb_file_name)

              CALL iter_dbase_txt
              CALL rename_disk_file(iterdb_file_name)
          ELSEIF(statefile_type == 1 .AND. switch_iterdb_output ==0)THEN
            ! here the input and output files are in netcdf form
              iterdb_file_name = ADJUSTL(gcnmp_iterdb_filename(1:LEN_TRIM(gcnmp_iterdb_filename)))//'.nc'
              CALL append_time_to_filename (time,iterdb_file_name)

              CALL iter_dbase_nc
              CALL rename_disk_file(iterdb_file_name)
          ELSEIF(statefile_type == -1)THEN

             ! here there was no input statefile. write the output according to
             ! setting of write_iterdb_txt
             IF(write_iterdb_txt )THEN
               iterdb_file_name = ADJUSTL(gcnmp_iterdb_filename(1:LEN_TRIM(gcnmp_iterdb_filename)))//'.txt'
               CALL append_time_to_filename (time,iterdb_file_name)

               CALL iter_dbase_txt
               CALL rename_disk_file(iterdb_file_name)
             ELSE
               iterdb_file_name = ADJUSTL(gcnmp_iterdb_filename(1:LEN_TRIM(gcnmp_iterdb_filename)))//'.nc'

               CALL append_time_to_filename (time,iterdb_file_name)
 
               CALL iter_dbase_nc   ! file rw_iterdb_netcdf.F90

               CALL rename_disk_file(iterdb_file_name)
 
             ENDIF
          ENDIF
      ENDIF


      IF (iterdb == -1 .AND.  ifixshap == 1 .AND.  mhdmethd .NE. 'tdem'.AND.  create_GCNMP_input)THEN
          write_iterdb_txt = .FALSE.
          iterdb_file_name = ADJUSTL(gcnmp_iterdb_filename(1:LEN_TRIM(gcnmp_iterdb_filename)))//'.nc'
          CALL append_time_to_filename (time,iterdb_file_name)
          CALL set_gcnmp_vars

          CALL iter_dbase_nc
          CALL rename_disk_file(iterdb_file_name)
      ENDIF

      IF ( iterdb == -1  .AND.  mhdmethd == 'tdem' .AND.  create_GCNMP_input)THEN
          write_iterdb_txt = .FALSE.
          iterdb_file_name = ADJUSTL(gcnmp_iterdb_filename(1:LEN_TRIM(gcnmp_iterdb_filename)))//'_tdem'//'.nc'
          CALL append_time_to_filename (time,iterdb_file_name)
          CALL set_gcnmp_vars

          CALL iter_dbase_nc
          CALL rename_disk_file(iterdb_file_name)
      ENDIF
 
      IF (mhdmethd  ==  'tdem' .AND.iterdb  ==   1  .AND.  create_GCNMP_input  )THEN

          CALL set_gcnmp_vars

          IF(write_iterdb_txt)THEN
             iterdb_file_name = ADJUSTL(gcnmp_iterdb_filename(1:LEN_TRIM(gcnmp_iterdb_filename)))//'_tdem'//'.txt'

             CALL append_time_to_filename (time,iterdb_file_name)

             CALL iter_dbase_txt
             CALL rename_disk_file(iterdb_file_name)
          ELSE
             iterdb_file_name = ADJUSTL(gcnmp_iterdb_filename(1:LEN_TRIM(gcnmp_iterdb_filename)))//'.nc'
             CALL append_time_to_filename (time,iterdb_file_name)

             CALL iter_dbase_nc
             CALL rename_disk_file(iterdb_file_name)
          ENDIF
 
      ENDIF

      RETURN
   END SUBROUTINE Statefile_proc


    SUBROUTINE read_statefile
! -------------------------------------------------------------------------------
! ---
! --- read text or netcdf file for Onetwo initial/boundary conditions.
! --- 
! ---------------------------------------------------------HSJ--6/20/07----------

    USE nrtype,                      ONLY : DP,I4B
    USE iterdbmd,                    ONLY : statefile_name,statefile_type
    USE iterdbmd_gcnmp,              ONLY : irwflag,iterdb_file_name
    USE numbrs,                      ONLY : nj
    USE set_12_gcnmp_vars,           ONLY : set_onetwo_vars
    USE solcon,                      ONLY : time,time0
    USE gcnmp_input,                 ONLY : write_iterdb_txt,switch_iterdb_output

        IMPLICIT NONE

        irwflag = 1                              ! set for read acces
        iterdb_file_name = statefile_name(1:LEN_TRIM(statefile_name))

!        switch_iterdb_output = 1              ! default output is  .nc


!       read statefile,text or netcdf form:
        IF(INDEX(statefile_name,'.nc') == 0)THEN ! regular text file
            CALL iter_dbase_txt                  ! read file using gcnmp variables
            statefile_type = 0                   ! require this value for resart option see Statefile_proc above
           IF(write_iterdb_txt)switch_iterdb_output=0
        ELSE                                     ! input is in netcdf form
                                                 
!          iterdb =3
           CALL iter_dbase_nc                    ! read file using gcnmp variables
           statefile_type = 1                    ! require this value for resart option see Statefile_proc above

           IF(write_iterdb_txt)switch_iterdb_output=1
        ENDIF

           CALL set_onetwo_vars                  ! translate from statefile variables to Onetwo variables (set_12_gcnmp_vars.f90)
           !statefile startup uses time in inone file
!           time = time0


           irwflag = 0                           !assume next access is writting

 
        RETURN

   END     SUBROUTINE read_statefile
