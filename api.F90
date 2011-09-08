 PROGRAM api

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! API                                                             !
! Altitude and Pressure Interpolator                              !
!                                                                 !
! This program is based on                                        !
! ------------------------                                        !
! p_interp v1.0                                                   !
! http://www.mmm.ucar.edu/wrf/src/p_interp.tar.gz                 !
! September 2008 - Cindy Bruyere [NCAR, USA]                      !
!                                                                 !
! Modifications                                                   !
! -------------                                                   !
! - Presentation of the program (cosmetics)                       !
! - Additional routine and arguments in existing routines         !
! - Improve memory management                                     ! 
! - Martian constants + Addition of 'grav' variable               !
! - Interpolation to height above areoid (MOLA zero datum)        !
! October 2009 - Aymeric Spiga [The Open University, UK]          !
! - Interpolation to height above surface                         !
! - Rotated winds (e.g. for polar projections)                    !
! November 2009 - AS                                              !
!                                                                 !
! Purpose                                                         !
! -------                                                         !
!  Program to read wrfout data (possibly several files)           !
!     and interpolate to pressure or height levels                !
!  A new NETCDF file is generated with appropriate suffix         !
!  The program reads namelist.api                                 !
!                                                                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      !
      ! VARIABLES
      !
      CHARACTER (LEN=500)                                :: path_to_input
      CHARACTER (LEN=500)                                :: path_to_output
      CHARACTER (LEN=500)                                :: input_name
      CHARACTER (LEN=500)                                :: output_name
      CHARACTER (LEN=20)                                 :: process
      CHARACTER (LEN=2000)                               :: fields
      REAL, DIMENSION(299)                               :: interp_levels
      INTEGER                                            :: interp_method=1
      INTEGER                                            :: extrapolate=0
      LOGICAL                                            :: debug=.FALSE.
      LOGICAL                                            :: unstagger_grid=.FALSE.
      LOGICAL                                            :: bit64=.FALSE.
      LOGICAL                                            :: oldvar=.TRUE.              

      INTEGER                                            :: funit,ios
      LOGICAL                                            :: is_used

      !
      ! NAMELISTS 
      !
      NAMELIST /io/ path_to_input, input_name, path_to_output, output_name, &
                    process, fields, debug, bit64, oldvar
      NAMELIST /interp_in/ interp_levels, interp_method, extrapolate, unstagger_grid
 
      !
      ! DEFAULT VALUES for VARIABLES
      !
      path_to_input   = './'
      path_to_output  = './'
      output_name     = ' '
      interp_levels   = -99999.
      process         = 'list' !!'all'
      

      !
      ! READ NAMELIST
      !
        DO funit=10,100
           INQUIRE(unit=funit, opened=is_used)
           IF (.not. is_used) EXIT
        END DO
        OPEN(funit,file='namelist.api',status='old',form='formatted',iostat=ios)
        IF ( ios /= 0 ) STOP "ERROR opening namelist.api"
        READ(funit,io)
        READ(funit,interp_in)
        CLOSE(funit)

      !!! MAIN CALL
      CALL api_main ( path_to_input, input_name, path_to_output, output_name, &
                                 process, fields, debug, bit64, oldvar, &
                                 interp_levels, interp_method, extrapolate, unstagger_grid, -99999. ) 

 END PROGRAM api

 SUBROUTINE api_main ( path_to_input, input_name, path_to_output, output_name, &
                       process, fields, debug, bit64, oldvar, &
                       interp_levels, interp_method, extrapolate, unstagger_grid, onelevel )

      IMPLICIT NONE
      INCLUDE 'netcdf.inc'

      !!
      !! EARTH CONSTANTS
      !! 
      !REAL, PARAMETER :: Rd  = 287.04   ! gas constant m2 s-2 K-1
      !REAL, PARAMETER :: Cp  = 7.*Rd/2. ! r=8.314511E+0*1000.E+0/mugaz
      !REAL, PARAMETER :: RCP = Rd/Cp
      !REAL, PARAMETER :: p0  = 100000.
      !REAL, PARAMETER :: grav = 9.81
      !REAL, PARAMETER :: tpot0 = 300.

      !
      ! MARS CONSTANTS
      ! 
      !REAL, PARAMETER :: Rd  = 192.     ! gas constant m2 s-2 K-1
      !REAL, PARAMETER :: Cp  = 844.6    ! r= 8.314511E+0*1000.E+0/mugaz
      REAL, PARAMETER :: Rd  = 191.0
      REAL, PARAMETER :: Cp  = 744.5
      REAL, PARAMETER :: RCP = Rd/Cp
      REAL, PARAMETER :: p0  = 610.
      REAL, PARAMETER :: grav = 3.72
      REAL, PARAMETER :: tpot0 = 220.   ! reference potential temperature

      !
      ! VARIABLES
      !
      CHARACTER,         ALLOCATABLE, DIMENSION(:,:,:,:) :: text
      CHARACTER (LEN=31),ALLOCATABLE, DIMENSION(:)       :: dnamei, dnamej
      CHARACTER(LEN=500),ALLOCATABLE, DIMENSION(:)       :: input_file_names
      CHARACTER(LEN=500),ALLOCATABLE, DIMENSION(:)       :: output_file_names
      DOUBLE PRECISION,  ALLOCATABLE, DIMENSION(:,:,:,:) :: ddata1, ddata2
      REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: data1, data2, data3
      REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: pres_field, pres_out
      REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: pres_stagU, pres_stagV
      REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: ght, phb, qv, tk, rh, tpot, u, v, umet, vmet 
      REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: psfc
      REAL,              ALLOCATABLE, DIMENSION(:,:)     :: ter
      REAL,              ALLOCATABLE, DIMENSION(:,:)     :: longi, lati
      REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: longit, latit
      REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: interm1, interm2
      INTEGER,           ALLOCATABLE, DIMENSION(:)       :: dvali, dvalj
      INTEGER,           ALLOCATABLE, DIMENSION(:,:,:,:) :: idata1, idata2
      INTEGER,                        DIMENSION(4)       :: start_dims = 1
      INTEGER,                        DIMENSION(4)       :: dims_in, dims_out
      INTEGER,                        DIMENSION(6)       :: ishape, jshape
      CHARACTER (LEN=500)                                :: cval !80
      CHARACTER (LEN=31)                                 :: cname, test_dim_name
      CHARACTER (LEN=500)                                :: input_file, output_file, att_text !80
      CHARACTER (LEN=500)                                :: path_to_input
      CHARACTER (LEN=500)                                :: path_to_output
      CHARACTER (LEN=500)                                :: input_name
      CHARACTER (LEN=500)                                :: output_name, tmp_name
      CHARACTER (LEN=10)                                 :: option
      CHARACTER (LEN=132)                                :: command
      CHARACTER (LEN=20)                                 :: process, dummy
      CHARACTER (LEN=2000)                               :: fields, process_these_fields
      REAL, DIMENSION(299)                               :: interp_levels 
      REAL                                               :: rval
      REAL                                               :: MISSING=1.e36
      REAL                                               :: truelat1, truelat2, stand_lon
      INTEGER                                            :: map_proj
      INTEGER                                            :: LINLOG = 1
      INTEGER                                            :: interp_method!=1
      INTEGER                                            :: extrapolate!=0
      INTEGER                                            :: ncid, mcid, rcode
      INTEGER                                            :: idm, ndims, nvars, natt, ngatts
      INTEGER                                            :: nunlimdimid
      INTEGER                                            :: i, ii, j, jj, ix, iy, iweg, isng, ibtg
      INTEGER                                            :: ivar, jvar
      INTEGER                                            :: times_in_file
      INTEGER                                            :: ilen, itype, ival, na, funit, ios
      INTEGER                                            :: num_metgrid_levels
      INTEGER                                            :: number_of_input_files
      INTEGER                                            :: ierr, loop, loslen, strlen, lent
      INTEGER                                            :: is_there
      INTEGER                                            :: kk
      LOGICAL                                            :: is_used
      LOGICAL                                            :: debug!=.FALSE.
      LOGICAL                                            :: interpolate=.FALSE.
      LOGICAL                                            :: unstagger_grid!=.FALSE.
      LOGICAL                                            :: fix_meta_stag=.FALSE.
      LOGICAL                                            :: bit64!=.FALSE.
      LOGICAL                                            :: first=.TRUE.
      LOGICAL                                            :: oldvar!=.FALSE.                   

      REAL :: onelevel
      if ( onelevel .ne. -99999. ) then
        interp_levels(1) = onelevel
        interp_levels(2:) = -99999.
      endif

      !
      ! INPUT FILE NAMES
      !
        lent = len_trim(path_to_input)
         !IF ( path_to_input(lent:lent) /= "/" ) THEN
         !   path_to_input = TRIM(path_to_input)//"/"
         !ENDIF
         !lent = len_trim(path_to_output)
         !IF ( path_to_output(lent:lent) /= "/" ) THEN
         !   path_to_output = TRIM(path_to_output)//"/"
         !ENDIF
        input_name = TRIM(path_to_input)//TRIM(input_name)
        !
        ! BUILD a UNIX COMMAND based on "ls" of all of the input files
        !  
        loslen = LEN ( command )
        CALL all_spaces ( command , loslen ) 
        WRITE ( command , FMT='("ls -1 ",A," > .foo")' ) TRIM ( input_name )

      !
      ! NUMBER OF INPUTS FILES
      !  
        !  We stuck all of the matching files in the ".foo" file.  Now we place the 
        !  number of the those file (i.e. how many there are) in ".foo1". 
        CALL SYSTEM ( TRIM ( command ) ) 
        CALL SYSTEM ( '( cat .foo | wc -l > .foo1 )' )
        !  Read the number of files.
        OPEN (FILE   = '.foo1'       , &
              UNIT   = 112           , &
              STATUS = 'OLD'         , &
              ACCESS = 'SEQUENTIAL'  , &
              FORM   = 'FORMATTED'     )
        READ ( 112 , * ) number_of_input_files
        CLOSE ( 112 )
        !  If there are zero files, we are toast.
        IF ( number_of_input_files .LE. 0 ) THEN
           print*, ' Oops, we need at least ONE input file for the program to read.'
           print*, '       Make sure you have the path, and name file(s) correct,'
           print*, '       including wild characters if needed.'
           STOP
        END IF

      !
      ! GET FILENAMES IN THE PROGRAM  
      !
        !  Allocate space for this many files.
        ALLOCATE (  input_file_names(number_of_input_files) , STAT=ierr )
        ALLOCATE ( output_file_names(number_of_input_files) , STAT=ierr )
        !  Did the allocate work OK?
        IF ( ierr .NE. 0 ) THEN
           print*, ' tried to allocate ', number_of_input_files, ' input files, (look at ./foo)'
           STOP
        END IF
        !  Initialize all of the file names to blank.
          input_file_names = '                                                  ' // &
                             '                                                  ' // &
                             '                                '
         output_file_names = '                                                  ' // &
                             '                                                  ' // &
                             '                                '
        !  Open the file that has the list of filenames.
        OPEN (FILE   = '.foo'        , &
              UNIT   = 111           , &
              STATUS = 'OLD'         , &
              ACCESS = 'SEQUENTIAL'  , &
              FORM   = 'FORMATTED'     )
        !  Read all of the file names and store them - define output file name.
        LINLOG = interp_method
        DO loop = 1 , number_of_input_files
           READ ( 111 , FMT='(A)' ) input_file_names(loop)
           IF ( output_name == ' ' ) THEN
              ilen = INDEX(TRIM(input_file_names(loop)),'/',.TRUE.) 
              output_file_names(loop) = TRIM(path_to_output)//input_file_names(loop)(ilen+1:)
              IF ( LINLOG .le. 2 ) output_file_names(loop) = TRIM(output_file_names(loop))//"_p"
              IF ( LINLOG  ==  3 ) output_file_names(loop) = TRIM(output_file_names(loop))//"_z"
              IF ( LINLOG  ==  4 ) output_file_names(loop) = TRIM(output_file_names(loop))//"_zabg"
           ELSE
              IF ( number_of_input_files == 1 ) THEN
                 output_file_names(loop) = TRIM(path_to_output)//TRIM(output_name)
              ELSE
                 write(tmp_name,'(A,A,"_",I4.4)') TRIM(path_to_output), TRIM(output_name), loop
                 output_file_names(loop) = tmp_name
              ENDIF
           ENDIF
        END DO
        CLOSE ( 111 )  !!AS: don't know why but it was written 112... fixed the bug.
        print*, " " 
        !   We clean up our own messes.
        CALL SYSTEM ( '/bin/rm -f .foo'  )
        CALL SYSTEM ( '/bin/rm -f .foo1' )

      !
      ! LIST OF FIELD WANTED for OUTPUT  
      !
        ! Do we have a list of field that we want on output?
!        process_these_fields = ','
!!!! ICI CHAMPS SORTIS PAR DEFAULT
!!!! ICI CHAMPS SORTIS PAR DEFAULT
!!!! ICI CHAMPS SORTIS PAR DEFAULT 
        process_these_fields = ',Times,XLAT,XLONG,' !! les premiere et derniere virgules sont importantes
!!!! ICI CHAMPS SORTIS PAR DEFAULT
!!!! ICI CHAMPS SORTIS PAR DEFAULT
!!!! ICI CHAMPS SORTIS PAR DEFAULT
        IF ( INDEX(process,'list') /= 0) THEN
          DO i = 1 , len(fields)
            IF (fields(i:i) /= ' ' ) THEN
              process_these_fields = trim(process_these_fields)//fields(i:i)
            ENDIF
          END DO
          process_these_fields = trim(process_these_fields)//","
        END IF

      !
      ! HELLO to WORLD  
      !
      write(6,*) 
      write(6,*) "##############################################################"
      write(6,'(A,i4,A)') " RUNNING api on ", number_of_input_files, " file(s)."
      write(6,*) 
      write(6,'(A,$)') " A) INTERPOLATION METHOD: "
      IF ( LINLOG == 1 ) write(6,*) "  Pressure - linear in p"
      IF ( LINLOG == 2 ) write(6,*) "  Pressure - linear in log p"
      IF ( LINLOG == 3 ) write(6,*) "  Height above areoid (MOLA zero datum)"
      IF ( LINLOG == 4 ) write(6,*) "  Height above surface"
      write(6,'(A,$)') " B) VERTICAL EXTRAPOLATION: "
      IF (extrapolate == 0) write(6,*) "  BELOW GROUND and ABOVE model top will be set to missing values "
      IF (extrapolate == 1) write(6,*) "  BELOW GROUND will be extrapolated and ABOVE model top will be set to values at model top"
      write(6,'(A,$)') " C) HORIZONTAL GRID: "
      IF (.not. unstagger_grid) write(6,*) "  Data will be output on C-grid " 
      IF (unstagger_grid) write(6,*) "  Data will be output on unstaggered grid "

      !
      ! GET LEVELS (pressure or altitude) TO INTERPOLATE TO
      !
      write(6,*)
      IF ( LINLOG .le. 2 ) write(6,*) "INTERPOLATING TO PRESSURE LEVELS: "
      IF ( LINLOG  ==  3 ) write(6,*) "INTERPOLATING TO ALTITUDE LEVELS: "
      IF ( LINLOG  ==  4 ) write(6,*) "INTERPOLATING TO ABG LEVELS: "
      num_metgrid_levels = 0
      DO
        IF (interp_levels(num_metgrid_levels+1) == -99999.) EXIT
        num_metgrid_levels = num_metgrid_levels + 1
        if (mod(num_metgrid_levels,8) /= 0 )write(6,'(f8.3,$)') interp_levels(num_metgrid_levels)
        if (mod(num_metgrid_levels,8) == 0 )write(6,'(f8.3)') interp_levels(num_metgrid_levels)
        IF ( LINLOG .le. 2 ) interp_levels(num_metgrid_levels) = interp_levels(num_metgrid_levels) * 100.0   !!! Pa
        IF ( LINLOG .gt. 2 ) interp_levels(num_metgrid_levels) = interp_levels(num_metgrid_levels) * 1000.0  !!! km
      END DO
      write(6,*)
      write(6,*)

      !
      ! LOOP on FILES
      !
      DO loop = 1, number_of_input_files

        IF (debug) write(6,*) "##############################################################"
        input_file  = input_file_names(loop)
        output_file = output_file_names(loop)

        IF ( .not. debug ) write(6,*) " Output will be written to: ",trim(output_file)

        IF (debug) THEN
          write(6,*) " INPUT FILE:         ",trim(input_file)
          write(6,*) " OUTPUT FILE:        ",trim(output_file)
          write(6,*) "  "
        ENDIF

      !
      ! OPEN INPUT AND OUTPUT FILE
      !
        rcode = nf_open(input_file, 0, ncid)
        if (rcode .ne. nf_noerr) call handle_err(rcode)
        if (bit64) then
          rcode = nf_create(output_file, NF_64BIT_OFFSET, mcid)
        else
          rcode = nf_create(output_file, 0, mcid)
        endif
        if (rcode .ne. nf_noerr) call handle_err(rcode)

      !
      ! GET BASIC INFORMATION ABOUT THE FILE
      ! most important 
      !   ndims:  number of dimensions
      !   nvars:  number of variables
      !   ngatts: number of global attributes
        rcode = nf_inq(ncid, ndims, nvars, ngatts, nunlimdimid)
        if (rcode .ne. nf_noerr) call handle_err(rcode)
        IF (debug) THEN
          write(6,*) ' INPUT file has = ',ndims, ' dimensions, '
          write(6,*) '                  ',nvars, ' variables, and '      
          write(6,*) '                  ',ngatts,' global attributes '
          write(6,*) "  "
        ENDIF
        rcode = nf_get_att_int (ncid, nf_global, 'WEST-EAST_GRID_DIMENSION', iweg)
        rcode = nf_get_att_int (ncid, nf_global, 'SOUTH-NORTH_GRID_DIMENSION', isng)
        rcode = nf_get_att_int (ncid, nf_global, 'BOTTOM-TOP_GRID_DIMENSION', ibtg)

      !  
      ! ALLOCATE SOME VARIABLES
      !
        IF (ALLOCATED(dnamei)) deallocate(dnamei)
            ALLOCATE (dnamei(20))
        IF (ALLOCATED(dnamej)) deallocate(dnamej)
            ALLOCATE (dnamej(20))
        IF (ALLOCATED(dvali)) deallocate(dvali)
            ALLOCATE (dvali(20))
        IF (ALLOCATED(dvalj)) deallocate(dvalj)
            ALLOCATE (dvalj(20))

      !      
      ! READ ALL DIMS FROM INPUT FILE AND CREATE SOME DIMS FOR OUTPUT FILE
      !
        j = 0
        DO i = 1, ndims
          rcode = nf_inq_dim(ncid, i, dnamei(i), dvali(i))
          IF ( dnamei(i) == "Time" ) THEN
            j = j + 1
            dnamej(j) = dnamei(i)
            dvalj(j) = dvali(i)
            rcode = nf_def_dim(mcid, dnamej(j), NF_UNLIMITED, j)
            times_in_file  = dvali(i)
          ENDIF 
        ENDDO
        !!! Create the new height/pressure dims
        j = j + 1
        IF ( LINLOG .le. 2 ) THEN 
              dnamej(j) = 'pressure' 
              dvalj(j) = num_metgrid_levels
        ELSE IF ( LINLOG .eq. 3 ) THEN 
              dnamej(j) = 'altitude' !'bottom_top' !'altitude'
              dvalj(j) = num_metgrid_levels
        ELSE 
              dnamej(j) = 'altitude_abg' !'bottom_top' !'altitude_abg'
              dvalj(j) = num_metgrid_levels
        ENDIF
        rcode = nf_def_dim(mcid, dnamej(j), dvalj(j), j)

      !
      ! DEALING WITH THE GLOBAL ATTRIBUTES
      !
        IF (debug) THEN
          write(6,*) 
          write(6,*) " OUTPUT FILE attributes:"
        ENDIF
        do i = 1, ngatts
          rcode = nf_inq_attname(ncid, nf_global, i,    cname)
          rcode = nf_inq_atttype(ncid, nf_global, cname, itype)
          rcode = nf_inq_attlen (ncid, nf_global, cname, ilen)
         if ( itype .eq. 2 ) then        ! characters
            rcode = nf_get_att_text (ncid, nf_global, cname, cval)
            if(cname(1:5) .eq. 'TITLE') then
               IF ( LINLOG .le. 2 ) THEN
               cval = cval(1:ilen)//" - ON PRES LEVELS"
               ELSE
               cval = cval(1:ilen)//" - ON z LEVELS"
               ENDIF
               ilen = len_trim(cval)
            endif
            IF (debug) &
            write(6,'("     i = ",i2," : ",A," = ",A)') &
                    i,cname,cval(1:ilen)
            rcode = nf_put_att_text(mcid, nf_global, cname, ilen,&
                      cval(1:ilen))
          elseif ( itype .eq. 4 ) then     ! integers
            rcode = nf_get_att_int (ncid, nf_global, cname, ival)
            IF ( INDEX(cname,'BOTTOM-TOP_PATCH') == 0 ) THEN
               IF (cname .eq. 'BOTTOM-TOP_GRID_DIMENSION') ival = num_metgrid_levels
               IF (debug) &
                   write(6,'("     i = ",i2," : ",A," = ",i7)') &
                       i,cname,ival        
               rcode = nf_put_att_int(mcid, nf_global, cname, itype,&
                         ilen, ival)
             ENDIF
            IF (cname .eq. 'MAP_PROJ') map_proj = ival 
           elseif ( itype .eq. 5 ) then    ! real
            rcode = nf_get_att_real (ncid, nf_global, cname, rval)
            IF (debug) &
               write(6,'("     i = ",i2," : ",A," = ",G18.10E2)') &
                    i,cname,rval
               rcode = nf_put_att_real(mcid, nf_global, cname, itype,&
                      ilen, rval)
            IF (cname .eq. 'TRUELAT1') truelat1 = rval
            IF (cname .eq. 'TRUELAT2') truelat2 = rval
            IF (cname .eq. 'STAND_LON') stand_lon = rval
          end if
        enddo
       rcode = nf_enddef(mcid)

      ! 
      ! WE NEED SOME BASIC FIELDS
      !
      ! --> P_TOP
      !
      IF ( LINLOG .le. 2 ) THEN
         IF ( DEBUG ) PRINT *, 'P_TOP'     
         IF (ALLOCATED(data1)) deallocate(data1)
         allocate (data1(times_in_file,1,1,1))
         rcode = nf_inq_varid    ( ncid, "P_TOP", i )
         rcode = nf_get_var_real ( ncid, i, data1 )
         IF ( first ) THEN
            IF ( extrapolate == 1 .AND. &
                 (data1(1,1,1,1)-interp_levels(num_metgrid_levels)) > 0.0 ) THEN
               write(6,*)
               write(6,*) "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
               write(6,*) " WARNING: Highest requested pressure level is above PTOP."
               write(6,'(A,F7.2,A)') "           Use all pressure level data above", data1(1,1,1,1)*.01, " mb"
               write(6,*) "          with caution."
               write(6,*) "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
            ENDIF
            first = .FALSE.
         ENDIF
         deallocate (data1)
      ENDIF   
      !
      ! --> INTERPOLATION LEVELS
      !
         IF (ALLOCATED(pres_out)) deallocate(pres_out)
         allocate (pres_out(iweg-1, isng-1, num_metgrid_levels, times_in_file))
         do i = 1, num_metgrid_levels
           pres_out (:,:,i,:) = interp_levels(i)
         enddo
      !
      ! --> LONGITUDES + everything for ROTATED WINDS
      !
      IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'uvmet') /= 0 ) THEN   
        IF (ALLOCATED(longi)) deallocate(longi)
           IF (ALLOCATED(longit)) deallocate(longit)
           allocate (longi(iweg-1, isng-1))
           allocate (longit(iweg-1, isng-1, times_in_file))
           IF ( DEBUG ) PRINT *, 'LONGITUDE' 
           rcode = nf_inq_varid    ( ncid, "XLONG", i )
           rcode = nf_get_var_real ( ncid, i, longit )
           longi = longit(:,:,1)
           deallocate(longit)
        IF (ALLOCATED(lati)) deallocate(lati)
           IF (ALLOCATED(latit)) deallocate(latit)
           allocate (lati(iweg-1, isng-1))
           allocate (latit(iweg-1, isng-1, times_in_file))
           IF ( DEBUG ) PRINT *, 'LATITUDE' 
           rcode = nf_inq_varid    ( ncid, "XLAT", i )
           rcode = nf_get_var_real ( ncid, i, latit )
           lati = latit(:,:,1)
           deallocate(latit)
        IF (ALLOCATED(u))       deallocate(u)
        IF (oldvar) THEN
           allocate (u      (iweg,   isng-1, ibtg-1, times_in_file ))
           IF ( DEBUG ) PRINT *, 'U COMPONENT'
           rcode = nf_inq_varid    ( ncid, "U", i )
           rcode = nf_get_var_real ( ncid, i, u )
        ELSE
           allocate (u      (iweg-1, isng-1, ibtg-1, times_in_file ))
           IF ( DEBUG ) PRINT *, 'UAVE COMPONENT'
           rcode = nf_inq_varid    ( ncid, "UAVE", i )
           rcode = nf_get_var_real ( ncid, i, u )
        ENDIF
        IF (ALLOCATED(v))       deallocate(v)
           IF (oldvar) THEN
           allocate (v      (iweg-1, isng,   ibtg-1, times_in_file ))
           IF ( DEBUG ) PRINT *, 'V COMPONENT'
           rcode = nf_inq_varid    ( ncid, "V", i )
           rcode = nf_get_var_real ( ncid, i, v )
        ELSE
           allocate (v      (iweg-1, isng-1, ibtg-1, times_in_file ))
           IF ( DEBUG ) PRINT *, 'VAVE COMPONENT'
           rcode = nf_inq_varid    ( ncid, "VAVE", i )
           rcode = nf_get_var_real ( ncid, i, v )
        ENDIF
        IF (ALLOCATED(umet))    deallocate(umet)
           allocate (umet   (iweg-1, isng-1, ibtg-1, times_in_file ))
        IF (ALLOCATED(vmet))    deallocate(vmet)
           allocate (vmet   (iweg-1, isng-1, ibtg-1, times_in_file ))
        IF (ALLOCATED(interm1)) deallocate(interm1)
           allocate (interm1(iweg-1, isng-1, ibtg-1 ))
        IF (ALLOCATED(interm2)) deallocate(interm2)
           allocate (interm2(iweg-1, isng-1, ibtg-1 ))
        IF ( DEBUG ) PRINT *, 'ROTATE WINDS'
        write(6,*) "  Data will be output on unstaggered grid "
        do kk = 1, times_in_file
         !IF ( DEBUG ) print *, kk
         IF (oldvar) THEN
           interm1(1:iweg-1,:,:) = ( u(1:iweg-1,:,:,kk) + u(2:iweg,:,:,kk) ) * .5
           interm2(:,1:isng-1,:) = ( v(:,1:isng-1,:,kk) + v(:,2:isng,:,kk) ) * .5
         ELSE
           interm1(1:iweg-1,:,:) = u(1:iweg-1,:,:,kk) 
           interm2(:,1:isng-1,:) = v(:,1:isng-1,:,kk) 
         ENDIF
         if (unstagger_grid) map_proj=3 !!! AS: unstagger_grid=T makes the program to put unstaggered U and V in uvmet 
         CALL calc_uvmet( interm1, interm2, umet(:,:,:,kk), vmet(:,:,:,kk),  &
                    truelat1, truelat2, stand_lon, map_proj, longi, lati, iweg-1, isng-1, ibtg-1 )
        enddo
        deallocate(lati)
        deallocate(longi)
        deallocate(u) 
        deallocate(v)
        deallocate(interm1)
        deallocate(interm2)
      ENDIF

      !   
      ! --> GEOPOTENTIAL HEIGHT (in mode 1-2, no need if geopotential is not requested)
      !
      IF ( LINLOG .gt. 2 .OR. INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'GHT') /= 0 ) THEN
         IF ( DEBUG ) PRINT *, 'GEOPOTENTIAL'
         !IF (ALLOCATED(ght)) deallocate(ght)
         !allocate (ght(iweg-1, isng-1, ibtg-1, times_in_file))
         IF (ALLOCATED(phb)) deallocate(phb)
         allocate (phb(iweg-1, isng-1, ibtg, times_in_file ))
                   !         IF (ALLOCATED(data1)) deallocate(data1)
                   !         allocate (data1(iweg-1, isng-1, ibtg, times_in_file ))
                   !         rcode = nf_inq_varid    ( ncid, "PH", i )
                   !         rcode = nf_get_var_real ( ncid, i, data1 )
                   !         rcode = nf_inq_varid    ( ncid, "PHB", i )
                   !         rcode = nf_get_var_real ( ncid, i, phb )
                   !         data1 = (data1 + phb) 
                   !         ght(:,:,1:ibtg-1,:) = ( data1(:,:,1:ibtg-1,:) + data1(:,:,2:ibtg,:) )*.5
                   !         deallocate (data1)
         rcode = nf_inq_varid    ( ncid, "PHTOT", i )
         rcode = nf_get_var_real ( ncid, i, phb )
         !ght(:,:,1:ibtg-1,:) = phb(:,:,1:ibtg-1,:)   !! costly in memory !!
         !deallocate (phb)
      ENDIF
      !   
      ! --> PRESSURE (in mode 3-4, no need if regular temperature is not requested) 
      !
      IF ( LINLOG .le. 2 .OR. INDEX(process,'all') /= 0 &
                .OR. INDEX(process_these_fields,'tk') /= 0 &
                .OR. INDEX(process_these_fields,'tpot') /= 0 ) THEN
         IF ( DEBUG ) PRINT *, 'PRESSURE'     
         IF (ALLOCATED(pres_field)) deallocate(pres_field)
         allocate (pres_field(iweg-1, isng-1, ibtg-1, times_in_file ))
                !         IF (ALLOCATED(data1)) deallocate(data1)
                !         allocate (data1(iweg-1, isng-1, ibtg-1, times_in_file ))
                !         rcode = nf_inq_varid    ( ncid, "P", i )
                !         rcode = nf_get_var_real ( ncid, i, pres_field )
                !         rcode = nf_inq_varid    ( ncid, "PB", i )
                !         rcode = nf_get_var_real ( ncid, i, data1 )
                !         pres_field = pres_field + data1
         rcode = nf_inq_varid    ( ncid, "PTOT", i )
         rcode = nf_get_var_real ( ncid, i, pres_field )
                !         deallocate (data1)
       !
       ! --> REGULAR and POTENTIAL TEMPERATURE
       ! 
         IF ( DEBUG ) PRINT *, 'TEMP'
         IF (ALLOCATED(tpot)) deallocate(tpot)
         allocate (tpot(iweg-1, isng-1, ibtg-1, times_in_file ))
         IF (oldvar) THEN 
            rcode = nf_inq_varid    ( ncid, "T", i )
            rcode = nf_get_var_real ( ncid, i, tpot )
            tpot = tpot + tpot0
         ELSE
            rcode = nf_inq_varid    ( ncid, "TAVE", i )
            rcode = nf_get_var_real ( ncid, i, tpot )
         ENDIF
         IF ( LINLOG .le. 2 .OR. INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'tk') /= 0 ) THEN      
           IF (ALLOCATED(tk)) deallocate(tk)
           allocate (tk(iweg-1, isng-1, ibtg-1, times_in_file ))
           tk = tpot * ( pres_field / p0 )**RCP
         ENDIF
         IF (INDEX(process_these_fields,'tpot') == 0 .AND. INDEX(process,'all') == 0) deallocate(tpot)  
       ENDIF  
       !
       ! --> SURFACE PRESSURE
       ! 
       IF ( LINLOG .le. 2 ) THEN
         IF ( DEBUG ) PRINT *, 'PSFC'      
         IF (ALLOCATED(psfc)) deallocate(psfc)
         allocate (psfc(iweg-1, isng-1, times_in_file ))
         IF (ALLOCATED(data1)) deallocate(data1)
         allocate (data1(iweg-1, isng-1, 1, times_in_file ))
         rcode = nf_inq_varid    ( ncid, "PSFC", i )
         rcode = nf_get_var_real ( ncid, i, data1 )
         psfc(:,:,:) = data1(:,:,1,:)
         deallocate (data1)
       ENDIF
       !
       ! --> TOPOGRAPHY
       ! 
       IF ( LINLOG .ne. 3 ) THEN
         IF ( DEBUG ) PRINT *, 'TOPO'      
         IF (ALLOCATED(ter)) deallocate(ter)
         allocate (ter(iweg-1, isng-1))
         IF (ALLOCATED(data1)) deallocate(data1)
         allocate (data1(iweg-1, isng-1, 1, times_in_file ))
         rcode = nf_inq_varid    ( ncid, "HGT", i )
         rcode = nf_get_var_real ( ncid, i, data1 )
         ter(:,:) = data1(:,:,1,1)
         IF ( LINLOG .ne. 4 ) deallocate (data1)
       ENDIF  

                 !         IF (ALLOCATED(qv)) deallocate(qv)
                 !         allocate (qv(iweg-1, isng-1, ibtg-1, times_in_file ))
                 !         rcode = nf_inq_varid    ( ncid, "QVAPOR", i )
                 !         rcode = nf_get_var_real ( ncid, i, qv )
                 !
                 ! NOT SUITABLE for MARS
                 !
                 !         IF (ALLOCATED(rh)) deallocate(rh)
                 !         allocate (rh(iweg-1, isng-1, ibtg-1, times_in_file ))
                 !         IF (ALLOCATED(data1)) deallocate(data1)
                 !         IF (ALLOCATED(data2)) deallocate(data2)
                 !         allocate (data1(iweg-1, isng-1, ibtg-1, times_in_file ))
                 !         allocate (data2(iweg-1, isng-1, ibtg-1, times_in_file ))
                 !         data1 = 10.*0.6112*exp(17.67*(tk-273.16)/(TK-29.65))
                 !         data2 = 0.622*data1/(0.01 * pres_field -  (1.-0.622)*data1)
                 !         rh    = 100.*AMAX1(AMIN1(qv/data2,1.0),0.0)
                 !         deallocate (data1)
                 !         deallocate (data2)
       !   
       ! --> USEFUL for INTERPOLATION
       !
           !!! at that point pres_field is not used for calculation anymore
           !!! it is used to set the initial coordinate for interpolation
           IF ( LINLOG == 3) THEN
                   IF (ALLOCATED(pres_field)) deallocate(pres_field)
                   allocate ( pres_field(iweg-1, isng-1, ibtg-1, times_in_file) )
                   pres_field = phb(:,:,1:ibtg-1,:) / grav
           ENDIF
           IF ( LINLOG == 4) THEN
                   IF (ALLOCATED(pres_field)) deallocate(pres_field)
                   allocate ( pres_field(iweg-1, isng-1, ibtg-1, times_in_file))
                DO kk = 1, ibtg-1
                   pres_field(:,:,kk,:) = - data1(:,:,1,:) + phb(:,:,kk,:) / grav
                ENDDO
                   deallocate (data1)
                   !PRINT *, pres_field(10,10,:,1)
           ENDIF
 

       IF (ALLOCATED(pres_stagU)) deallocate(pres_stagU)
       IF (ALLOCATED(pres_stagV)) deallocate(pres_stagV)
       IF ( DEBUG ) PRINT *, 'STAGGERED COORDINATES'
           allocate (pres_stagU(iweg, isng-1, ibtg-1, times_in_file ))
           allocate (pres_stagV(iweg-1, isng, ibtg-1, times_in_file ))
           pres_stagU(1,:,:,:)        =  pres_field(1,:,:,:)
           pres_stagU(iweg,:,:,:)     =  pres_field(iweg-1,:,:,:)
           pres_stagU(2:iweg-1,:,:,:) = (pres_field(1:iweg-2,:,:,:) + pres_field(2:iweg-1,:,:,:))*.5
           pres_stagV(:,1,:,:)        =  pres_field(:,1,:,:)
           pres_stagV(:,isng,:,:)     =  pres_field(:,isng-1,:,:)
           pres_stagV(:,2:isng-1,:,:) = (pres_field(:,1:isng-2,:,:) + pres_field(:,2:isng-1,:,:))*.5

       !-----------------------------
       !-----------------------------   
       ! TRAIN FILE 
       !-----------------------------
       !-----------------------------   
        IF (debug) THEN
          write(6,*) 
          write(6,*) 
          write(6,*) "FILE variables:"
        ENDIF
        jvar = 0
        loop_variables : DO ivar = 1, nvars

          rcode = nf_inq_var(ncid, ivar, cval, itype, idm, ishape, natt)

          !!! Do we want this variable ?
                                   !          IF ( trim(cval) == 'P'  .OR. trim(cval) == 'PB' )  CYCLE loop_variables
                                   !          IF ( trim(cval) == 'PH' .OR. trim(cval) == 'PHB' ) CYCLE loop_variables
          !IF ( trim(cval) == 'PTOT' )                         CYCLE loop_variables ! PTOT: no
          IF ( trim(cval) == 'PHTOT' )                        CYCLE loop_variables ! PHTOT: no
          IF ( trim(cval) == 'T' )                            CYCLE loop_variables ! T: no (tk and tpot calculated above)
          IF ( unstagger_grid .AND. (INDEX(cval,'_U') /= 0) ) CYCLE loop_variables !!! no sense in keeping these
          IF ( unstagger_grid .AND. (INDEX(cval,'_V') /= 0) ) CYCLE loop_variables !!! no sense in keeping these
          IF ( INDEX(process,'all') == 0 ) THEN
             !!! Only want some variables - see which
             dummy = ","//trim(cval)//","
             is_there = INDEX(process_these_fields,trim(dummy))
             IF ( is_there == 0 ) THEN
                IF ( debug ) print*,"NOTE: ", trim(cval), " - Not requested"
                CYCLE loop_variables !!! don't want this one
             ENDIF
          ENDIF

          IF ( idm >= 4 .AND. itype == 4 ) THEN
            print*,"NOTE: We cannot deal with 3D integers - maybe later"
            CYCLE loop_variables
          ENDIF
          IF ( itype == 6 ) THEN
            print*,"NOTE: We cannot deal with double precision data - maybe later"
            CYCLE loop_variables
          ENDIF
          IF ( itype == 2 .OR. itype == 4 .OR. itype == 5 ) THEN
            !!! OK I know what to do this this
          ELSE
            print*,"NOTE: Do not understand this data type ", itype, " skip field."
            CYCLE loop_variables
          ENDIF

          !IF ( trim(cval) == 'U' )   cval = 'UU'
          !IF ( trim(cval) == 'V' )   cval = 'VV'
          !IF ( trim(cval) == 'T' )   cval = 'TT'

          !!! OK - we want this - lets continue
          jvar = jvar + 1
          jshape = 0
          interpolate = .FALSE.
          fix_meta_stag = .FALSE.
          rcode = nf_redef(mcid)
          DO ii = 1, idm 
            test_dim_name = dnamei(ishape(ii))
            IF ( test_dim_name == 'bottom_top' .OR. test_dim_name == 'bottom_top_stag' ) THEN
                 IF ( test_dim_name == 'bottom_top_stag' ) fix_meta_stag = .TRUE.
                 IF (LINLOG .le. 2) test_dim_name = 'pressure'
                 IF (LINLOG .eq. 3) test_dim_name = 'altitude' !'bottom_top' !'altitude' 
                 IF (LINLOG .eq. 4) test_dim_name = 'altitude_abg' !'bottom_top' !'altitude_abg'
                 interpolate = .TRUE.
            ENDIF
            IF ( unstagger_grid .AND. test_dim_name == 'west_east_stag' )   THEN
               test_dim_name = 'west_east'
               fix_meta_stag = .TRUE.
            ENDIF
            IF ( unstagger_grid .AND. test_dim_name == 'south_north_stag' ) THEN
               test_dim_name = 'south_north'
               fix_meta_stag = .TRUE.
            ENDIF
            DO jj = 1,j
              IF ( test_dim_name == dnamej(jj) ) THEN
                jshape(ii) = jj
              ENDIF
            ENDDO

            IF ( jshape(ii) == 0 ) THEN
              j = j + 1
              jshape(ii) = j
              dnamej(j) = dnamei(ishape(ii))
              dvalj(j) = dvali(ishape(ii))
              rcode = nf_def_dim(mcid, dnamej(j), dvalj(j), j)
            ENDIF
          ENDDO
          rcode = nf_def_var(mcid, cval, itype, idm, jshape, jvar)
          !print *, 'mcid',   mcid
          !print *, 'cval',   cval
          !print *, 'itype',  itype
          !print *, 'idm',    idm
          !print *, 'jshape', jshape
          !print *, 'jvar',   jvar  

          DO na = 1, natt
             rcode = nf_inq_attname(ncid, ivar, na, cname)
             IF ( fix_meta_stag .AND. trim(cname) == 'stagger' ) THEN
                att_text = "-"
                ilen = len_trim(att_text)
                rcode = nf_put_att_text(mcid, jvar, cname, ilen, att_text(1:ilen) )
             ELSEIF ( fix_meta_stag .AND. trim(cname) == 'coordinates' ) THEN
                att_text = "XLONG XLAT"
                ilen = len_trim(att_text)
                rcode = nf_put_att_text(mcid, jvar, cname, ilen, att_text(1:ilen) )
             ELSE
                rcode = nf_copy_att(ncid, ivar, cname, mcid, jvar)
             ENDIF
          ENDDO

          IF ( extrapolate == 0 ) THEN
            rcode = nf_put_att_real(mcid, jvar, 'missing_value', NF_FLOAT, 1, MISSING )
          ENDIF

          rcode = nf_enddef(mcid)

!          
! GET THE DIMS FOR INPUT AND OUTPUT FROM THE SHAPE
!
          dims_in  = 1
          dims_out = 1
          DO ii = 1,idm
            dims_in(ii)  = dvali(ishape(ii))
            dims_out(ii) = dvalj(jshape(ii))
          ENDDO
          IF (debug) write(6,*) 'VAR: ',trim(cval)
          IF (debug) THEN
            write(6,*) '     DIMS  IN: ',dims_in
            write(6,*) '     DIMS OUT: ',dims_out
          ENDIF

!
! ALLOCATE THE INPUT AND OUTPUT ARRAYS
! READ THE DATA FROM INPUT FILE
!
  IF     (itype == 2) THEN          ! character
    allocate (text(dims_in(1), dims_in(2), dims_in(3), dims_in(4)))
    rcode = nf_get_var_text(ncid, ivar, text)
    rcode = nf_put_vara_text (mcid, jvar, start_dims, dims_in, text)
    IF (debug) write(6,*) '     SAMPLE VALUE = ',text(:,1,1,1)
    deallocate (text)
 
  ELSEIF (itype == 4) THEN          ! integer
    allocate (idata1(dims_in(1), dims_in(2), dims_in(3), dims_in(4)))
    rcode = nf_get_var_int(ncid, ivar, idata1)
    rcode = nf_put_vara_int (mcid, jvar, start_dims, dims_in, idata1)
    IF (debug) write(6,*) '     SAMPLE VALUE = ',idata1(dims_in(1)/2,dims_in(2)/2,1,1)
    deallocate (idata1)
  
  ELSEIF (itype == 5) THEN          ! real
     allocate (data1(dims_in(1), dims_in(2), dims_in(3), dims_in(4)))
     allocate (data2(dims_out(1),dims_out(2),dims_out(3),dims_out(4)))
     rcode = nf_get_var_real(ncid, ivar, data1)
  
            IF (idm >= 4 .AND. interpolate) THEN  
               IF (debug) write(6,*) '     THIS IS A FIELD WE NEED TO INTERPOLATE'       

               IF ( dims_in(1) == iweg .AND. .not. unstagger_grid ) THEN
                  CALL interp (data2, data1, pres_stagU, interp_levels, psfc, ter, tk, qv,   &
                                 iweg, isng-1, ibtg-1, dims_in(4), &
                                 num_metgrid_levels, LINLOG, extrapolate, .FALSE., MISSING)
               ELSEIF ( dims_in(2) == isng .AND. .not. unstagger_grid ) THEN
                  CALL interp (data2, data1, pres_stagV, interp_levels, psfc, ter, tk, qv,   &
                                 iweg-1, isng, ibtg-1, dims_in(4), &
                                 num_metgrid_levels, LINLOG, extrapolate, .FALSE., MISSING)
               ELSEIF ( dims_in(1) == iweg .AND. unstagger_grid ) THEN
                  allocate (data3(iweg-1, isng-1, ibtg-1, dims_in(4)))
                  data3(1:iweg-1,:,:,:) = (data1(1:iweg-1,:,:,:) + data1(2:iweg,:,:,:)) * .5
                  CALL interp (data2, data3, pres_field, interp_levels, psfc, ter, tk, qv,   &
                                 iweg-1, isng-1, ibtg-1, dims_in(4), &
                                 num_metgrid_levels, LINLOG, extrapolate, .FALSE., MISSING)
                  deallocate(data3)
               ELSEIF ( dims_in(2) == isng .AND. unstagger_grid ) THEN
                  allocate (data3(iweg-1, isng-1, ibtg-1, dims_in(4)))
                  data3(:,1:isng-1,:,:) = (data1(:,1:isng-1,:,:) + data1(:,2:isng,:,:)) * .5
                  CALL interp (data2, data3, pres_field, interp_levels, psfc, ter, tk, qv,   &
                                 iweg-1, isng-1, ibtg-1, dims_in(4), &
                                 num_metgrid_levels, LINLOG, extrapolate, .FALSE., MISSING)
                  deallocate(data3)
               ELSEIF ( dims_in(3) == ibtg ) THEN
                  allocate (data3(iweg-1, isng-1, ibtg-1, dims_in(4)))
                  IF (debug) PRINT *, 'VERTICAL UNSTAGGERING'
                  data3(:,:,1:ibtg-1,:) = (data1(:,:,1:ibtg-1,:) + data1(:,:,2:ibtg,:)) * .5
                  CALL interp (data2, data3, pres_field, interp_levels, psfc, ter, tk, qv,   &
                                 dims_in(1), dims_in(2), ibtg-1, dims_in(4), &
                                 num_metgrid_levels, LINLOG, extrapolate, .FALSE., MISSING)
                  deallocate(data3)
               ELSE
                  CALL interp (data2, data1, pres_field, interp_levels, psfc, ter, tk, qv,   &
                                 dims_in(1), dims_in(2), dims_in(3), dims_in(4), &
                                 num_metgrid_levels, LINLOG, extrapolate, .FALSE., MISSING)
               END IF

               IF (debug) write(6,*) '     SAMPLE VALUE IN  = ',data1(dims_in(1)/2,dims_in(2)/2,1,1)
               IF (debug) write(6,*) '     SAMPLE VALUE OUT = ',data2(dims_out(1)/2,dims_out(2)/2,1,1)

            ELSEIF (idm == 3 .AND. unstagger_grid ) THEN  
               IF ( dims_in(1) == iweg ) THEN
                  data2(1:iweg-1,:,:,:) = (data1(1:iweg-1,:,:,:) + data1(2:iweg,:,:,:)) * .5
               ELSEIF ( dims_in(2) == isng ) THEN
                  data2(:,1:isng-1,:,:) = (data1(:,1:isng-1,:,:) + data1(:,2:isng,:,:)) * .5
               ELSE
                  data2 = data1
               ENDIF
               IF (debug) write(6,*) '     SAMPLE VALUE  = ',data1(dims_in(1)/2,dims_in(2)/2,1,1)

            ELSE
                  data2 = data1
               IF (debug) write(6,*) '     SAMPLE VALUE  = ',data1(dims_in(1)/2,dims_in(2)/2,1,1)

            ENDIF
            
            rcode = nf_put_vara_real (mcid, jvar, start_dims, dims_out, data2)
            !print *, 'mcid ', mcid
            !print *, 'jvar ', jvar
            !print *, 'start_dims ', start_dims
            !print *, 'dims_out ', dims_out

            deallocate (data1)
            deallocate (data2)
  
          ENDIF
  
        ENDDO loop_variables

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! ADDITIONAL VARIABLES !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF ( debug ) print*," "
IF ( debug ) print*,"Calculating some diagnostics"
interpolate = .TRUE.

                              !!! NB: what follows is not useful because we'd like diagnostics for each history timestep
                                      jshape = 0
                                      DO ii = 1, 4
                                        IF ( ii == 1 ) test_dim_name = 'west_east'
                                        IF ( ii == 2 ) test_dim_name = 'south_north'
                                        IF (( ii == 3 ) .and. (LINLOG .le. 2))  test_dim_name = 'pressure'
                                        IF (( ii == 3 ) .and. (LINLOG .eq. 3))  test_dim_name = 'altitude' !'bottom_top' !'altitude'
                                        IF (( ii == 3 ) .and. (LINLOG .eq. 4))  test_dim_name = 'altitude_abg' !'bottom_top' !'altitude_abg'
                                        IF ( ii == 4 ) test_dim_name = 'Time'
                                        DO jj = 1,j
                                          IF ( test_dim_name == dnamej(jj) ) THEN
                                            jshape(ii) = jj
                                          ENDIF
                                        ENDDO
                                        IF ( jshape(ii) == 0 ) THEN
                                          j = j + 1
                                          jshape(ii) = j
                                          dnamej(j) = dnamei(ishape(ii))
                                          dvalj(j) = dvali(ishape(ii))
                                          rcode = nf_def_dim(mcid, dnamej(j), dvalj(j), j)
                                        ENDIF
                                      ENDDO
                                      dims_in  = 1
                                      dims_out = 1
                                      DO ii = 1,4
                                        dims_out(ii) = dvalj(jshape(ii))
                                        !print *, dims_out(ii)
                                      ENDDO
                                      !!! NB: what follows is useful because we'd like diagnostics for each history timestep
                                      dims_in = dims_out

        !IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'PRES') /= 0 ) THEN
        !   !!! PRES
        !   jvar = jvar + 1
        !   CALL def_var (mcid, jvar, "PRES", 5, 4, jshape, "XZY", "Pressure           ", "Pa", "-", "XLONG XLAT")
        !   IF (debug) THEN
        !     write(6,*) 'VAR: PRES'
        !     write(6,*) '     DIMS OUT: ',dims_out
        !   ENDIF
        !   rcode = nf_put_vara_real (mcid, jvar, start_dims, dims_out, pres_out)
        !   IF (debug) write(6,*) '     SAMPLE VALUE OUT = ',pres_out(dims_out(1)/2,dims_out(2)/2,1,1)
        !ENDIF

!
! EN FAIT IL FAUDRAIT QUE LE jshape SOIT BIEN DEFINI POUR CE QUI SUIT
!

        !
        ! OUTPUT REGULAR TEMPERATURE 
        !
        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'tk') /= 0 ) THEN
           jvar = jvar + 1
           CALL def_var (mcid, jvar, "tk  ", 5, 4, jshape, "XZY", "Temperature         ", "K ", "-", "XLONG XLAT", MISSING)
           IF (debug) THEN
             write(6,*) 'VAR: tk'
             write(6,*) '     DIMS OUT: ',dims_out
           ENDIF
           allocate (data2(dims_out(1),dims_out(2),dims_out(3),dims_out(4)))
           CALL interp (data2, tk, pres_field, interp_levels, psfc, ter, tk, qv,   &
                          iweg-1, isng-1, ibtg-1, dims_in(4), &          
                          num_metgrid_levels, LINLOG, extrapolate, .FALSE., MISSING)
           rcode = nf_put_vara_real (mcid, jvar, start_dims, dims_out, data2)
           IF (debug) write(6,*) '     SAMPLE VALUE OUT = ',data2(dims_out(1)/2,dims_out(2)/2,1,1)
           deallocate(data2)
        ENDIF

        !
        ! OUTPUT ROTATED WINDS
        !
        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'uvmet') /= 0 ) THEN
           jvar = jvar + 1
           CALL def_var (mcid, jvar, "Um  ", 5, 4, jshape, "XZY", "U rotated wind      ", "K ", "-", "XLONG XLAT", MISSING)
           IF (debug) THEN
             write(6,*) 'VAR: u (rotated)'
             write(6,*) '     DIMS OUT: ',dims_out
           ENDIF
           allocate (data2(dims_out(1),dims_out(2),dims_out(3),dims_out(4)))
           CALL interp (data2, umet, pres_field, interp_levels, psfc, ter, tk, qv,  &
                          iweg-1, isng-1, ibtg-1, dims_in(4), &
                          num_metgrid_levels, LINLOG, extrapolate, .FALSE., MISSING)
           rcode = nf_put_vara_real (mcid, jvar, start_dims, dims_out, data2)
           IF (debug) write(6,*) '     SAMPLE VALUE OUT=',data2(dims_out(1)/2,dims_out(2)/2,1,1)
           deallocate(data2)
           
           jvar = jvar + 1
           CALL def_var (mcid, jvar, "Vm  ", 5, 4, jshape, "XZY", "V rotated wind      ", "K ", "-", "XLONG XLAT", MISSING)
           IF (debug) THEN
             write(6,*) 'VAR: v (rotated)'
             write(6,*) '     DIMS OUT: ',dims_out
           ENDIF
           allocate (data2(dims_out(1),dims_out(2),dims_out(3),dims_out(4)))
           CALL interp (data2, vmet, pres_field, interp_levels, psfc, ter, tk, qv,  &
                          iweg-1, isng-1, ibtg-1, dims_in(4), &
                          num_metgrid_levels, LINLOG, extrapolate, .FALSE., MISSING)
           rcode = nf_put_vara_real (mcid, jvar, start_dims, dims_out, data2)
           IF (debug) write(6,*) '     SAMPLE VALUE OUT=',data2(dims_out(1)/2,dims_out(2)/2,1,1)
           deallocate(data2)
        ENDIF

        !
        ! OUTPUT POTENTIAL TEMPERATURE
        !
        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'tpot') /= 0 ) THEN
           jvar = jvar + 1
           CALL def_var (mcid, jvar, "tpot", 5, 4, jshape, "XZY", "Potential Temperature ", "K ", "-", "XLONG XLAT", MISSING)
           IF (debug) THEN
             write(6,*) 'VAR: tpot'
             write(6,*) '     DIMS OUT: ',dims_out
           ENDIF
           allocate (data2(dims_out(1),dims_out(2),dims_out(3),dims_out(4)))
           CALL interp (data2, tpot, pres_field, interp_levels, psfc, ter, tk, qv,  &
                          iweg-1, isng-1, ibtg-1, dims_in(4), &
                          num_metgrid_levels, LINLOG, extrapolate, .FALSE., MISSING)
           rcode = nf_put_vara_real (mcid, jvar, start_dims, dims_out, data2)
           IF (debug) write(6,*) '     SAMPLE VALUE OUT =',data2(dims_out(1)/2,dims_out(2)/2,1,1)
           deallocate(data2)
        ENDIF

        !
        ! OUTPUT GEOPOTENTIAL HEIGHT
        !
        IF (LINLOG .le. 2) THEN
        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'GHT') /= 0 ) THEN
           jvar = jvar + 1
           CALL def_var (mcid, jvar, "GHT ", 5, 4, jshape, "XZY", "Geopotential Height", "m ", "-", "XLONG XLAT", MISSING)
           IF (debug) THEN
             write(6,*) 'VAR: GHT'
             write(6,*) '     DIMS OUT: ',dims_out
           ENDIF
           allocate (data2(dims_out(1),dims_out(2),dims_out(3),dims_out(4)))
           CALL interp (data2, phb(:,:,1:ibtg-1,:), pres_field, interp_levels, psfc, ter, tk, qv,   &
                          iweg-1, isng-1, ibtg-1, dims_in(4), &          
                          num_metgrid_levels, LINLOG, extrapolate, .TRUE., MISSING)
           data2 = data2/grav
           rcode = nf_put_vara_real (mcid, jvar, start_dims, dims_out, data2)
           IF (debug) write(6,*) '     SAMPLE VALUE OUT = ',data2(dims_out(1)/2,dims_out(2)/2,1,1)
           deallocate(data2)
        ENDIF
        ELSE
           PRINT *, 'Geopotential height is actually vertical coordinate'     
        ENDIF

        !        IF ( INDEX(process,'all') /= 0 .OR. INDEX(process_these_fields,'RH') /= 0 ) THEN
        !           !!! RH
        !           jvar = jvar + 1
        !           CALL def_var (mcid, jvar, "RH  ", 5, 4, jshape, "XZY", "Relative Humidity  ", "% ", "-", "XLONG XLAT")
        !           IF (debug) THEN
        !             write(6,*) 'VAR: RH'
        !             write(6,*) '     DIMS OUT: ',dims_out
        !           ENDIF
        !           allocate (data2(dims_out(1),dims_out(2),dims_out(3),dims_out(4)))
        !           CALL interp (data2, rh, pres_field, interp_levels, psfc, ter, tk, qv,   &
        !                          iweg-1, isng-1, ibtg-1, dims_in(4), &          
        !                          num_metgrid_levels, LINLOG, extrapolate, .FALSE., MISSING)
        !           WHERE ( rh < 0.0 ) 
        !              rh = 0.0
        !           ENDWHERE
        !           WHERE ( rh > 100.0 )
        !              rh = 100.0
        !           ENDWHERE
        !           rcode = nf_put_vara_real (mcid, jvar, start_dims, dims_out, data2)
        !           IF (debug) write(6,*) '     SAMPLE VALUE OUT = ',data2(dims_out(1)/2,dims_out(2)/2,1,1)
        !           deallocate(data2)
        !        ENDIF


!=====================================================================================
!=====================================================================================
!
! VERTICAL COORDINATES
!
jvar = jvar + 1
jshape(:)=0
jshape(1)=2
CALL def_var (mcid, jvar, "vert", 5, 1, jshape, " Z ", "Vert. coord.        ","m ", "-", "XLONG XLAT", MISSING)
start_dims(1)=1
start_dims(2)=1
start_dims(3)=1
start_dims(4)=1
dims_out(1)=dims_out(3)
dims_out(2)=1
dims_out(3)=1
dims_out(4)=1
allocate (data2(dims_out(1),dims_out(2),dims_out(3),dims_out(4)))
do kk=1,dims_out(1)
  data2(kk,1,1,1) = interp_levels(kk)
enddo
rcode = nf_put_vara_real (mcid, jvar, start_dims, dims_out, data2)
deallocate(data2)
!=====================================================================================
!=====================================================================================



      ! 
      ! CLOSE FILES ON EXIT
      !
        rcode = nf_close(ncid)
        rcode = nf_close(mcid)
        write(6,*) 
      ENDDO 

      !!!! forgotten in the initial program
      DEALLOCATE (  input_file_names )
      DEALLOCATE ( output_file_names )

      !
      ! WELL...
      !
      write(6,*) "I used Cp, R : ",Cp,Rd 
      write(6,*)
      write(6,*) "##########################################"
      write(6,*) " END of API PROGRAM - SUCCESS so far"      
      write(6,*) "##########################################"

END SUBROUTINE
! END PROGRAM api
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
 SUBROUTINE handle_err(rcode)
    INTEGER rcode
    write(6,*) 'Error number ',rcode
    stop
 END SUBROUTINE
!---------------------------------------------------------------------
 SUBROUTINE interp (data_out, data_in, pres_field, interp_levels, psfc, ter, tk, qv, ix, iy, iz, it, &
                     num_metgrid_levels, LINLOG, extrapolate, GEOPT, MISSING)

     INTEGER                                          :: ix, iy, iz, it
     INTEGER                                          :: num_metgrid_levels, LINLOG
     REAL, DIMENSION(ix, iy, num_metgrid_levels, it)  :: data_out
     REAL, DIMENSION(ix, iy, iz, it)                  :: data_in, pres_field, tk, qv
     REAL, DIMENSION(ix, iy, it)                      :: psfc
     REAL, DIMENSION(ix, iy)                          :: ter
     REAL, DIMENSION(num_metgrid_levels)              :: interp_levels

     INTEGER                                          :: i, j, itt, k, kk, kin
     REAL, DIMENSION(num_metgrid_levels)              :: data_out1D
     REAL, DIMENSION(iz)                              :: data_in1D, pres_field1D
     INTEGER                                          :: extrapolate
     REAL                                             :: MISSING
     REAL, DIMENSION(ix, iy, num_metgrid_levels, it)  :: N
     REAL                                             :: sumA, sumN, AVE_geopt
     LOGICAL                                          :: GEOPT

     N = 1.0

     do itt = 1, it
     !PRINT *, 'TIME... ', itt, ' in ', iz, ' out ', num_metgrid_levels
        do j = 1, iy
        do i = 1, ix
           data_in1D(:)    = data_in(i,j,:,itt)
           pres_field1D(:) = pres_field(i,j,:,itt)
           IF ( LINLOG .le. 2 ) THEN
             CALL int1D (data_out1D, data_in1D, pres_field1D, interp_levels, iz, num_metgrid_levels, LINLOG, MISSING)
           ELSE  
             CALL interp_1d (data_in1D, pres_field1D, iz, data_out1D, interp_levels, num_metgrid_levels, 'z', MISSING)
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            !interp_1d (data_in1D, pres_field1D, num_metgrid_levels, data_out1D, interp_levels, iz, 'z')
                            !CALL interp_1d( data_in_1d, z_data_1d, bottom_top_dim, data_out_1d, z_levs, number_of_zlevs, vertical_type)
                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ENDIF
           data_out(i,j,:,itt) = data_out1D(:)
        end do
        end do
     end do
     !PRINT *, 'ok'

     ! Fill in missing values
     IF ( extrapolate == 0 ) RETURN       !! no extrapolation - we are out of here



IF ( LINLOG .ge. 0 ) RETURN   !! TEMPORARY: no extrapolation 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CAUTION BELOW: METHOD NOT ADAPTED TO MARS 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !!! MARS MARS MARS
     !expon=287.04*.0065/9.81
     expon=192.*.0045/3.72

     ! First find where about 400 hPa is located
     kk = 0
     find_kk : do k = 1, num_metgrid_levels
        kk = k
        if ( interp_levels(k) <= 40000. ) exit find_kk
     end do find_kk

     
     IF ( GEOPT ) THEN     !! geopt is treated different below ground

        do itt = 1, it
           do k = 1, kk
              do j = 1, iy
              do i = 1, ix
                 IF ( data_out(i,j,k,itt) == MISSING .AND. interp_levels(k) < psfc(i,j,itt) ) THEN

!                We are below the first model level, but above the ground 

!!! MARS MARS
                    data_out(i,j,k,itt) = ((interp_levels(k) - pres_field(i,j,1,itt))*ter(i,j)*3.72 +  &
                                           (psfc(i,j,itt) - interp_levels(k))*data_in(i,j,1,itt) ) /   &
                                          (psfc(i,j,itt) - pres_field(i,j,1,itt))

                 ELSEIF ( data_out(i,j,k,itt) == MISSING ) THEN

!                We are below both the ground and the lowest data level.

!                First, find the model level that is closest to a "target" pressure
!                level, where the "target" pressure is delta-p less that the local
!                value of a horizontally smoothed surface pressure field.  We use
!                delta-p = 150 hPa here. A standard lapse rate temperature profile
!                passing through the temperature at this model level will be used
!                to define the temperature profile below ground.  This is similar
!                to the Benjamin and Miller (1990) method, except that for
!                simplicity, they used 700 hPa everywhere for the "target" pressure.
!                Code similar to what is implemented in RIP4


!!! MARS MARS MARS
                    ptarget = (psfc(i,j,itt)*.01) - 1.50 !150.
                    dpmin=1.e2 !1.e4
                    kupper = 0
                    loop_kIN : do kin=iz,1,-1
                       kupper = kin
                       dp=abs( (pres_field(i,j,kin,itt)*.01) - ptarget )
                       if (dp.gt.dpmin) exit loop_kIN
                       dpmin=min(dpmin,dp)
                    enddo loop_kIN

                    pbot=max(pres_field(i,j,1,itt),psfc(i,j,itt))
                    zbot=min(data_in(i,j,1,itt)/3.72,ter(i,j))   !!MARS MARS

                    tbotextrap=tk(i,j,kupper,itt)*(pbot/pres_field(i,j,kupper,itt))**expon
                    tvbotextrap=virtual(tbotextrap,qv(i,j,1,itt))

                    data_out(i,j,k,itt) = (zbot+tvbotextrap/.0045*(1.-(interp_levels(k)/pbot)**expon))*3.72   !!MARS MARS
               
                 ENDIF
              enddo
              enddo
           enddo
        enddo


        !!! Code for filling missing data with an average - we don't want to do this
        !!do itt = 1, it
           !!loop_levels : do k = 1, num_metgrid_levels
              !!sumA = SUM(data_out(:,:,k,itt), MASK = data_out(:,:,k,itt) /= MISSING)
              !!sumN = SUM(N(:,:,k,itt), MASK = data_out(:,:,k,itt) /= MISSING)
              !!IF ( sumN == 0. ) CYCLE loop_levels
              !!AVE_geopt = sumA/sumN
              !!WHERE ( data_out(:,:,k,itt) == MISSING )
                 !!data_out(:,:,k,itt) = AVE_geopt
              !!END WHERE
           !!end do loop_levels
        !!end do

     END IF
     
     !!! All other fields and geopt at higher levels come here
     do itt = 1, it
        do j = 1, iy
        do i = 1, ix
          do k = 1, kk
             if ( data_out(i,j,k,itt) == MISSING ) data_out(i,j,k,itt) = data_in(i,j,1,itt)
          end do
          do k = kk+1, num_metgrid_levels
             if ( data_out(i,j,k,itt) == MISSING ) data_out(i,j,k,itt) = data_in(i,j,iz,itt)
          end do
        end do
        end do
     end do
       

 END SUBROUTINE interp 
!------------------------------------------------------------------------------
!--------------------------------------------------------
 SUBROUTINE int1D(xxout, xxin, ppin, ppout, npin, npout, LINLOG, MISSING)

! Modified from int2p - NCL code
! routine to interpolate from one set of pressure levels
! .   to another set  using linear or ln(p) interpolation
!
! NCL: xout = int2p (pin,xin,pout,linlog)
! This code was originally written for a specific purpose.
! .   Several features were added for incorporation into NCL's
! .   function suite including linear extrapolation.
!
! nomenclature:
!
! .   ppin   - input pressure levels. The pin can be
! .            be in ascending or descending order
! .   xxin   - data at corresponding input pressure levels
! .   npin   - number of input pressure levels >= 2
! .   ppout  - output pressure levels (input by user)
! .            same (ascending or descending) order as pin
! .   xxout  - data at corresponding output pressure levels
! .   npout  - number of output pressure levels
! .   linlog - if abs(linlog)=1 use linear interp in pressure
! .            if abs(linlog)=2 linear interp in ln(pressure)
! .   missing- missing data code. 

!                                                ! input types
      INTEGER   :: npin,npout,linlog,ier
      real      :: ppin(npin),xxin(npin),ppout(npout)
      real      :: MISSING       
     logical                                          :: AVERAGE
!                                                ! output
      real      :: xxout(npout)
      INTEGER   :: j1,np,nl,nin,nlmax,nplvl
      INTEGER   :: nlsave,np1,no1,n1,n2,nlstrt
      real      :: slope,pa,pb,pc

! automatic arrays
      real      :: pin(npin),xin(npin),p(npin),x(npin)
      real      :: pout(npout),xout(npout)


      xxout = MISSING
      pout  = ppout
      p     = ppin
      x     = xxin
      nlmax = npin

! exact p-level matches
      nlstrt = 1
      nlsave = 1
      do np = 1,npout
          xout(np) = MISSING
          do nl = nlstrt,nlmax
              if (pout(np).eq.p(nl)) then
                  xout(np) = x(nl)
                  nlsave = nl + 1
                  go to 10
              end if
          end do
   10     nlstrt = nlsave
      end do

      if (LINLOG.eq.1) then
          do np = 1,npout
              do nl = 1,nlmax - 1
                  if (pout(np).lt.p(nl) .and. pout(np).gt.p(nl+1)) then
                      slope = (x(nl)-x(nl+1))/ (p(nl)-p(nl+1))
                      xout(np) = x(nl+1) + slope* (pout(np)-p(nl+1))
                  end if
              end do
          end do
      elseif (LINLOG.eq.2) then
          do np = 1,npout
              do nl = 1,nlmax - 1
                  if (pout(np).lt.p(nl) .and. pout(np).gt.p(nl+1)) then
                      pa = log(p(nl))
                      pb = log(pout(np))
! special case: in case someone inadvertently enter p=0.
                      if (p(nl+1).gt.0.d0) then
                          pc = log(p(nl+1))
                      else
                          pc = log(1.d-4)
                      end if

                      slope = (x(nl)-x(nl+1))/ (pa-pc)
                      xout(np) = x(nl+1) + slope* (pb-pc)
                  end if
              end do
          end do
      end if


! place results in the return array;
      xxout = xout

 END SUBROUTINE int1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE interp_1d( a, xa, na, b, xb, nb, vertical_type, MISSING)

  implicit none

 ! Arguments
  integer, intent(in)              :: na, nb
  real, intent(in), dimension(na)  :: a, xa
  real, intent(in), dimension(nb)  :: xb
  real, intent(out), dimension(nb) :: b
  character (len=1)                :: vertical_type

  !! Local variables
  real                             :: MISSING                                
  integer                          :: n_in, n_out
  real                             :: w1, w2
  logical                          :: interp

  IF ( vertical_type == 'p' ) THEN

    DO n_out = 1, nb

      b(n_out) = MISSING 
      interp = .false.
      n_in = 1

      DO WHILE ( (.not.interp) .and. (n_in < na) )
        IF( (xa(n_in)   >= xb(n_out)) .and. &
            (xa(n_in+1) <= xb(n_out))        ) THEN
          interp = .true.
          w1 = (xa(n_in+1)-xb(n_out))/(xa(n_in+1)-xa(n_in))
          w2 = 1. - w1
          b(n_out) = w1*a(n_in) + w2*a(n_in+1)
        END IF
        n_in = n_in +1
      ENDDO

    ENDDO

  ELSE

    DO n_out = 1, nb
  
      b(n_out) = MISSING 
      interp = .false.
      n_in = 1
  
      DO WHILE ( (.not.interp) .and. (n_in < na) )
        IF( (xa(n_in)   <= xb(n_out)) .and. &
            (xa(n_in+1) >= xb(n_out))        ) THEN
          interp = .true.
          w1 = (xa(n_in+1)-xb(n_out))/(xa(n_in+1)-xa(n_in))
          w2 = 1. - w1
          b(n_out) = w1*a(n_in) + w2*a(n_in+1)
        END IF
        n_in = n_in +1
      ENDDO

    ENDDO

  END IF

  END SUBROUTINE interp_1d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------
 FUNCTION virtual (tmp,rmix)
!      This function returns virtual temperature in K, given temperature
!      in K and mixing ratio in kg/kg.

     real                              :: tmp, rmix, virtual

     virtual=tmp*(0.622+rmix)/(0.622*(1.+rmix))

 END FUNCTION virtual

!------------------------------------------------------------------------------
 SUBROUTINE all_spaces ( command , length_of_char ) 

      IMPLICIT NONE

      INTEGER                                       :: length_of_char
      CHARACTER (LEN=length_of_char)                :: command
      INTEGER                                       :: loop

      DO loop = 1 , length_of_char
         command(loop:loop) = ' '
      END DO

 END SUBROUTINE all_spaces

!------------------------------------------------------------------------------
 SUBROUTINE def_var (mcid, jvar, cval, itype, idm, jshape, order, desc, units, stag, coord, missing )

      IMPLICIT NONE

      INCLUDE 'netcdf.inc'

      INTEGER              :: mcid, jvar
      CHARACTER (LEN =  4) :: cval
      INTEGER              :: itype, idm
      REAL, DIMENSION(6)   :: jshape
      CHARACTER (LEN =  3) :: order
      CHARACTER (LEN = 19) :: desc
      CHARACTER (LEN =  2) :: units
      CHARACTER (LEN =  1) :: stag
      CHARACTER (LEN = 10) :: coord

      INTEGER            :: rcode, ilen
      CHARACTER (LEN=30) :: att_text

      REAL               :: missing


      IF ( itype == 5 ) THEN
         rcode = nf_redef(mcid)
         rcode = nf_def_var(mcid, trim(cval), NF_REAL, idm, jshape, jvar)
         rcode = nf_put_att_int(mcid, jvar, "FieldType", NF_INT, 1, 104)
      ENDIF

      att_text = order
      ilen = len_trim(att_text)
      rcode = nf_put_att_text(mcid, jvar, "MemoryOrder", ilen, att_text(1:ilen) )

      att_text = desc
      ilen = len_trim(att_text)
      rcode = nf_put_att_text(mcid, jvar, "description", ilen, att_text(1:ilen) )

      att_text = units
      ilen = len_trim(att_text)
      rcode = nf_put_att_text(mcid, jvar, "units", ilen, att_text(1:ilen) )

      att_text = stag
      ilen = len_trim(att_text)
      rcode = nf_put_att_text(mcid, jvar, "stagger", ilen, att_text(1:ilen) )

      att_text = coord
      ilen = len_trim(att_text)
      rcode = nf_put_att_text(mcid, jvar, "coordinates", ilen, att_text(1:ilen) )

      rcode = nf_put_att_real(mcid, jvar, "missing_value", NF_FLOAT, 1, MISSING )

      rcode = nf_enddef(mcid)

 END SUBROUTINE def_var
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!! Diagnostics: U & V on earth coordinates
  SUBROUTINE calc_uvmet( UUU,           VVV,             &
                         UUUmet,        VVVmet,          &
                         truelat1,      truelat2,        &
                         stand_lon,     map_proj,        &
                         longi,         lati,            &
                         west_east_dim, south_north_dim, bottom_top_dim)

  IMPLICIT NONE

  !Arguments
  integer        :: west_east_dim, south_north_dim, bottom_top_dim
  real, dimension(west_east_dim,south_north_dim,bottom_top_dim)    :: UUU
  real, dimension(west_east_dim,south_north_dim,bottom_top_dim)    :: VVV
  real, dimension(west_east_dim,south_north_dim,bottom_top_dim)    :: UUUmet
  real, dimension(west_east_dim,south_north_dim,bottom_top_dim)    :: VVVmet
  real, dimension(west_east_dim,south_north_dim)                   :: longi, lati
  real :: truelat1, truelat2, stand_lon
  integer :: map_proj

  !Local
  integer                                         :: i, j, k
  real                                            :: cone
  real, dimension(west_east_dim,south_north_dim)  :: diff, alpha

   real, parameter :: PI = 3.141592653589793
   real, parameter :: DEG_PER_RAD = 180./PI
   real, parameter :: RAD_PER_DEG = PI/180.

  IF ( map_proj .ge. 3 ) THEN                         ! No need to rotate
    !PRINT *, 'NO NEED TO ROTATE !!!! equivalent to output U,V with unstagger_grid'
    UUUmet(:,:,:) = UUU
    VVVmet(:,:,:) = VVV
  ELSE

  cone = 1.                                          !  PS
  IF ( map_proj .eq. 1) THEN                         !  Lambert Conformal mapping
    IF (ABS(truelat1-truelat2) .GT. 0.1) THEN
       cone=(ALOG(COS(truelat1*RAD_PER_DEG))-            &
             ALOG(COS(truelat2*RAD_PER_DEG))) /          &
       (ALOG(TAN((90.-ABS(truelat1))*RAD_PER_DEG*0.5 ))- &
        ALOG(TAN((90.-ABS(truelat2))*RAD_PER_DEG*0.5 )) )
    ELSE
       cone = SIN(ABS(truelat1)*RAD_PER_DEG )
    ENDIF
  END IF


  diff = longi - stand_lon
  DO i = 1, west_east_dim
  DO j = 1, south_north_dim
    IF ( diff(i,j) .gt. 180. ) THEN
      diff(i,j) = diff(i,j) - 360.
    END IF
    IF ( diff(i,j) .lt. -180. ) THEN
      diff(i,j) = diff(i,j) + 360.
    END IF
  END DO
  END DO


  DO i = 1, west_east_dim
  DO j = 1, south_north_dim
     IF ( lati(i,j) .lt. 0. ) THEN
       alpha(i,j) = - diff(i,j) * cone * RAD_PER_DEG
     ELSE
       alpha(i,j) = diff(i,j) * cone * RAD_PER_DEG
     END IF
  END DO
  END DO


  DO k = 1,bottom_top_dim
    UUUmet(:,:,k) = VVV(:,:,k)*sin(alpha) + UUU(:,:,k)*cos(alpha)
    VVVmet(:,:,k) = VVV(:,:,k)*cos(alpha) - UUU(:,:,k)*sin(alpha)
  END DO
  END IF

  END SUBROUTINE calc_uvmet

