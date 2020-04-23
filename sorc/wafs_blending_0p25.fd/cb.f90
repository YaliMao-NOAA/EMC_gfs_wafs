module cb
! ABSTRACT: This program reads 2D UPP fields and converts to CB in
!           GRIB2
!
! PROGRAM HISTORY LOG:
! 2020-04-21  Y Mao
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
! values of product def template needed to read data from GRIB 2 file
  type pdt_t
     integer :: npdt   ! number of template 4
     integer :: icat   ! catogory
     integer :: iprm   ! parameter
     integer :: ilev   ! type of level (code table 4.5)
  end type pdt_t
! PDT parameters in the input GRIB2 file (template 4 number, category, parameter, type of level)
  type(pdt_t), parameter :: &
       pdt_conv_pres_bot = pdt_t(0, 3, 0, 242), &
       pdt_conv_pres_top = pdt_t(0, 3, 0, 243), &
       pdt_conv_pcp_rate = pdt_t(8, 1, 196, 1)

! values used to write data to GRIB 2 file
  type gparms_t
     integer :: npdt   ! number of template 4
     integer :: icat   ! catogory
     integer :: iprm   ! parameter
     integer :: ilev   ! type of level (code table 4.5)
     integer :: stat   ! TYPE OF STATISTICAL PROCESSING
     !
     integer :: ndrt   ! number of template 5
     integer :: drt2   ! Binary scale factor
     integer :: drt3   ! Decimal scale factor
     integer :: drt4   ! Number of bits to hold data
     !
     real(kind=r_kind) :: msng   ! missing data (below surface)
     logical :: bitmap           ! whether to use bitmap for sparse data
  end type gparms_t

  type(gparms_t), parameter :: &
       cbcov_gparms = gparms_t(0,6,25,10,-1,40,0,3,10,"yes"),&
       cbbot_gparms = gparms_t(0,3,3,11,-1,40,13,5,16,"yes"),&
       cbtop_gparms = gparms_t(0,3,3,12,-1,40,14,5,16,"yes")

contains

!----------------------------------------------------------------------------
  subroutine cb_algo(gfile1, gfile2)
! reads input data
    implicit none
    character(*), intent(in) :: gfile1,gfile2

    integer :: ifl1,ifl2

    integer :: iret, nxy
    type(gribfield) :: gfld,gfld1,gfld2
    real, dimension(:,:), allocatable :: cbtop, cbbot, cbcov
    logical*1,allocatable,target :: bmap(:)

    call getlun90(ifl1,1)
    call getlun90(ifl2,1)

    print *, "CB file handles=",ifl1,ifl2

    ! no bit-map
    call get_grib2(ifl1, pdt_conv_pcp_rate, 0, gfld, nxy, iret)
    allocate(cbcov(nxy))
    cbcov = gfld%fld
    call cb_cover(cbcov)

    ! with bit-map and UPP CLDRAD.f -50000.
    call get_grib2(ifl1, pdt_conv_pres_bot, 0, gfld1, nxy, iret)
    allocate(cbbot(nxy))
    cbbot = gfld%fld

    ! with bit-map and UPP CLDRAD.f -50000.
    call get_grib2(ifl1, pdt_conv_pres_top, 0, gfld2, nxy, iret)
    allocate(cbtop(nxy))
    cbtop = gfld%fld

    allocate(bmap(nxy))
    bmap = .true.

    where(cbcov <= 0.0)
       bmap = .false.
    end where

    where(.not. gfld1%ibmap)
       bmap = .false.
    elsewhere(cbbot <= 0.)
       bmap = .false.
    end where

    where(.not. gfld2%ibmap)
       bmap = .false.
    elsewhere(cbtop <= 0.)
       bmap = .false.
    end where

    where(bmap)
       where(.not. (cbtop < 400.*100. .and. &
             cbbot - cbtop > 300.*100.))
          bmap = .false.
       else
          cbbot = P2H(cbbot)
          cbtop = P2H(cbtop)
       end where
    end where

    call put_grib2(ifl2,cbcov_gparms,0,gfld, 0,bmap,cbcov,iret) 
    call put_grib2(ifl2,cbbot_gparms,0,gfld1,0,bmap,cbbot,iret) 
    call put_grib2(ifl2,cbtop_gparms,0,gfld2,0,bmap,cbtop,iret) 

    deallocate(bmap)
    deallocate(cbbot)
    deallocate(cbtop)
    deallocate(cbcov)
  end subroutine cb_algo

  SUBROUTINE GETLUN90(LUN,OPTN)
!* THIS PROGRAM GETS UNIQUE LOGICAL UNIT NUMBERS FOR OPFILE
!* OR RETURNS THEM TO THE POOL FOR CLFILE
    IMPLICIT NONE
    INTEGER, PARAMETER :: CNCT=1,DSCT=2
    INTEGER :: LUN,OPTN,I
    INTEGER :: NUM(80)=(/ &
                  99,98,97,96,95,94,93,92,91,90, &
                  89,88,87,86,85,84,83,82,81,80, &
                  79,78,77,76,75,74,73,72,71,70, &
                  69,68,67,66,65,64,63,62,61,60, &
                  59,58,57,56,55,54,53,52,51,50, &
                  49,48,47,46,45,44,43,42,41,40, &
                  39,38,37,36,35,34,33,32,31,30, &
                  29,28,27,26,25,24,23,22,21,20 /)
!* START
    IF(OPTN == CNCT) THEN
       DO I=1,80
          IF(NUM(I)>0) THEN
             LUN=NUM(I)
             NUM(I)=-NUM(I)
             return
          ENDIF
       END DO
       PRINT*, 'NEED MORE THAN 80 UNIT NUMBERS'
    ELSE IF(OPTN == DSCT) THEN
!* MAKE THE NUMBER AVAILABLE BY SETTING POSITIVE
       DO I=1,80
          IF(LUN == -NUM(I)) NUM(I)=ABS(NUM(I))
       ENDDO
    END IF

    RETURN
  END SUBROUTINE GETLUN90

  ELEMENTAL FUNCTION P2H(p)
    implicit none
    real, intent(in) :: p
    real :: P2H
!   To convert pressure levels (Pa) to geopotantial heights
!   Uses ICAO standard atmosphere parameters as defined here:
!      https://www.nen.nl/pdfpreview/preview_29424.pdf
    real, parameter :: lapse = 0.0065
    real, parameter :: surf_temp = 288.15
    real, parameter :: gravity = 9.80665
    real, parameter :: moles_dry_air = 0.02896442
    real, parameter :: gas_const = 8.31432
    real, parameter :: surf_pres = 1.01325e5
    real, parameter :: power_const = (gravity * moles_dry_air) &
                                       / (gas_const * lapse)
    real, parameter :: strat_temp=216.65
    real, parameter :: strat_pres=22631.7
    real, parameter :: con_rd = 2.8705e+2      ! gas constant air
    real, parameter :: alpha_strat = -con_rd*strat_temp/gravity
    if(p >= strat_pres) then
       P2H = (surf_temp/lapse)*(1-(p/surf_pres)**(1/power_const))
    else
       P2H = 11000. + alpha_strat * log(p/strat_pres)
    end if
  END FUNCTION P2H

  ! Calculate CB coverage by using fuzzy logic
  ! Evaluate membership of val in a fuzzy set fuzzy.
  ! Assume f is in x-log scale
  subroutine cb_cover(cbcov)
    implicit none
    real, intent(inout) :: cbcov

    ! x - convective precipitation [1.0e6*kg/(m2s)]
    ! y - cloud cover fraction, between 0 and 1
    ! These are original values from Slingo (Table 1):
    ! c = -.006 + 0.125*log(p)
    ! x = 1.6 3.6 8.1 18.5 39.0 89.0 197.0 440.0 984.0 10000.0
    ! y = 0.0 0.1 0.2  0.3  0.4  0.5   0.6   0.7   0.8     0.8
    integer, parameter :: NP=10
    real, parameter :: x(NP) = &
         (/ 1.6,3.6,8.1,18.5,39.0,89.0,197.0,440.0,984.0,10000.0 /)   
    real, parameter :: y(NP) = &
         (/ 0.0,0.1,0.2, 0.3, 0.4, 0.5,  0.6,  0.7,  0.8,    0.8 /)
    real, xlog(NP)

    xlog = log(x)

    where(cbcov <= 0. )
       cbcov = 0.0
    elsewhere
       cbcov = fuzzy_member(NP,xlog,y,log(1.0e6*cbcov))
    end where
  end subroutine cb_cover

  elemental function fuzzy_member(NP,x,y,val)
    real :: fuzzy_member
    integer,intent(in)::NP
    real,intent(in) :: x(NP)
    real,intent(in) :: y(NP)
    real,intent(in) :: val

    integer :: i
    real :: delta, mem

    if (val <= x(1)) then
        mem = 0.0
    else if (val >= x(NP)) then
        mem = 0.0
    else
        do i = 2, NP
            if (val < x(i)) then
                delta = x(i) -  x(i-1)
                if (delta <= 0.0) then
                    mem = y(i-1)
                else
                    mem = y(i) * (val-x(i-1)) + &
                          y(i-1) * (x(i)-val)) / delta
                end if
                exit
            end if
        end do
    end if
    fuzzy_member = mem
  end function fuzzy_member

!----------------------------------------------------------------------------
  subroutine get_grib2(iunit,pdt, pres_level, gfld, nxy, iret)
    implicit none
    integer, intent(in) :: iunit
    type(pdt_t), intent(in) :: pdt
    integer, intent(in) :: pres_level ! pressure level in Pa
    type(gribfield), intent(out) :: gfld
    integer, intent(out) :: nxy
    integer, intent(out) :: iret

    integer j,jdisc,jpdtn,jgdtn
    integer,dimension(200) :: jids,jpdt,jgdt
    logical :: unpack

    integer :: i

    j        = 0          ! search from 0
    jdisc    = 0          ! for met field:0 hydro: 1, land: 2
    jids(:)  = -9999
    !-- set product defination template 4
    jpdtn    = pdt%npdt   ! number of product defination template 4
    jpdt(:)  = -9999
    jpdt(1)  = pdt%icat   ! category 
    jpdt(2)  = pdt%iprm   ! parameter number
    jpdt(10) = pdt%ilev   ! type of level (code table 4.5)
    jpdt(12) = pres_level ! level value
    !-- set grid defination template/section 3
    jgdtn    = -1  
    jgdt(:)  = -9999
    unpack=.true.
    ! Get field from file
    if(jpdtn == 8) then
       do i = 1, 6 ! Bucket precip accumulation time up to 6 hour
          jpdt(27) = i
          call getgb2(glob_lu_in, 0, j, jdisc, jids, jpdtn, jpdt, &
               jgdtn, jgdt, unpack, j, gfld, iret)
          if( iret == 0) then
             print *,'call get_grib2, iret=',iret, pdt,"at bucket accumulation time=",i
             exit
          endif
       end do
    else
       call getgb2(glob_lu_in, 0, j, jdisc, jids, jpdtn, jpdt, &
            jgdtn, jgdt, unpack, j, gfld, iret)
       if( iret /= 0) then
          print *,'call get_grib2, iret=',iret, pdt,"on level=",pres_level 
       endif
    end if
    nxy = gfld%igdtmpl(8) * gfld%igdtmpl(9)

  end subroutine get_grib2

!----------------------------------------------------------------------------
  subroutine put_grib2(ifl,parms, nlevel, gfld, ibmap,bmap,fld, iret)
! basically the same as putgb2, but with flexible template 4 and template 5
! writes calculated values for one field at all pressure levels
    implicit none
    integer, intent(in) :: ifl
    type(gparms_t), intent(in) :: parms    ! grib2 parameters of template 4 & 5
    integer, intent(in) :: nlevel          ! pressure level in Pa, integer
    type(gribfield), intent(in) :: gfld    ! a sample input carrying information
    integer, intent(in) :: ibmap ! indicator whether to use bitmap
    logical, intent(in) :: bmap(:)
    real(4), intent(in) :: fld(:)     ! the data to be written
    integer, intent(out) :: iret           ! return status code  

    CHARACTER(LEN=1),ALLOCATABLE,DIMENSION(:) :: CGRIB
    integer(4) :: lcgrib, lengrib
    integer :: listsec0(2)
    integer :: igds(5)
    real    :: coordlist=0.0
    integer :: ilistopt=0
    ! flexible arrays of template 4, 5
    integer, allocatable :: ipdtmpl(:), idrtmpl(:)

    character(len=*), parameter :: myself = 'put_grib2(): '

!   ALLOCATE ARRAY FOR GRIB2 FIELD
    lcgrib=gfld%ngrdpts*4
    allocate(cgrib(lcgrib),stat=iret)
    if ( iret/=0 ) then
       print *, myself, iret
       iret=2
    endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  CREATE NEW MESSAGE
    listsec0(1)=gfld%discipline
    listsec0(2)=gfld%version
    if ( associated(gfld%idsect) ) then
       call gribcreate(cgrib,lcgrib,listsec0,gfld%idsect,iret)
       if (iret /= 0) then
          write(*,*) myself, ' ERROR creating new GRIB2 field = ',iret
       endif
    else
       print *, myself, ' No Section 1 info available. '
       iret=10
       deallocate(cgrib)
       return
    endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  ADD GRID TO GRIB2 MESSAGE (Grid Definition Section 3)
    igds(1)=gfld%griddef    ! Source of grid definition (see Code Table 3.0)
    igds(2)=gfld%ngrdpts    ! Number of grid points in the defined grid.
    igds(3)=gfld%numoct_opt ! Number of octets needed for each additional grid points definition
    igds(4)=gfld%interp_opt ! Interpretation of list for optional points definition (Code Table 3.11)
    igds(5)=gfld%igdtnum    ! Grid Definition Template Number (Code Table3.1)
    if ( associated(gfld%igdtmpl) ) then
       call addgrid(cgrib, lcgrib, igds, gfld%igdtmpl, gfld%igdtlen,&
                   ilistopt, gfld%num_opt, iret)
       if (iret/=0) then
          write(*,*) myself, ' ERROR adding grid info = ',iret
       endif
    else
       print *, myself, ' No GDT info available. '
       iret=11
       deallocate(cgrib)
       return
    endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  ADD DATA FIELD TO GRIB2 MESSAGE
    ! template 4
    if( parms%npdt == 0 .or. parms%npdt == 7) then
       allocate(ipdtmpl(15))
    else
       allocate(ipdtmpl(18))
       ipdtmpl(16) = parms%stat
       ipdtmpl(17) = 3
       ipdtmpl(18) = 1
    endif
    ipdtmpl(1:15) = gfld%ipdtmpl(1:15)
    ipdtmpl(1)    = parms%icat
    ipdtmpl(2)    = parms%iprm
    ipdtmpl(10)   = parms%ilev
    ipdtmpl(12)   = nlevel
    ! template 5
    if( parms%ndrt == 40) then
       allocate(idrtmpl(7))
    endif
    idrtmpl(1) = 0 ! Any value. Will be overwritten
    idrtmpl(2) = parms%drt2
    idrtmpl(3) = parms%drt3
    idrtmpl(4) = parms%drt4
    idrtmpl(5) = 0
    idrtmpl(6) = 0
    idrtmpl(7) = 255
    ! call addfield
    call addfield(cgrib, lcgrib, parms%npdt, ipdtmpl, & 
                  size(ipdtmpl), coordlist, gfld%num_coord, &
                  parms%ndrt, idrtmpl, size(idrtmpl), &
                  fld, gfld%ngrdpts, ibmap, bmap, iret)
    if (iret /= 0) then
       write(*,*) myself, 'ERROR adding data field = ',iret
    endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  CLOSE GRIB2 MESSAGE AND WRITE TO FILE
    call gribend(cgrib, lcgrib, lengrib, iret)
    call wryte(ifl, lengrib, cgrib)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    deallocate(cgrib)
    deallocate(ipdtmpl)
    deallocate(idrtmpl)
    RETURN
  end subroutine put_grib2

end module cb

program main

  use cb

  implicit none

  character(60) :: gfile1,gfile2

  INTEGER :: NARG

  call start()

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  GET ARGUMENTS
  NARG=IARGC()
  IF(NARG /= 2) THEN
     CALL ERRMSG('cb_algo:  Incorrect usage')
     CALL ERRMSG('Usage: wafs_grib2_cb0p25 grib2file1 grib2file2')
     CALL ERREXIT(2)
  ENDIF

  CALL GETARG(1,gfile1)
  CALL GETARG(2,gfile2)

  call cb_algo(trim(gfile1),trim(gfile2))

end program main

