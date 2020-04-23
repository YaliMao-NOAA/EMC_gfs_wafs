program  relabel
!$$$  MAIN PROGRAM DOCUMENTATION BLOCK
!
! MAIN PROGRAM:  copied from setmissing.f90
!   2020-03-17
!
! ABSTRACT: This program is designged for WAFS products at 0.25 degree.
!           Should not be used after UPP WAFS output on levels of exact numbers.
!           which needs to relabel ICAO standard reference pressure level to 
!           exact level. WGRIB2 has an option '-set_lev' but it will change
!           14 15 elements of TEMPLATE 4.0, 'Scale factor of second fixed surface'
!           and 'Scaled value of second fixed surfaces'respectively.
!           TOCGRIB2 uses all elements of TEMPLATE 4.0. And DEGRIB2 crashes.
!
!           The program reads GRIB2 file and simply relabel the pressure levels.
!
! PROGRAM HISTORY LOG:
! 2020-03-17  Y Mao
!
! USAGE:
!   INPUT FILES:
!     UNIT 10  - Input GRIB file
!
!   OUTPUT FILES:
!     UNIT 20  - Output GRIB file
!
! USAGE:
! COMMAND LINE:
!     relabel  inputfile outputfile
!   OUTPUT FILES:  (INCLUDING SCRATCH FILES)
!     6        - STANDARD FORTRAN PRINT FILE
!
!   SUBPROGRAMS CALLED: (LIST ALL CALLED FROM ANYWHERE IN CODES)
!     LIBRARY:
!       G2LIB    - GB_INFO, GT_GETFLD
!       W3LIB    - GBYTE, SKGB
!       BACIO    - BAOPENR, BAREAD, BACLOSE
!       SYSTEM   - IARGC   FUNCTION RETURNS NUMBER OF ARGUMENT ON
!                          COMMAND LINE
!                - GETARG  ROUTINE RETURNS COMMAND LINE ARGUMENT
!
!   EXIT STATES:
!     COND =   0 - SUCCESSFUL RUN
!
! REMARKS: COMMAND LINE CAN HAVE ONE FILE NAME.
!     
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE:  IBM 
!

  use grib_mod
  use params

  implicit none

  integer, parameter :: msk1=32000,msk2=4000
  CHARACTER(len=1),allocatable,dimension(:) :: cgrib
  integer :: listsec0(3),listsec1(13)
  integer :: currlen=0
  logical :: unpack,expand
  type(gribfield) :: gfld
  integer :: itot, icount, iseek,lskip,lgrib,lengrib,j
  integer :: numfields,numlocal,maxlocal

  character(100) :: inputfile,outputfile,avar
  INTEGER :: NARG,IARGC,n,i
  integer :: IFL1,IFL2,iret

  call start()
  unpack=.true.
  expand=.true.
      
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  GET ARGUMENTS
  NARG=IARGC()
  IF(NARG /= 2) THEN
     CALL ERRMSG('relabel:  Incorrect usage')
     CALL ERRMSG('Usage: relabel inputfile outputfile')
     CALL ERREXIT(2)
  ENDIF

  IFL1=10
  IFL2=20

  !
  CALL GETARG(1,inputfile)
  CALL BAOPENR(ifl1,trim(inputfile),iret)
  if(iret/=0)print*,'cant open ',trim(inputfile)
  !
  CALL GETARG(2,outputfile)
  call baopenw(IFL2, outputfile, iret)
  if(iret/=0)print*,'cant open ',trim(outputfile)

  itot=0
  icount=0
  iseek=0
  do
     call skgb(ifl1,iseek,msk1,lskip,lgrib)
     if (lgrib==0) exit    ! end loop at EOF or problem
     if (lgrib>currlen) then
        if (allocated(cgrib)) deallocate(cgrib)
        allocate(cgrib(lgrib),stat=iret)
        currlen=lgrib
     endif
     call baread(ifl1,lskip,lgrib,lengrib,cgrib)
     if (lgrib/=lengrib) then
        print *,' relabel: IO Error.'
        call errexit(9)
     endif
     iseek=lskip+lgrib
     icount=icount+1
     PRINT *
     PRINT *,'GRIB MESSAGE ',icount,' starts at',lskip+1
     PRINT *

! Unpack GRIB2 field
     call gb_info(cgrib,lengrib,listsec0,listsec1,&
                  numfields,numlocal,maxlocal,iret)
     if (iret/=0) then
        write(6,*) ' ERROR querying GRIB2 message = ',iret
        stop 10
     endif
     itot=itot+numfields

     do n=1,numfields
        call gf_getfld(cgrib,lengrib,n,unpack,expand,gfld,iret)
        if (iret/=0) then
           write(6,*) ' ERROR extracting field = ',iret
           cycle
        endif

     enddo

     if(iret==0)then
        if(gfld%ipdtmpl(12) == 10000) gfld%ipdtmpl(12)=10040
        if(gfld%ipdtmpl(12) == 12500) gfld%ipdtmpl(12)=12770
        if(gfld%ipdtmpl(12) == 15000) gfld%ipdtmpl(12)=14750
        if(gfld%ipdtmpl(12) == 17500) gfld%ipdtmpl(12)=17870
        if(gfld%ipdtmpl(12) == 20000) gfld%ipdtmpl(12)=19680
        if(gfld%ipdtmpl(12) == 22500) gfld%ipdtmpl(12)=22730
        if(gfld%ipdtmpl(12) == 27500) gfld%ipdtmpl(12)=27450
        if(gfld%ipdtmpl(12) == 30000) gfld%ipdtmpl(12)=30090
        if(gfld%ipdtmpl(12) == 35000) gfld%ipdtmpl(12)=34430
        if(gfld%ipdtmpl(12) == 40000) gfld%ipdtmpl(12)=39270
        if(gfld%ipdtmpl(12) == 45000) gfld%ipdtmpl(12)=44650
        if(gfld%ipdtmpl(12) == 50000) gfld%ipdtmpl(12)=50600
        if(gfld%ipdtmpl(12) == 60000) gfld%ipdtmpl(12)=59520
        if(gfld%ipdtmpl(12) == 70000) gfld%ipdtmpl(12)=69680
        if(gfld%ipdtmpl(12) == 75000) gfld%ipdtmpl(12)=75260
        if(gfld%ipdtmpl(12) == 80000) gfld%ipdtmpl(12)=81200
        if(gfld%ipdtmpl(12) == 85000) gfld%ipdtmpl(12)=84310
!        if(gfld%ipdtmpl(12) == 10040) gfld%ipdtmpl(12)=10000
!        if(gfld%ipdtmpl(12) == 12770) gfld%ipdtmpl(12)=12500
!        if(gfld%ipdtmpl(12) == 14750) gfld%ipdtmpl(12)=15000
!        if(gfld%ipdtmpl(12) == 17870) gfld%ipdtmpl(12)=17500
!        if(gfld%ipdtmpl(12) == 19680) gfld%ipdtmpl(12)=20000
!        if(gfld%ipdtmpl(12) == 22730) gfld%ipdtmpl(12)=22500
!        if(gfld%ipdtmpl(12) == 27450) gfld%ipdtmpl(12)=27500
!        if(gfld%ipdtmpl(12) == 30090) gfld%ipdtmpl(12)=30000
!        if(gfld%ipdtmpl(12) == 34430) gfld%ipdtmpl(12)=35000
!        if(gfld%ipdtmpl(12) == 39270) gfld%ipdtmpl(12)=40000
!        if(gfld%ipdtmpl(12) == 44650) gfld%ipdtmpl(12)=45000
!        if(gfld%ipdtmpl(12) == 50600) gfld%ipdtmpl(12)=50000
!        if(gfld%ipdtmpl(12) == 59520) gfld%ipdtmpl(12)=60000
!        if(gfld%ipdtmpl(12) == 69680) gfld%ipdtmpl(12)=70000
!        if(gfld%ipdtmpl(12) == 75260) gfld%ipdtmpl(12)=75000
!        if(gfld%ipdtmpl(12) == 81200) gfld%ipdtmpl(12)=80000
!        if(gfld%ipdtmpl(12) == 84310) gfld%ipdtmpl(12)=85000
     end if

     call putgb2(IFL2,gfld,iret)
     call gf_free(gfld)

  enddo

  call BACLOSE(IFL1,iret)
  call BACLOSE(IFL2,iret)

  print *," "
  print *, ' Total Number of Fields Found = ',itot
  print*,'sending the data for output'
  stop
end program relabel
