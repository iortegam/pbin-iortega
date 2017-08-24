      subroutine fib(fname,iyrst,imst,idst,iyrnd,imnd,idnd,sst,ice)
      integer*4 iyrst,imst,idst,iyrnd,imnd,idnd 
      real*4 sst(360,180)
      integer*4 ice(360,180)
      integer*4 itag(360,180)
c      dimension sst(360,180)
c      dimension ice(360,180)
c      dimension itag(360,180)

      character fname*64

!f2py character intent(in) :: fname
!f2py real(4) intent(out) :: sst
!f2py integer intent(out) ::ice
!f2py integer intent(out) ::itag
!f2py integer intent(out) :: iyrst, imst, idst, iyrnd, imnd, idnd

c----------------------------------------------------------------------------
c  this program reads the version 2 monthly oi SST data
c
c  VARIABLE LIST
c  -------------
c  sst   - sea surface temperature array (deg C)
c  ice   - ice concentration array (%)  (0-100,   >100 = land or coast)
c  iyrst - year of start date of analysis
c  imst  - month of start date of analysis
c  idst  - day of start date of analysis
c  iyrnd - year of end date of analysis
c  imnd  - month of end date of analysis
c  idnd  - day of end date of analysis
c  ndays - number of days in analysis (start date thru enddate)
c  index - analysis version for reference
c
c
c  GRID ORIENTATION
c  ----------------
c  The first gridbox of each array is centered on 0.5E, 89.5S.  The points
c  move eastward to 359.5E, then northward to 89.5N.
c
c  So, the orientation of the sst and cice data arrays is :
c
c  sst(1,1)     =   0.5E, 89.5S
c  sst(2,1)     =   1.5E, 89.5S
c  sst(360,1)   = 359.5E, 89.5S
c  sst(1,2)     =   0.5E, 88.5S
c  sst(1,180)   =   0.5E, 89.5N
c  sst(360,180) = 359.5E, 89.5N
c
c
c----------------------------------------------------------------------------


c  uncomment out the lines below if you want to read the land-sea mask
c
c
c     Read in land sea tags (1 for ocean; 0 for land)
c      open (11,file=ftag,form='formatted')
c      READ (11,'(360i1)') itag
c      close(11)

c      write(*,'("Enter Input File Name : ",$)')
c      read(*,'(a)') fname
c      fname = '/home/ivan/Documents/pytest/oiv2mon.1981.asc'

      open(10,file=fname)


  10  read(10,6,end=99) iyrst,imst,idst,iyrnd,imnd,idnd,ndays,indx
  6   format(7I5,i10)
      write(*,*) iyrst,imst,idst,iyrnd,imnd,idnd,ndays,indx                                                             
      read(10,8) sst
      read(10,9) ice
  8   format(20f4.1)
  9   format(26i3)
      goto 10

  99  continue

      end subroutine fib