      subroutine fib(fname,iyrst,imst,idst,iyrnd,imnd,idnd,sst,ice)
      integer*4 iyrst(12),imst(12),idst(12),iyrnd(12),imnd(12),idnd(12) 
      real*4 sst(360,180,12)
      integer*4 ice(360,180,12)
      character fname*64


!f2py character intent(in) :: fname
!f2py real(4) intent(out) :: sst
!f2py integer intent(out) ::ice
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

      real*4 sst_i(360,180)
      integer*4 ice_i(360,180)

      open(10,file=fname)

      do K=1,12

      read(10,6) iyrst_i,imst_i,idst_i,iyrnd_i,imnd_i,idnd_i,ndays_i,indx_i
  6   format(7I5,i10)
c      write(*,*) iyrst_i,imst_i,idst_i,iyrnd_i,imnd_i,idnd_i,ndays_i,indx_i
      iyrst(K)  = iyrst_i
      imst(K)   = imst_i
      idst(K)   = idst_i
      iyrnd(K)  = iyrnd_i
      imnd(K)   = imnd_i
      idnd(K)   = idnd_i
      read(10,8) sst_i
      read(10,9) ice_i
  8   format(20f4.1)
  9   format(26i3)
c      write(*,*) sst_i
      sst(:,:,K) = sst_i
      ice(:,:,K) = ice_i
c      goto 10
      ENDDO

  

      end subroutine fib