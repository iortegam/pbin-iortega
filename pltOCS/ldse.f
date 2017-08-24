      subroutine ldse(ftag,itag)
      character ftag*64
      integer*4 itag(360,180)

!f2py character intent(in) ::ftag
!f2py integer intent(out) ::itag

c---------------------------------------------------------------------------
c  This program reads the ascii land-sea mask
c
c
c
c The OI analysis is done over all ocean areas and the Great Lakes.
c There is no analysis over land.  The land values are filled by a Cressman
c interpolation to produce a complete grid for possible interpolation to
c other grids.  The ocean and land areas are defined by a land sea mask.
c This data set is a binary, direct access file, lstags.ondeg.dat, which
c is included in the same directory.  The spatial grid is defined identically
c to the grid for the SST arrays.  The values in ls.dat are set to 1 over
c the ocean and 0 over land. It can be read by the following fortran code:
c
c
c
c  GRID ORIENTATION
c  ----------------
c  The first gridbox of each array is centered on 0.5E, 89.5S.  The points
c  move eastward to 359.5E, then northward to 89.5N.
c
c  So, the orientation of the sst and cice data arrays is :
c
c  itagls(1,1)     =   0.5E, 89.5S
c  itagls(2,1)     =   1.5E, 89.5S
c  itagls(360,1)   = 359.5E, 89.5S
c  itagls(1,2)     =   0.5E, 88.5S
c  itagls(1,180)   =   0.5E, 89.5N
c  itagls(360,180) = 359.5E, 89.5N
c
c
c----------------------------------------------------------------------------



c      dimension itag(360,180)
      open (10,file=ftag,form='formatted')
c
c     Read in land sea tags (1 for ocean; 0 for land)
c
      READ (10,'(360i1)') itag

      END