	program writeleak
c	
c	
c       Usage:  writeleak vis=    leak=    
c	Example: writeleak vis=3c273 leak=leaklsb
c		where vis is standard MIRIAD visibility file and
c		      leak is a textfile with the following format.
c         Format of leak:
c	      48 numbers on one line separated by commas; with 4 numbers 
c	      as the leakage for each antenna
c	
c	Modified:
c       01 Jan 06 dpm - take variable number of arguments
c	
	integer MAXANT
	parameter (MAXANT=23)

	integer i,j,tIn,iostat,nleaks,nant
	character vis*64,buf*256
	real leaks(92)
	complex D(2,maxant)

	call keyini 
	call keya('vis',vis,' ')
	if(vis.eq.' ')call bug('f','Input data-set must be given')
	call mkeyr('leak',leaks,MAXANT*4,nleaks)
	if(nleaks.eq.0)call bug('f','Input leak must be given')
        call hopen(tIn,vis,'old',iostat)
	call keyfin

	if(mod(nleaks,4).ne.0) call bug('f',
	1    '4 leakages required per antenna')
	nant=nleaks/4
c	write (buf,'(a,i3)')'Antennas found: ',nant
c	call output(buf)
c	call output('                   Dx                Dy')
c	do i=1,nant
c	   write (buf,'(a,i3,a,4f9.5)') '  Ant',i,': ',
c	1	(leaks(j),j=(i-1)*4+1,i*4)
c	   call output(buf)
c	end do
	do i=1,nant
	   D(1,i)=leaks(4*(i-1)+1)*(1.,0.)+leaks(4*(i-1)+2)*(0.,1.)
	   D(2,i)=leaks(4*(i-1)+3)*(1.,0.)+leaks(4*(i-1)+4)*(0.,1.)
	end do
	call LeakTab(tIn,D,nant)

	call hclose(tIn)
	stop
	end

c************************************************************************
        subroutine LeakTab(tIn,D,nants)
c
        implicit none
        integer tIn,nants
        complex D(2,nants)
c
c  Save the polarisation gains in an item in the calibrator file.
c  This uses some dirty tricks. Firstly it uses wrhdc to create the
c  item (because there is no way in FORTRAN to make the correct header).
c  Second it uses hwriter, because there is no hwritec.
c
c  Input:
c    tIn        Handle of the calibrator file.
c    D          Leakage parameters.
c    nants      Number of antennae.
c------------------------------------------------------------------------
        integer item,iostat
c
        call wrhdc(tIn,'leakage',(1.,0.))
        call haccess(tIn,item,'leakage','append',iostat)
        if(iostat.ne.0)then
          call bug('w','Error opening the output leakage table')
          call bugno('f',iostat)
        endif
c
c  The item alignment requirement requires that we start writing at
c  byte 8. Bytes 0-3 are the header, and bytes 4-7 are unused.
c
        call hwriter(item,D,8,2*8*nants,iostat)
        if(iostat.ne.0)then
          call bug('w','Error writing to the leakage table')
          call bugno('f',iostat)
        endif
        call hdaccess(item,iostat)
        if(iostat.ne.0)call bugno('f',iostat)
        end
c************************************************************************
