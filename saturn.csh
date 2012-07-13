#!/bin/csh
# saturn.csh
# Dick Plambeck, 31-may-2012

# ---------------------------------------------------------------------------------------#
# user-settable parameters
#
set RAW = raw
set REFANT = 1             # MUST be a 10-m telescope!
set SRC = Saturn              # source name
set PASSCAL = 3c279          # passband calibrator name
set GAINCAL = 3c279          # gain calibrator name
set GAINFLUX = 20         # flux density of gain calibrator, if you know it
set LEAKFILE = "None"      # visibility data file with leak corrections, if available
set MAP = Saturn              # continuum map name
set I_RMS = 0.9            # noise in total intensity map
set QU_RMS = 0.05           # measured noise in Stokes Q,U,V maps
set REGION = 'arcsec,box(20,-18,-20,18)'	# region to plot
set LSBLEAK = /o/chat/c_data/leakages/20120417-current.219.leak
set USBLEAK = /o/chat/c_data/leakages/20120417-current.229.leak
set DSBLEAK = /o/chat/c_data/leakages/20120417-current.224.leak

# ---------------------------------------------------------------------------------------#

goto start

    rm -r raw ct002.230Saturn.1.miriad
    tar xvf ct002.230Saturn.1.miriad.tar.gz 
    mv ct002.230Saturn.1.miriad raw
    uvflag vis=raw select='source(3c279o6),time(00:00,01:10)' flagval=flag
      # catalog position was messed up for 3c279o6 until 01:10 UT

  # ------------------------------------------------------------------------------------- #
  # step 1: 'xyphase' calibration
  #
  # The 'xyphase' is the LCP-RCP phase difference.  A direct calibration of the xyphase 
  # is possible only on the 10-m telescopes by observing a polarized noise source.  This
  # calibration is performed automatically every ~45 minutes by the standard observing
  # script.  Xyphase calibration data are labeled 'purpose(P)' in the Miriad dataset.
  # ------------------------------------------------------------------------------------- #

    xyauto vis=$RAW select='purpose(P)' 
        # ... fits purpose(P) data, stores phase correction as a bandpass
 
    smagpplt vis=$RAW options=bandpass,nofit,wrap device=/xs yrange=-180,180 \
      xaxis=chan yaxis=phase nxy=5,3
        # ... examine result; all phases will be zero except for C1-C6 LCP

    rm -r wide.xy
    uvcat vis=$RAW out=wide.xy options=nopol select='-source(noise),-auto,-purpose(P)' 
        # ... rewrite data to apply the correction

  # ------------------------------------------------------------------------------------- #
  # step 2: passband correction
  #
  # At this point phase(LCP) = phase(RCP) on the 10-m, but not the 6-m, telescopes.
  # It is essential to choose a 10-m telescope as the REFANT when fitting the passband 
  # in order to transfer the xyphase calibration to the 6-m telescopes.
  # ------------------------------------------------------------------------------------- #

    uvflag vis=wide.xy select='time(12jun23:00:37:18,12jun23:00:37:20)' flagval=flag
    uvflag vis=wide.xy select='time(12jun23:02:59:20,12jun23:02:59:30)' flagval=flag
        # ... flag 2 data points with crazy gains

    mfcal vis=wide.xy select='source('$PASSCAL')' interval=0.1 refant=$REFANT 
    smagpplt vis=wide.xy options=bandpass,nofit,wrap device=/xs yrange=-180,180 \
      xaxis=chan yaxis=phase
        # ... examine the passband solution

    rm -r wide.pb
    uvcat vis=wide.xy out=wide.pb options=nocal
        # ... rewrite to apply passband, but not gains 

  # ------------------------------------------------------------------------------------- #
  # step 3: gain correction
  #
  # Now fit instrumental amp and phase variations vs time
  # Use C8 as reference antenna because it was near the center of the array
  # ------------------------------------------------------------------------------------- #

    mfcal vis=wide.pb select='source('$PASSCAL')' interval=0.1 refant=8
    smagpplt vis=wide.pb options=bandpass,nofit,wrap device=/xs yrange=-180,180 \
      xaxis=chan yaxis=phase
        # ... confirm that passbands are flat, from previous step
    gpplt vis=wide.pb yaxis=phase options=wrap yrange=-180,180 device=/xs nxy=5,3
	    # ... examine the gain phases vs time
    gpplt vis=wide.pb yaxis=amp device=/xs nxy=5,3 #yrange=0,3
	    # ... examine the gain amplitudes vs time

    gpaver vis=wide.pb interval=10 options=scalar
        # ... time average the gains to smooth them
    gpplt vis=wide.pb yaxis=phase options=wrap yrange=-180,180 device=/xs nxy=5,3
	    # ... examine the gain phases vs time
    gpplt vis=wide.pb yaxis=amp device=/xs nxy=5,3 yrange=0,3
	    # ... examine the gain amplitudes vs time

    rm -r lsb.cal
    uvcat vis=wide.pb out=lsb.cal options=nopol,nopass,unflagged select='win(1,3,5,7)' 
    rm -r usb.cal
    uvcat vis=wide.pb out=usb.cal options=nopol,nopass,unflagged select='win(9,11,13,15)' 
        # ... break data into separate lsb and usb files; rewriting applies gains vs time

  # ------------------------------------------------------------------------------------- #
  # step 4: leakage calibration
  # 
  # Leakage corrections compensate for cross-coupling between the LCP and RCP channels.
  # Fit leakages to the calibrator if it was observed over wide parallactic angle range,
  # or copy leakages from standard tables.
  # ------------------------------------------------------------------------------------- #
  
    uvplt vis=wide.av select='source('$GAINCAL'),pol(LL)' axis=time,parang device=/xs \
	  options=nobase,nopol
        # ... plot parallactic angle coverage of gain calibrator

    if ($LEAKFILE == "None") then
      gpcal vis=lsb.cal options=circular,qusolve,noxy,nopass flux=$GAINFLUX interval=0.5 \
        refant=$REFANT select='source('$GAINCAL')' 
      gpcal vis=usb.cal options=circular,qusolve,noxy,nopass flux=$GAINFLUX interval=0.5 \
        refant=$REFANT select='source('$GAINCAL')' 
          # ... compute LSB and USB leakages from 3c279 data
    else
      gpcopy vis=$LSBLEAK out=lsb.cal options=nocal,nopass
      gpcopy vis=$USBLEAK out=usb.cal options=nocal,nopass
          # ... use standard leakage tables
    endif

  # ------------------------------------------------------------------------------------- #
  # step 5: make Saturn maps
  #
  # At this point the data are fully calibrated against 3c279. Make the maps.
  # ------------------------------------------------------------------------------------- #
  
    uvflag vis=lsb.cal select='time(12jun23:02:40,12jun23:02:41:20)' flagval=flag
    uvflag vis=usb.cal select='time(12jun23:02:40,12jun23:02:41:20)' flagval=flag
       # ... drop a few more data points with crazy gains

    uvplt vis=lsb.cal,usb.cal line=chan,1,1,96 axis=uvdist,amp options=nobase,nopol,nopass device=/xs \
	   select='source(saturn),dra(-1,1),pol(LL,RR)' nxy=1,1 
          # ... hope to see a bessel function here

    rm -r  $MAP.I.mp $MAP.Q.mp $MAP.U.mp $MAP.V.mp $MAP.bm
    invert vis=lsb.cal,usb.cal line=chan,4,1,24,24 \
      map=$MAP.I.mp,$MAP.Q.mp,$MAP.U.mp,$MAP.V.mp beam=$MAP.bm stokes=I,Q,U,V sup=0 \
      'select=source(SATURN),dra(-1,1),uvrange(0,100)' options=mfs,systemp,nopass,nocal cell=0.5 imsize=128
    cgdisp in=Saturn.I.mp device=/xs options=full labtyp=arcsec
    $<

    rm noiseList
    foreach MP ($MAP.I $MAP.Q $MAP.U $MAP.V)
      rm -r $MP.sl
      clean map=$MP.mp beam=$MAP.bm out=$MP.sl niters=5000
      rm -r $MP.cm
      restor map=$MP.mp beam=$MAP.bm model=$MP.sl out=$MP.cm 
      cgdisp in=$MP.cm device=/xs labtyp=arcsec 
      echo " " >> noiseList
      echo $MP".cm" >> noiseList
      imlist options=stat in=$MP.cm region='arcsec,box(30,-30,-30,-15)' | tail -2 >> noiseList
        # ... measure actual noise in a box offset from the center; change region if
        # ... source extends into it
      $<
    end

    # --- normally we would stop here, but Saturn allows for iterative phase-only selfcal --- #
    imhist in=Saturn.I.sl
      # ... choose a clipping level from the histogram
    selfcal vis=lsb.cal select='source(Saturn),dra(-1,1),pol(LL,RR)' model=Saturn.I.sl clip=1. interval=0.1 refant=8  
    gpplt vis=lsb.cal device=/xs yaxis=phase yrange=-180,180 options=wrap
    selfcal vis=usb.cal select='source(Saturn),dra(-1,1),pol(LL,RR)' model=Saturn.I.sl clip=1. interval=0.1 refant=8  
    gpplt vis=usb.cal device=/xs yaxis=phase yrange=-180,180 options=wrap

    # --- remake maps after 1 selfcal --- #
    rm -r  $MAP.I.mp $MAP.Q.mp $MAP.U.mp $MAP.V.mp $MAP.bm
    invert vis=lsb.cal,usb.cal line=chan,4,1,24,24 \
      map=$MAP.I.mp,$MAP.Q.mp,$MAP.U.mp,$MAP.V.mp beam=$MAP.bm stokes=I,Q,U,V sup=0 \
      'select=source(SATURN),dra(-1,1),uvrange(0,100)' options=mfs,systemp,nopass cell=0.5 imsize=128
    cgdisp in=Saturn.I.mp device=/xs options=full labtyp=arcsec
    $<

    rm noiseList
    foreach MP ($MAP.I $MAP.Q $MAP.U $MAP.V)
      rm -r $MP.sl
      clean map=$MP.mp beam=$MAP.bm out=$MP.sl niters=5000
      rm -r $MP.cm
      restor map=$MP.mp beam=$MAP.bm model=$MP.sl out=$MP.cm 
      cgdisp in=$MP.cm device=/xs labtyp=arcsec 
      echo " " >> noiseList
      echo $MP".cm" >> noiseList
      imlist options=stat in=$MP.cm region='arcsec,box(30,-30,-30,-15)' | tail -2 >> noiseList
        # ... measure actual noise in a box offset from the center; change region if
        # ... source extends into it
      $<
    end
    more noiseList
goto end
start:
    # --- plot Q and U on grayscale I image --- #
    cgdisp in=$MAP.I.cm,$MAP.Q.cm,$MAP.U.cm type=pixel,contour,contour options=full \
      region=$REGION labtyp=hms,dms cols1=2 cols2=4 slev=a,$QU_RMS,a,$QU_RMS \
      levs1=-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,3,4,5,6,7,8,9,10,11,12,13,14,15 \
      levs2=-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,3,4,5,6,7,8,9,10,11,12,13,14,15 \
      line=1,3,3 device=/xs
        # ... >3 sigma detection of Stokes Q and/or U needed for a polarization detection
        # ... remember that Q and U can be positive or negative
	$<

    rm -r $MAP.poli.cm $MAP.polm.cm $MAP.pa.cm
    impol poli=$MAP.poli.cm polm=$MAP.polm.cm pa=$MAP.pa.cm sigma=$QU_RMS,$I_RMS \
      in=$MAP.Q.cm,$MAP.U.cm,$MAP.I.cm sncut=3,2
        # ... derive pol intensity and PA maps from Stokes parameters     

    cgdisp in=$MAP.I.cm,$MAP.polm.cm,$MAP.pa.cm type=pixel,amp,angle \
      region=$REGION options=full labtyp=hms,dms vecfac=1.2,4,4 beamtyp=b,l,4 \
      lines=1,8 cols1=1 slev=a,$I_RMS \
      levs1=-6,-5,-4,-3,3,4,5,6,8,10,15,20,25,30,35,40,45,50,55 \
      device=/xs range=0,50,lin,3 
        # ... plot polarization vectors on total intensity map

goto end

  # ------------------------------------------------------------------------------------- #
  # step 6: measure primary beam polarization from offset 3c279 maps, and PA on 3c286
  #
  # getPA.csh source offset interval
  # ------------------------------------------------------------------------------------- #

    getPA.csh  3C279O1 -36.26,-0.02  0.1
    getPA.csh  3C279O2 -18.36,-0.02  0.1
    getPA.csh  3C279O3  -0.45,-0.02  0.1
    getPA.csh  3C279O4  -0.45,8.98   0.1
    getPA.csh  3C279O5  17.46,-0.02  0.1
    getPA.csh  3C279O6  35.36,-0.02  0.1
    getPA.csh  3c286    0.,0.        10. 
goto end


  # ------------------------------------------------------------------------------------- #
  # cleanup
  # ------------------------------------------------------------------------------------- #

    rm -r wide.xy
    rm -r wide.pb
    rm -r $MAP.*.sl $MAP.*.mp

end:
