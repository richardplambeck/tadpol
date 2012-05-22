#!/bin/csh
# reduce.csh  - 16 may 2012 version

# designed to reduce TADPOL data with 3 wideband correlator sections,
#  and 1 narrowband section with CO in USB, SiO in LSB

# inputs: $RAW1, $RAW2, etc -> concatenate into $RAW
# fit xyphase to $RAW -> rewrite to 'wide.xy', 'win7.xy', 'win15.xy'
# fit passband to 'wide.xy' -> rewrite to 'wide.pb'
# segregate MWC349 data, if it exists -> 'MWC349.pb' 
# ... for possible later analysis of H30alpha recomb line polarization
# average wide.pb into 6 channels -> 'wide.av'
# phaseonly selfcal on wide.av, then use bootflux to get calibrator fluxes
# edit FluxSource.cat, use amplitude selfcal to identify really bad gains,
#   flag data as necessary
# break wide.av into wide.lsb and wide.usb (leakages differ); fit leakages
#   on gain calibrator if possible, send to Chat to include in master list

# =========== begin user-supplied parameters ============ #

setenv MIRFLUXTAB /fringe2/plambeck/pol/FluxSource.cat
  # ... my own private version of FluxSource.cat

set RAW1 = c0888I.22D_230DR21OH.1.miriad
# set RAW2 = ...
set RAW = raw

set REFANT = 2
  # ... must be one of the 10-m telescopes!
set SRC = DR21
  # ... actual source name in the data file
set ALLCALS = ( BLLAC )
  # ... list of sources to process with bootflux
set PASSCAL = BLLAC	    
set FLUXCAL = MWC349
set LEAKCAL = BLLAC
set LEAKFLUX = 7.1 
set INOISE = 0.01
set QUNOISE = 0.00060
set REGION = 'arcsec,box(12,-12,-12,12)'
set MAP = DR21OH
set BADANTS = 23
  # ... antenna numbers not to use in flux calibration
  # ... C23 is just a placeholder, meaning use all antennas

set LSBLEAK = /o/chat/c_data/leakages/20120417-20120514.219.leak
set USBLEAK = /o/chat/c_data/leakages/20120417-20120514.229.leak
set DSBLEAK = /o/chat/c_data/leakages/20120417-20120514.224.leak

set CO_I_SIGMA = 0.2
set CO_QU_SIGMA = 0.06
set SIO_I_SIGMA = 0.04
set SIO_QU_SIGMA = 0.017

# =========== end user-supplied parameters ============ #

goto start

  # Untar and (if necessary) concatenate data sets
    foreach FILE ($RAW1)
      rm -r $FILE
      tar xvf $FILE.tar.gz
    end
    rm -r $RAW
    mv $RAW1 $RAW
    # uvcat vis=$RAW1,$RAW2 out=$RAW options=nocal,nopass,nopol

  # Initial data inspection
    uvindex vis=$RAW log=uvindex.log
    view uvindex.log
goto end

# >>> edit source and calibrator names <<< #

  # Plot tsys; must separate into LL and RR or risk overflowing buffers 
    rm -r win1LL
    uvcat vis=$RAW out=win1LL select='win(1),pol(LL),-source(noise),-purpose(P)' \
      options=nocal,nopol,nopass
    rm -r win1RR
    uvcat vis=$RAW out=win1RR select='win(1),pol(RR),-source(noise),-purpose(P)' \
      options=nocal,nopol,nopass
    varplt vis=win1LL yaxis=systemp device=/xs yrange=0,1000 nxy=5,3
    $<
    varplt vis=win1RR yaxis=systemp device=/xs yrange=0,1000 nxy=5,3
    $<
    varplt vis=win1LL yaxis=rmspath device=/xs
    $<
    varplt vis=win1LL yaxis=tau230 device=/xs
goto end

# >>> flag data as appropriate <<< #

  # =================================
  # ====== XYphase calibration ====== 
  # =================================

  # Fit grid data, denoted by purpose(P), examine passband fit
    xyauto vis=$RAW select='purpose(P)'
    smagpplt vis=$RAW options=bandpass,nofit,wrap device=/xs yrange=-180,180 \
      xaxis=chan yaxis=phase
	  # ... all phases will be 0 except LCP for ants 1-6

  # Plot grid data, confirm that corrected phases are near 0
    uvspec vis=$RAW select='auto,-ant(7,8,9,10,11,12,13,14,15),purpose(P)' \
      axis=chan,phase device=/xs nxy=3,2 yrange=-50,50 interval=10

  # Rewrite to apply xyphase calibration; segregate wides and narrows here
	rm -r wide.xy
    uvcat vis=$RAW out=wide.xy options=nopol \
	  select='-source(noise),win(1,3,5,9,11,13),-auto,-purpose(P)' 
    rm -r win7.xy
    uvcat vis=$RAW out=win7.xy options=nopol \
	  select='-source(noise),win(7),-auto,-purpose(P)' 
    rm -r win15.xy
    uvcat vis=$RAW out=win15.xy options=nopol \
	  select='-source(noise),win(15),-auto,-purpose(P)' 

      # ... note: at this point phase(LL) = phase(RR) on 10-m, but not 6-m, telescopes
      # ... we have not yet solved for a passband - there will be band-to-band offsets

  # ===========================================
  # ====== Wideband passband calibration ====== 
  # ===========================================

    mfcal vis=wide.xy select='source('$PASSCAL'),-purpose(P)' interval=0.1 refant=$REFANT 
    smagpplt vis=wide.xy options=bandpass,nofit,wrap device=/xs yrange=-180,180 \
      xaxis=chan yaxis=phase
	$<
    smagpplt vis=wide.xy options=bandpass device=/xs yrange=0.5,1.5 \
      xaxis=chan yaxis=amp nxy=5,3
	$<
	gpplt vis=wide.xy device=/xs yaxis=phase yrange=-180,180 options=wrap nxy=5,3
	$<
	gpplt vis=wide.xy device=/xs yaxis=amp yrange=0,5 nxy=5,3
goto end

# >>> flag any data (wide and narrow) with crazy gains, refit passbands if necessary <<< #

  # Rewrite wide data to apply passband 
    rm -r wide.pb
    uvcat vis=wide.xy out=wide.pb options=nocal,nopol 

  # Check the passband in 2 hr blocks on passband and gain calibrators
    selfcal vis=wide.pb interval=0.1 refant=$REFANT select='-source('$SRC'),pol(LL,RR)'
    uvspec vis=wide.pb axis=chan,phase select='-source('$SRC,$FLUXCAL'),pol(LL,RR)' \
      interval=120 device=/xs yrange=-180,180 nxy=4,4

  # Segregate MWC349 channel data, if present, for later analysis of H30alpha polarization
	if ($FLUXCAL == MWC349) then
	  rm -r MWC349.pb
	  uvcat vis=wide.pb select='source(MWC349)' out=MWC349.pb
	  selfcal vis=MWC349.pb interval=0.1 select='pol(LL,RR)' options=amp,noscale flux=2.0 \
		refant=$REFANT line=chan,1,1,235
	  uvspec vis=MWC349.pb axis=chan,amp interval=100 device=/xs options=avall,nobase \
	    hann=3 nxy=1,1 yrange=-2,40
	endif

  # Condense wide to just 6 channels to speed further analysis
    rm -r wide.av
    uvaver vis=wide.pb out=wide.av line=chan,6,1,47,47 options=nocal


  # ==============================
  # ====== Flux calibration ====== 
  # ==============================

  # Inspect amplitude gains on flux calibrator (planet or MWC349) --- #	
	if ($FLUXCAL != MWC349) then
		selfcal vis=wide.av select='source('$FLUXCAL'),pol(LL,RR)' interval=0.1 \
		  refant=$REFANT options=amplitude,apriori,noscale line=chan,1,1,6
    	uvplt vis=wide.av axis=uvdist,amp device=/xs options=nobase nxy=1,1 \
		  select='source('$FLUXCAL')'
    else
	    selfcal vis=wide.av interval=0.1 select='source(MWC349),pol(LL,RR)' \
		  options=amp,noscale flux=2.0 refant=$REFANT line=chan,1,1,5
		    # ... omit window(6) = original win(13) because it covers the H30alpha recomb line
		gpplt vis=wide.av yaxis=amp device=/xs yrange=0,3 nxy=5,3
    endif
	gplist vis=wide.av
goto end

# >>> if desired, edit BADANTS to exclude any antennas with crazy gains on fluxcal <<< #

  # Phaseonly selfcal on all calibrators, then use bootflux to get fluxes 
    selfcal vis=wide.av select='-source('$SRC'),pol(LL,RR)' interval=0.1 refant=$REFANT \
      line=chan,1,1,6
    gpplt vis=wide.av yaxis=phase options=wrap yrange=-180,180 device=/xs nxy=5,3
    $<
    foreach CAL ( $ALLCALS ) 
	  if ($FLUXCAL == "MWC349" ) then
        bootflux vis=wide.av line=chan,1,1,5 primary=MWC349,1.9 \
	 	  select='source(MWC349,'$CAL'),pol(LL,RR),-ant('$BADANTS')' taver=5
	  else
        bootflux vis=wide.av line=chan,1,1,5 primary=$FLUXCAL \
	 	  select='source('$FLUXCAL','$CAL'),pol(LL,RR),-ant('$BADANTS')' taver=5
	  endif
      $<
    end
goto end

# ... MWC349 = 2.0 Jy -> 1.15 Jy/K 
# ... then S(BLLAC) = 7.1 Jy
	
# >>> edit /fringe2/plambeck/pol/FluxSource.cat with new fluxes <<< #
# >>> edit LEAKFLUX <<< #

  # Final check for any crazy antenna gains, using source fluxes in FluxSource.cat
	selfcal vis=wide.av refant=$REFANT interval=0.1 options=amplitude,apriori,noscale \
	  select='-source('$SRC',MWC349),pol(LL,RR)' 
    gpplt vis=wide.av yaxis=amp device=/xs yrange=0,3 nxy=5,3 nxy=5,3
    $<
    gpaver vis=wide.av options=scalar interval=15
    gpplt vis=wide.av yaxis=amp device=/xs yrange=0,3 nxy=5,3
goto end

# >>> flag any data with completely crazy gains; remember to flag win7 and win15 also <<< #

  # Segregate LSB and USB wideband data
    rm -r wide.lsb
    uvcat vis=wide.av out=wide.lsb select='win(1,2,3)' options=nopass,nocal,nopol
    rm -r wide.usb
    uvcat vis=wide.av out=wide.usb select='win(4,5,6)' options=nopass,nocal,nopol

  # =================================
  # ====== Leakage calibration ====== 
  # =================================

  # Check parallactic angle coverage of leakage calibrator 
    uvplt vis=wide.av select='source('$LEAKCAL'),pol(LL)' axis=time,parang device=/xs \
	  options=nobase,nopol
    $<

  # fit DSB leakages 
	gpcal vis=wide.av options=circular,qusolve,noxy,nopass flux=$LEAKFLUX interval=0.5 \
	  refant=$REFANT select='source('$LEAKCAL')' line=chan,1,1,6 > gpcal.DSB.log
    gpplt vis=wide.av options=polarization device=/xs nxy=2,2 yrange=0,0.2 
    uvflux vis=wide.av select='source('$LEAKCAL')' line=chan,1,1,6 stokes=I,Q,U,V \
	  options=uvpol,nopass
 
    # .., uvflux is informational only, to compare calibrator polarization derived
    # ... from DSB, LSB, USB (it's comforting if they are similar)

  # fit LSB leakages
	gpcal vis=wide.lsb options=circular,qusolve,noxy,nopass flux=$LEAKFLUX interval=0.5 \
	  refant=$REFANT select='source('$LEAKCAL')' line=chan,1,1,3 > gpcal.LSB.log
    gpplt vis=wide.lsb options=polarization device=/xs nxy=2,2 yrange=0,0.2
    uvflux vis=wide.lsb select='source('$LEAKCAL')' line=chan,1,1,3 stokes=I,Q,U,V \
	  options=uvpol,nopass

  # fit USB leakages
	gpcal vis=wide.usb options=circular,qusolve,noxy,nopass flux=$LEAKFLUX interval=0.5 \
	  refant=$REFANT select='source('$LEAKCAL')' line=chan,1,1,3 > gpcal.USB.log
    gpplt vis=wide.usb options=polarization device=/xs nxy=2,2 yrange=0,0.2
    uvflux vis=wide.usb select='source('$LEAKCAL')' line=chan,1,1,3 stokes=I,Q,U,V \
	  options=uvpol,nopass

# >>> Chat maintains master leakage list; send him uvindex.log and gpcal.log files <<< #

    # ... we think it is better to use the grand average leakages derived from many 
    # ... datasets, rather than fit a leakage to each individual data set

  # Copy leakages from master leak file
    gpcopy vis=$DSBLEAK out=wide.av options=nocal,nopass
    uvflux vis=wide.av select='source('$LEAKCAL')' line=chan,1,1,6 stokes=I,Q,U,V \
	  options=uvpol,nopass
    
    gpcopy vis=$LSBLEAK out=wide.lsb options=nocal,nopass
    uvflux vis=wide.lsb select='source('$LEAKCAL')' line=chan,1,1,3 stokes=I,Q,U,V \
	  options=uvpol,nopass
 
    gpcopy vis=$USBLEAK out=wide.usb options=nocal,nopass
    uvflux vis=wide.usb select='source('$LEAKCAL')' line=chan,1,1,3 stokes=I,Q,U,V \
	  options=uvpol,nopass

  # Once again fit gain vs time, smooth to 15 min time resolution, copy to lsb and usb
	selfcal vis=wide.av refant=$REFANT interval=0.1 options=amplitude,apriori,noscale \
	  select='-source('$SRC',MWC349),pol(LL,RR)' 
    gpaver vis=wide.av options=scalar interval=15
    gpplt vis=wide.av yaxis=amp device=/xs yrange=0,3 nxy=5,3
    gpcopy vis=wide.av out=wide.lsb options=nopass,nopol
    gpcopy vis=wide.av out=wide.usb options=nopass,nopol

  # this is an old Jedi mind trick, known by Tony Wong
    puthd in=wide.lsb/senmodel value='GSV' type=ascii
    puthd in=wide.usb/senmodel value='GSV' type=ascii

    # ... at this point wide.lsb and wide.usb are fully calibrated; xyphase and
    # ... passband corrections were applied in deriving the window-averaged data;
    # ... each file contains a gains item and a leakage item; senmodel=GSV supposedly
    # ... will induce Miriad to include the gains in the variance calculation done
    # ... for 'invert options=systemp' - this should downweight data with lousy pointing

  # ================================
  # ====== Dust continuum map ====== 
  # ================================

	rm -r  $MAP.I.mp $MAP.Q.mp $MAP.U.mp $MAP.V.mp $MAP.bm
	invert vis=wide.lsb,wide.usb line=chan,3,1,1 \
  	  map=$MAP.I.mp,$MAP.Q.mp,$MAP.U.mp,$MAP.V.mp beam=$MAP.bm stokes=I,Q,U,V sup=0 \
	  'select=source('$SRC')' options=mfs,systemp cell=0.25 imsize=512
	implot in=$MAP.I.mp device=/xs

	rm noiseList
	foreach MP ($MAP.I $MAP.Q $MAP.U $MAP.V)
	  rm -r $MP.sl
	  clean map=$MP.mp beam=$MAP.bm out=$MP.sl niters=3000
	  rm -r $MP.cm
	  restor map=$MP.mp beam=$MAP.bm model=$MP.sl out=$MP.cm 
	  cgdisp in=$MP.cm device=/xs region='arcsec,box(-15,-15,15,15)' labtyp=arcsec 
	  echo " " >> noiseList
	  echo $MP".cm" >> noiseList
	  imlist options=stat in=$MP.cm region='arcsec,box(20,-20,-20,-5)' | tail -2 >> noiseList
	  $<
	end
	tail -20 noiseList
goto end

# >>> edit INOISE and QUNOISE <<< #

  # Plot contour images of Q and U on gray scale image of I 
	cgdisp in=$MAP.I.cm,$MAP.Q.cm,$MAP.U.cm type=pixel,contour,contour options=full \
	  region=$REGION labtyp=hms,dms cols1=2 cols2=4 slev=a,$QUNOISE,a,$QUNOISE \
      levs1=-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,3,4,5,6,7,8,9,10,11,12,13,14,15 \
      levs2=-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,3,4,5,6,7,8,9,10,11,12,13,14,15 \
      device=/xs
	cgdisp in=$MAP.I.cm,$MAP.Q.cm,$MAP.U.cm type=pixel,contour,contour options=full \
	  region=$REGION labtyp=hms,dms cols1=2 cols2=4 slev=a,$QUNOISE,a,$QUNOISE \
      levs1=-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,3,4,5,6,7,8,9,10,11,12,13,14,15 \
      levs2=-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,3,4,5,6,7,8,9,10,11,12,13,14,15 \
      device=$MAP.IQU.ps/cps
	$<

	rm -r $MAP.poli.cm $MAP.polm.cm $MAP.pa.cm
	impol poli=$MAP.poli.cm polm=$MAP.polm.cm pa=$MAP.pa.cm sigma=$QUNOISE,$INOISE \
	  in=$MAP.Q.cm,$MAP.U.cm,$MAP.I.cm sncut=3,2

  # Polarization vectors on I contour map
	cgdisp in=$MAP.I.cm,$MAP.poli.cm,$MAP.pa.cm type=contour,amp,angle \
	  region=$REGION options=full,rot90 labtyp=hms,dms vecfac=1.2,4,4 beamtyp=b,l,4 \
      lines=1,1,10 cols1=1 slev=a,$INOISE \
      levs1=-6,-5,-4,-3,3,4,5,6,8,10,15,20,25,30,35,40,45,50,55 \
      device=/xs
	cgdisp in=$MAP.I.cm,$MAP.poli.cm,$MAP.pa.cm type=contour,amp,angle \
	  region=$REGION options=full,rot90 labtyp=hms,dms vecfac=1.2,4,4 beamtyp=b,l,4 \
      lines=1,1,10 cols1=1 slev=a,$INOISE \
      levs1=-6,-5,-4,-3,3,4,5,6,8,10,15,20,25,30,35,40,45,50,55 \
      device=$MAP.ps/cps

goto end

  # Cleanup
    rm -r raw
    rm -r wide.xy
    rm -r wide.pb
    rm -r $MAP.*.sl $MAP.*.mp

  # ============================================
  # ====== channel maps, preliminary steps ====== 
  # ============================================

  # once again do wideband phaseonly selfcal on calibrators
    selfcal vis=wide.av select='-source('$SRC'),pol(LL,RR)' interval=0.1 refant=$REFANT \
      line=chan,1,1,6

  # for later use, create wideband calibrator file with phases (but not amps) corrected
    rm -r wide.tmp
    uvcat vis=wide.av out=wide.tmp select='-source('$SRC','$FLUXCAL')'
	selfcal vis=wide.tmp refant=$REFANT interval=0.1 options=amplitude,apriori,noscale \
	  select='-source('$SRC','$FLUXCAL'),pol(LL,RR)' 
    gpplt vis=wide.tmp yaxis=phase options=wrap yrange=-180,180 device=/xs nxy=5,3
    gpplt vis=wide.tmp yaxis=amp device=/xs yrange=0,3 nxy=5,3
    gpaver vis=wide.tmp options=scalar interval=15
    gpplt vis=wide.tmp yaxis=amp device=/xs yrange=0,3 nxy=5,3

  # now time average the gains in wide.av
    gpaver vis=wide.av interval=15
 
    # ... at this point, gains in wide.av correct phase vs time (amps=1)
    # ... gains in wide.tmp correct amp vs time (phases=0)


  # =============================
  # ====== CO channel maps ====== 
  # =============================

  # Copy phase vs time from wide.av, rewrite data to apply
	gpcopy vis=wide.av out=win15.xy options=nopass,nopol
	rm -r win15.tmp
    uvcat vis=win15.xy out=win15.tmp options=nopass,nopol
    uvplt vis=win15.tmp axis=time,phase device=/xs line=chan,1,1,191 average=10 \
	  select='-source('$SRC')' options=nocal,nopass,nopol nxy=2,2
      # ... phases for LL,RR should now be constant vs time; RL and LR will change
      # ... with time if source is polarized; phs(LL)=phs(RR) on baselines 
      # ... between 10-m telescopes because of previous xyphase cal

# -----------------------------------------------------------------------------
# special for BLLAC - flag chans with foreground CO absorption
#  uvflag vis=win15.tmp select='source(BLLAC)' line=chan,5,82,1,1 flagval=flag
# -----------------------------------...........-------------------------------

  # Fit passband (and avg gain) to grand avg of win15 data
    mfcal vis=win15.tmp line=chan,10,1,19,19 select='source('$PASSCAL')' \
      interval=10000 refant=$REFANT
    smagpplt vis=win15.tmp options=bandpass,nofit,wrap device=/xs yrange=-180,180 \
      xaxis=chan yaxis=phase nxy=5,3
    smagpplt vis=win15.tmp options=bandpass,nofit,wrap device=/xs yrange=0,3 \
      xaxis=chan yaxis=amp nxy=5,3
	$<
	gpplt vis=win15.tmp yaxis=amp device=/xs yrange=0,3 nxy=5,3
	$<
	gpplt vis=win15.tmp yaxis=phase device=/xs yrange=-180,180 options=wrap nxy=5,3
	$<
goto end

  # Rewrite the data to apply the passband 
    rm -r win15.cal
    uvcat vis=win15.tmp out=win15.cal options=nopol

  # Copy amplitude gains vs time from wide.tmp, leakages from wide.USB
    gpcopy vis=wide.tmp out=win15.cal options=nopass,nopol
	gpcopy vis=wide.usb out=win15.cal options=nopass,nocal
    puthd in=win15.cal/senmodel value='GSV' type=ascii
    uvplt vis=win15.cal axis=time,phase device=/xs line=chan,1,1,191 average=10 \
	  select='-source('$SRC')'
    $<
    uvspec vis=win15.cal interval=10000 device=/xs axis=chan,phase yrange=-180,180 \
	  select='source('$PASSCAL')' nxy=4,4 line=chan,38,1,5,5

  # Make CO channel maps 
	puthd in=win15.cal/restfreq value=230.538 type=double
	rm -r co.I.mp co.Q.mp co.U.mp co.V.mp co.bm
	invert vis=win15.cal line=chan,38,1,5,5 \
	  map=co.I.mp,co.Q.mp,co.U.mp,co.V.mp beam=co.bm stokes=I,Q,U,V sup=0 \
	  'select=source('$SRC')' options=systemp cell=0.25 imsize=512
	cgdisp in=co.I.mp device=/xs region=arcsec,'box(-30,-30,30,30)' options=3pixel \
	  csize=0,1,1 nxy=3,3

	rm noiseList
	foreach MP (co.I co.Q co.U co.V)
	  rm -r $MP.sl
	  clean map=$MP.mp beam=$MAP.bm out=$MP.sl niters=1000
	  rm -r $MP.cm
	  restor map=$MP.mp beam=$MAP.bm model=$MP.sl out=$MP.cm
	  cgdisp in=$MP.cm device=/xs region=quarter labtyp=arcsec options=3pixel csize=0,1,1
	  echo " " >> noiseList
	  imlist options=stat in=$MP.cm region='arcsec,box(20,-20,-20,20)' log=junk
      cat junk >> noiseList
	  $<
	end
    view noiseList
goto end

# >>> edit CO_I_SIGMA, CO_QU_SIGMA <<< #

    cgdisp in=co.I.cm,co.Q.cm,co.U.cm type=pixel,contour,contour \
      device=/xs region='arcsec,box(-20,-20,20,20)'  \
	  options=3pixel csize=0,1,1 range=$CO_I_SIGMA,8,lin,1 \
      cols1=2 cols2=4 slev=a,$CO_QU_SIGMA,a,$CO_QU_SIGMA \
      levs1=-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,3,4,5,6,7,8,9,10,11,12,13,14,15 \
      levs2=-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,3,4,5,6,7,8,9,10,11,12,13,14,15 \
      nxy=2,2

  # Create red and blue outflow maps
	rm -r cored.cm
	moment in=co.I.cm out=cored.cm mom=0 clip=-1000,.5 region='image(11,18)'
    cgdisp in=cored.cm device=/xs region=quarter range=6,100,lin,1

	rm -r coblue.cm
	moment in=co.I.cm out=coblue.cm mom=0 clip=-1000,.5 region='image(22,36)'
    cgdisp in=coblue.cm device=/xs region=quarter range=6,100,lin,1

	cgdisp device=/xs in=$MAP.I.cm,$MAP.poli.cm,$MAP.pa.cm,cored.cm,coblue.cm \
      type=pixel,amp,angle,contour,contour \
	  region='arcsec,box(18,-14,-18,14)' options=full,rot90 labtyp=hms,dms \
	  vecfac=1.2,4,4 beamtyp=b,l,4 lines=1,1,1,10 cols1=2 cols2=4 \
	  slev=a,15,a,15 \
	  levs1=1,2,3,4,5,6,7,8,9,10,12,20,30,40,50 \
	  levs2=1,2,3,4,5,6,7,8,9,10,12,20,30,40,50 
goto end

  # Cleanup
    rm -r win15.xy
    rm -r win15.tmp 
    rm -r co.*.sl co.*.mp

  # ==============================
  # ====== SiO channel maps ====== 
  # ==============================
    
  # Copy phase vs time from wide.av, rewrite data to apply
	gpcopy vis=wide.av out=win7.xy options=nopass,nopol
	rm -r win7.tmp
    uvcat vis=win7.xy out=win7.tmp options=nopass,nopol
    uvplt vis=win7.tmp axis=time,phase device=/xs line=chan,1,1,191 average=10 \
	  select='-source('$SRC')' options=nocal,nopass,nopol nxy=2,2

  # Fit passband (and avg gain) to grand avg of win7 data
    mfcal vis=win7.tmp line=chan,10,1,19,19 select='source('$PASSCAL')' \
      interval=10000 refant=$REFANT
    smagpplt vis=win7.tmp options=bandpass,nofit,wrap device=/xs yrange=-180,180 \
      xaxis=chan yaxis=phase nxy=5,3
    smagpplt vis=win7.tmp options=bandpass,nofit,wrap device=/xs yrange=0,3 \
      xaxis=chan yaxis=amp nxy=5,3
	$<
	gpplt vis=win7.tmp yaxis=amp device=/xs yrange=0,3 nxy=5,3
	$<
	gpplt vis=win7.tmp yaxis=phase device=/xs yrange=-180,180 options=wrap nxy=5,3
	$<
goto end

  # Rewrite the data to apply the passband 
    rm -r win7.cal
    uvcat vis=win7.tmp out=win7.cal options=nopol

  # Copy amplitude gains vs time from wide.tmp, leakages from wide.lsb
    gpcopy vis=wide.tmp out=win7.cal options=nopass,nopol
	gpcopy vis=wide.lsb out=win7.cal options=nopass,nocal
    puthd in=win7.cal/senmodel value='GSV' type=ascii
    uvplt vis=win7.cal axis=time,phase device=/xs line=chan,1,1,191 average=10 \
	  select='-source('$SRC')'
    $<
    uvspec vis=win7.cal interval=10000 device=/xs axis=chan,phase yrange=-180,180 \
	  select='source('$PASSCAL')' nxy=4,4 line=chan,38,1,5,5

  # Make SiO channel maps
	puthd in=win7.cal/restfreq value=217.10498 type=double
	rm -r sio.I.mp sio.Q.mp sio.U.mp sio.V.mp sio.bm
	invert vis=win7.cal line=chan,38,1,5,5 \
	  map=sio.I.mp,sio.Q.mp,sio.U.mp,sio.V.mp beam=sio.bm stokes=I,Q,U,V sup=0 \
	  'select=source('$SRC')' options=systemp cell=0.25 imsize=512
	cgdisp in=sio.I.mp device=/xs region=arcsec,'box(-30,-30,30,30)' options=3value \
	  csize=0,1,1 nxy=3,3

	rm noiseList
	foreach MP (sio.I sio.Q sio.U sio.V)
	  rm -r $MP.sl
	  clean map=$MP.mp beam=$MAP.bm out=$MP.sl niters=1000
	  rm -r $MP.cm
	  restor map=$MP.mp beam=$MAP.bm model=$MP.sl out=$MP.cm
	  cgdisp in=$MP.cm device=/xs region=quarter labtyp=arcsec
	  echo " " >> noiseList
	  imlist options=stat in=$MP.cm region='arcsec,box(20,-20,-20,20)' log=junk
      cat junk >> noiseList
	  $<
	end
    view noiseList
goto end

# >>> edit SIO_I_SIGMA, SIO_QU_SIGMA <<< #

    cgdisp in=sio.I.cm,sio.Q.cm,sio.U.cm type=pixel,contour,contour \
      device=/xs region='arcsec,box(-20,-20,20,20)'  \
	  options=3pixel csize=0,1,1 range=$SIO_I_SIGMA,8,lin,1 \
      siols1=2 siols2=4 slev=a,$SIO_QU_SIGMA,a,$SIO_QU_SIGMA \
      levs1=-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,3,4,5,6,7,8,9,10,11,12,13,14,15 \
      levs2=-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,3,4,5,6,7,8,9,10,11,12,13,14,15 \
      nxy=2,2

  # Create red and blue outflow maps
	rm -r siored.cm
	moment in=sio.I.cm out=siored.cm mom=0 clip=-1000,.5 region='image(11,18)'
    cgdisp in=siored.cm device=/xs region=quarter range=6,100,lin,1

	rm -r sioblue.cm
	moment in=sio.I.cm out=sioblue.cm mom=0 clip=-1000,.5 region='image(22,36)'
    cgdisp in=sioblue.cm device=/xs region=quarter range=6,100,lin,1

	cgdisp device=/xs in=$MAP.I.cm,$MAP.poli.cm,$MAP.pa.cm,siored.cm,sioblue.cm \
      type=pixel,amp,angle,contour,contour \
	  region='arcsec,box(18,-14,-18,14)' options=full,rot90 labtyp=hms,dms \
	  vecfac=1.2,4,4 beamtyp=b,l,4 lines=1,1,1,10 siols1=2 siols2=4 \
	  slev=a,15,a,15 \
	  levs1=1,2,3,4,5,6,7,8,9,10,12,20,30,40,50 \
	  levs2=1,2,3,4,5,6,7,8,9,10,12,20,30,40,50 

  # Cleanup
    rm -r win7.xy
    rm -r win7.tmp
    rm -r sio.sl sio.mp

end:
