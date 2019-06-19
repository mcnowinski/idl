
; SPECFIT.PRO
;
; PURPOSE:
;    Given an asteroid spectrum reflectance and wavelength file, this
;    will use stock olivine/pyroxene/chromium models to model the
;    abundance of each in the asteroid. This is compositional
;    modeling. Many factors such as grain size and space weathering
;    are also included. This is based on the Skuratov modeling
;    procedure and is therefore different than band center and area
;    analysis. Shuratov uses the physics of light and reflectance to
;    solve for the mixture of composition.
;
; INPUT:
;    - file containing columns of wavelength and reflectance for asteroid
;    - initial guess for fit (floats)
;          [oliv%, opx%, cpx%, chr%]
;    - starting parameters (hard wired):
;          [grain_size, weathering_parameter,opacity, porosity, $
;          normalization_wavelength]
;
; OUTPUT:
;    - one file containting the parameter solutions in "./data/"
;    - spectral plot with fit solution in "./plot/"
;
; NOTES:
;    Hard wired parameters may be chanaged. However, take care in carrying all
;    changes throughout the code. For example, chromium is kept at a fixed
;    precentage. If you wan to unfix it, you need to do this at the ".fixed"
;    location and also indicate the change in the function "equalOne".
;
; HISTORY;
;    2014 11 01 // bjb -- acquisition, lots of function building
;    2014 05 13 // bjb -- added probability of solution determination (findProb)
;    2014 05 16 // bjb -- added approximation function to findProb for finer
;                         probability determination
;    2014 06 10 // bjb -- removed dependency on setup function (for
;                         intuitiveness)
;                         added residual to plot. some plot printing changes.
;    2014 07 11 // bjb -- added functionality to fit a specific range only
;    2014 07 14 // bjb -- cleaning of code.
;###################################################################################

;----------------------------------------------------
; READ THE COMPOSITION FILES
;----------------------------------------------------
FUNCTION readComps,theFile,dir,cutoff,xinter

  ; read the file for parsing
  readcol, dir+theFile, wav, f1, f2, /silent

  ; Xrange dictates the interpolation
  npoints = n_elements(xinter)
  good  = where(wav LE cutoff)
  data = fltarr(3,npoints)

  data(0,*) = interpol(wav(good), wav(good), xinter)
  data(1,*) = interpol(f1(good), wav(good),xinter)
  data(2,*) = interpol(f2(good), wav(good), xinter)

  RETURN, data
END

;----------------------------------------------------
; THE SHKURATOV MATH: SEE SHKURATOV 1999
;----------------------------------------------------
FUNCTION shkuratov,theData,theLen,opac

  ;theData(0 is the wavelength
  ;theData(1 and the Data(2 are the optical constants (n and k), the real and imaginary parts of the refractive index, respectively
  ro = (theData(1,*)-1)^2/(theData(1,*)+1)^2 ; Fresnel coef
  Re = ro+0.05
  Ri = 1.04-1/(theData(1,*)^2)               ; intern ref coef
  Rb = (0.28*TheData(1,*)-0.2)*Re            ; bckwds ref coef
  Rf = Re-Rb                                 ; forwds ref coef
  Te = 1-Re
  Ti = Te/(theData(1,*)^2)

  t = (4*!Pi*theData(2,*)*opac)/theData(0,*) ; optical density

  r_b = Rb+0.5*Te*Ti*Ri*exp(-2*t(*))/(1-Ri*exp(-t(*)))
  r_f = Rf+Te*Ti*exp(-t(*))+0.5*Te*Ti*Ri*exp(-2*t(*))/(1-Ri*exp(-t(*)))

  RETURN, [r_b, r_f]
END

;----------------------------------------------------
; DETERMINE THE FLUX TO NORMALIZE AT
;----------------------------------------------------
FUNCTION getNorm, waves, yflux, normAt

  near    = min(abs(waves-normAt), index) ; find x near normat
  normFac = yflux[index]

  RETURN, normFac
END

;----------------------------------------------------
; ABILITY TO WEIGHT THE FIT
;----------------------------------------------------
FUNCTION findWeights, waves, lowbrac, midbrac, highbrac

  ; Determin how to weight the fit
  ; min ---> lowbrac <---> midbrac <---> hibrac <--- max

  ; find index of bracket point
  dum  = min(abs(waves-lowbrac), indLow)
  dum  = min(abs(waves-highbrac), indHig)
  dum  = min(abs(waves-midbrac), indMid)

  ; Fill list with regional weights (0=lw, 1=hg)
  lowWeight = replicate(0.0,n_elements(waves[0:indLow-1]))                  ; beg-->low
  midWeight = replicate(0.6, n_elements(waves[indMid+1:indHig-1]))            ; mid-->hig
  higWeight = replicate(1.0, n_elements(waves[indLow:indMid]))              ; low-->mid
  ignWeight = replicate(0.0, n_elements(waves[indHig:n_elements(waves)-1])) ; hig-->end

  weightArray = [lowWeight, higWeight, midWeight,ignWeight]

  RETURN, weightArray
END

;----------------------------------------------------
; MAKE SUM OF COMPONENTS EQUAL 100%
;----------------------------------------------------
FUNCTION equalOne, inVals, fixed

  ; determine distance from 100%, divide remainder amoung
  ; components and be sure to stay > 0 and < 100
  offBy = 1-total(inVals)
  addTo = offBy/(n_elements(inVals)-total(fixed))
  FOR j=0,n_elements(inVals)-1 DO BEGIN
     IF fixed[j] NE 1 AND inVals[j]+addTo GT 0 THEN BEGIN
        inVals[j] = inVals[j] + addTo
     ENDIF
  ENDFOR

  ; do again a few times if not getting to 100%
  totalInVals = total(inVals)
  IF abs(1-totalInVals) GT 0.002 THEN BEGIN
     FOR b=0,10 DO BEGIN
        offBy = 1-total(inVals)
        addTo = offBy/(n_elements(inVals)-total(fixed))
        FOR j=0,n_elements(inVals)-1 DO BEGIN
           IF fixed[j] NE 1 AND inVals[j]+addTo GT 0 THEN BEGIN
              inVals[j] = inVals[j] + addTo
           ENDIF
        ENDFOR
     ENDFOR
  ENDIF

  RETURN, inVals
END

;----------------------------------------------------
; Approximation function for probability determination
;----------------------------------------------------
FUNCTION calcProb, ratio

  ; given a ratio solution, determine it probability of type
  outProb = 50.*exp(-3.23946*exp(-2.18829*ratio^0.88103))

  RETURN, outProb

END

;----------------------------------------------------
; Function to determine probabilities of solution
;----------------------------------------------------
FUNCTION findProb, inSol

  ; define regional boundaries and error
  ;(found empirically when comparing to dunn
  HL_bound   = 0.5427
  LLL_bound = 0.615625
  stdNum    = 0.05

  ; deal with H solutions
  IF inSol LE HL_bound THEN BEGIN
     type = 'H'
     toHLbound  = calcProb((HL_bound - inSol)/stdNum)
     toLLLbound = calcProb((LLL_bound - inSol)/stdNum)

     Hprob  = 50. + toHLbound
     Lprob  = toLLLbound - toHLbound
     LLprob = 50. - toLLLbound
  ENDIF

  ; deal with LL solutions
  IF inSol GE LLL_bound THEN BEGIN
     type = 'LL'
     toHLbound  = calcProb((inSol - HL_bound)/stdNum)
     toLLLbound = calcProb((inSol - LLL_bound)/stdNum)

     LLprob = 50. + toLLLbound
     Lprob  = toHLbound - toLLLbound
     Hprob  = 50. - toHLbound
  ENDIF

  ; finally deal with L solutions
  IF inSol GT HL_bound AND inSol LT LLL_bound THEN BEGIN
     type = 'L'
     toHLbound  = calcProb((inSol - HL_bound)/stdNum)
     toLLLbound = calcProb((LLL_bound - inSol)/stdNum)

     Hprob  = 50. - toHLbound
     Lprob  = toHLbound + toLLLbound
     LLprob = 50. - toLLLbound
  ENDIF


  probs = {typeProbs:[Hprob, Lprob, LLprob], letterType:type}
  RETURN, probs
END

;----------------------------------------------------
;Prepare for the fitting
;----------------------------------------------------
FUNCTION FIT_FLUX, X, P

  COMMON FIT, files, fracs, fixed, num_c, dir

  ; upper limit of fit
  cutoff = max(X)

  ; Make the solutions equal one
  fixedValues = fixed       ; 1=fixed -- This must be manually flagged
  newP = equalOne(P[0:num_c-1], fixedValues)
  P[0:num_c-1] = newP

  q=P(num_c+3)                        ; the porosity
  opacity = P(num_c)*P(num_c+2)           ; the opacity (gs*opac)
  normLoc = P(num_c+4)                ; normalization point

  ; Parse the data and interlope appropriately
  comps = []
  FOR I = 0, num_c - 1 DO BEGIN
    comp = readComps(files(I), dir, cutoff, X)
    comps = [comps, comp]
  ENDFOR 

  len    = n_elements(comps(0,*))
  ;-------------
  ; Shkuratov
  backs = []
  FOR I = 0, num_c - 1 DO BEGIN
    i1 = I*3
    i2 = I*3+2
    _back = shkuratov(comps[i1:i2, *], len, opacity)
    backs = [backs, _back]
  ENDFOR

  _numerator_b = 0
  _numerator_f = 0
  _denominator = 0
  FOR I = 0, num_c - 1 DO BEGIN
    i1 = I*2
    i2 = I*2+1
    r_b = backs(i1, *)
    r_f = backs(i2, *)
    _denominator = _denominator + P(I)
    _numerator_b = _numerator_b + r_b*P(I)
    _numerator_f = _numerator_f + r_f*P(I)    
  ENDFOR

  pb = q * _numerator_b / _denominator
  pf = q * _numerator_f / _denominator + 1 - q

  ; the reflectance
  Fa  = [(1+pb^2-pf^2)/(2*pb)-[((1+pb^2-pf(*)^2)/(2*pb))^2-1]^0.5]
  ;-------------

  ; ---
  ; Normalization: yes if you ignore the albedo; no if you know it and
  ; try to use the albedo as an additional constraint
  ; ---

  ; weathering equation
  FOR b=0,len-1 DO BEGIN
     Fa(b) = exp(-1.*P(num_c+1)/comps(0,b))*Fa(b)
  ENDFOR

  ; normalize the solution
  anorm = getNorm(comps(0,*),Fa,normLoc)
  Fa = Fa/anorm
  RETURN, transpose(Fa)
END
;----------------------------------------------------
;####################################################
; SPECFIT.PRO
; Input a spectrum and initial guess values. Fit the
; spectrum using Shkuratov. Output a plot of the fit
; and a file of fit parameters.
;####################################################
;----------------------------------------------------
PRO SPECFIT, specfile, _files, _fracs, _fixed

   COMMON FIT, files, fracs, fixed, num_c, dir
   
   files = _files
   fracs = _fracs
   fixed = _fixed

   ; specfile is reflectance input
   ; files is list of files containing the optical constants (n and k) for each constituent mineral
   ; fracs is the initial guess at a fractional abundance of each constituent. This should add to 1.
   ; fixed defines (by setting to 1) which of the fracs values cannot be modified by the code 

   ;
   ; parameters
   ;
   ; turn on verbose output
   VERBOSE = 1  
   ;
   ; the files containing the optical constants for the constituent minerals
   ; these files need to be in the ./ocs directory (relative to this script)
   ;files=['oliv70.dat','opx70.dat','cpx70.dat','chromite.dat']
   dir = ROUTINE_DIR() + '../ocs/'
   
   ; number of constituents
   num_c = N_ELEMENTS(files)
   
   ; ensure the number of fracs equals the number of files
   IF N_ELEMENTS(fracs) NE num_c THEN BEGIN   
     PRINT, "Error. Number of fractional values should be ", N_ELEMENTS(files), "."
     STOP
   ENDIF  

   ; ensure the number of fixes equals the number of files
   IF N_ELEMENTS(fixed) NE num_c THEN BEGIN
     PRINT, "Error. Number of fixed values should be ", N_ELEMENTS(files), "."
     STOP
   ENDIF
   
   ; ensure fractional values add to 1
   IF TOTAL(fracs) NE 1 THEN BEGIN
     PRINT, "Error. Fractional values do not add to 1."
     STOP    
   ENDIF
   
   ; read in the asteroid spectrum and only take indexes that have non-zero reflectances
   ; file should look likes this:
   ; 0.45 0.718309859
   ; 0.503125 0.718309859
   ; 0.53125 0.757042254
   ; 0.55 0.85915493
   ; 0.58125 0.862676056
   ; ...
   ;readcol,specfile,xarrIN,yarrIN,format='F', skipline=2,/silent
   ;
   ;read in target reflectance values
   ;
   readcol,specfile,xarrIN,yarrIN,format='F',/silent
   good = where(yarrIN GT 0)
   yarr = yarrIN(good)
   xarr = xarrIN(good)

   ; read in the optical constant files
   cmpfile = files(0)
   ; open first optical constant file
   readcol,dir+cmpfile,interwv,f,f, format='F', /silent

   ; interlop the data to match resolution
   ; only use comparison file values that are between the Xmin and Xmax of the input asteroid spectra
   interwv = interwv(where(interwv GE min(xarr) AND interwv LE max(xarr)))
   len     = n_elements(interwv)
   ; interpolate reflectance data over the set of optical constant wavelengths
   yarr    = interpol(yarr, xarr, interwv)
   xarr    = interwv

   ; find the y value to normalize at
   mcrn1 = min(abs(xarr-0.9), indmcrn1)
   nrm = min(abs(yarr-max(yarr[0:indmcrn1])), indnrm)

   ; decide the initial guesses for fitting0
   normAtWav = xarr[indnrm]
   P = [fracs,150,0.3,1.0,0.5,normAtWav];xarr[indnrm]]
   IF indnrm EQ 0.0  THEN BEGIN           ; if no vis data
       P = [fracs,150,0.3,1.0,0.5,normAtWav]
   ENDIF

   ; normalize the data and chance shape
   normYarr = getNorm(xarr, yarr, normAtWav)
   X = transpose(xarr)
   Y = transpose(yarr/normYarr)

   ; make up error for the data (the weights are what matter)
   anArray = replicate(1,n_elements(Y))
   errs    = 0.0001 * anArray
   weighting = findWeights(X,normAtWav,1.2,2.0);P[8],1.2,2.0)

   ;----------------------------------------------
   ; MPfit - the first time
   ; set up limits of input parameters for fitting
   ;----------------------------------------------

   ; define structer for input parameters limits for fitting
   nelem   = n_elements(P)
   pi = replicate({value:0,fixed:0,limited:[1,0],limits:[0D,0D]},nelem)

   ; 0 = parameters are free
   ;pi(0).fixed = 0              ;olivine                    ;0
   ;pi(1).fixed = 0              ;orthopyroxene              ;0
   ;pi(2).fixed = 0              ;chlinopyroxene             ;1
   ;pi(3).fixed = 0              ;chromium                   ;1
   FOR I = 0, num_c - 1 DO BEGIN
    pi(I).fixed = fixed(I)
   ENDFOR
   
   pi(num_c).fixed = 0              ;grain size                 ;0
   pi(num_c+1).fixed = 0              ;weathering parameter       ;0
   pi(num_c+2).fixed = 1              ;opacity                    ;1
   pi(num_c+3).fixed = 1              ;porosity                   ;1
   pi(num_c+4).fixed = 0              ;normalization wavelength   ;0

   ; our materials must obide by the limits below
   pi.limited = [1,1]

   ; each "material" can vary from 0 to 100 percent of the total
   ; if it is not fixed
   ;pi(0).limits = [0.D,1.D]     ;olivine
   ;pi(1).limits = [0.D,1.D]     ;orthopyroxene
   ;pi(2).limits = [0.D,1.D]     ;chlinopyroxene
   ;pi(3).limits = [0.D,1.D]     ;chromium
   FOR I = 0, num_c - 1 DO BEGIN
     pi(I).limits = [0.D,1.D]
   ENDFOR

   ; Grain size 50:500   --> less than 50 the model breaks down
   pi(num_c).limits = [50.D,1000.D]
   ; weather param
   pi(num_c+1).limits = [-1.D,1.D]
   pi(num_c+2).limits = [1./pi(num_c).limits[0],10.D] ; [1-->10]
   ; porosity range
   pi(num_c+3).limits = [0D,1D]
   ; normalization range
   pi(num_c+4).limits = [min(X),max(X)]

   ; only fist between P[8[ and max of range
   ; why, why, why?
   ;toFit  = where(X GE P[8] AND X LE max(X)) ; 0.6-->2.0
   ;toFitX = X(toFit)
   ;toFitY = Y(toFit)
   toFitX = X
   toFitY = Y

   ; fit the spectrum
   FIRStsolution = mpfitfun('FIT_FLUX', toFitX, toFitY, errs, P, PARINFO=pi, $
                           WEIGHTS=weighting, QUIET=1)

   ;----------------------------------------------
   ; MPfit - second pass
   ; scale the fit for better fit
   ;----------------------------------------------
   ;pi(0).fixed = 1              ; olivine                  ;1
   ;pi(1).fixed = 1              ; orthopyroxene            ;1
   ;pi(2).fixed = 1              ; chlinopyroxene           ;1
   ;pi(3).fixed = 1              ; chromium                 ;1
   FOR I = 0, num_c - 1 DO BEGIN
     pi(I).fixed = fixed(I)
   ENDFOR   
   
   pi(num_c).fixed = 1              ; grain size               ;1
   pi(num_c+1).fixed = 1              ; weathering parameter     ;1
   pi(num_c+2).fixed = 0              ; opacity                  ;0
   pi(num_c+3).fixed = 0              ; porosity                 ;0
   pi(num_c+4).fixed = 1              ; normalization wavelength ;1

   ; fit the spectrum (again)
   pi(num_c+1).limits = [-10.D,10.D]  ; let weatherin parameter try more space
   pi(num_c+3).limits = [-10D,10D]    ; let porosity try more space
   solution = mpfitfun('FIT_FLUX', toFitX, toFitY, errs, FIRSTsolution, PARINFO=pi, $
                       WEIGHTS=weighting, QUIET=1)

   ;-------------
   ; These numbers are empirically determined using
   ; the comparison of solutions to the Dunn solutions
   ;-------------
   slopeAdjust = 1.0/0.66379384
   zeroAdjust  = -0.30764479

   ; Imperically adjust the final fraction and some statistics
   ; fracAdjust = (solution[0]/(solution[0]+solution[1]+solution[2]))*slopespecfitAdjust+zeroAdjust
   fracAdjust = (solution[0]/(solution[0]+solution[1]+solution[2]))*slopeAdjust+zeroAdjust
   
   ; determing the fit deviation from the data
   avgDev = total(abs(Y[indnrm:*]-FIT_FLUX(X[indnrm:*],solution)))/n_elements(Y[indnrm:*])
   ; Get the probability for each type
   solProbs = findProb(fracAdjust)

   ;-------------
   ; print the solutions to screen
   ;-------------
   ;IF keyword_set(VERBOSE) THEN BEGIN

      print, 'TYPE  =',solProbs.letterType,     format='(A7,4X,A)'
      print, 'FRAC  =',fracAdjust,              format='(A7,1X,F9.4,1X,A2,F6.4)'
      print, 'AVGDEV=',avgDev,                  format='(A7,1X,F9.4)'
      print, ''
      
      ;print, 'Oliv  =',solution[0],             format='(A7,1X,F9.4,1X,A2,1X,F6.4)'
      ;print, 'Opx   =',solution[1],             format='(A7,1X,F9.4,1X,A2,1X,F6.4)'
      ;print, 'Cpx   =',solution[2],             format='(A7,1X,F9.4,1X,A2,1X,F6.4)'
      ;print, 'Chro  =',solution[3],             format='(A7,1X,F9.4,1X,A2,1X,F6.4)'
      FOR I = 0, num_c - 1 DO BEGIN
        print, files[I],solution[I],             format='(A12,1X,F9.4,1X,A2,1X,F6.4)'
      ENDFOR

      print, 'Grain Size',solution[num_c],             format='(A12,1X,F9.4)'
      print, 'efGS',solution[num_c]*solution[num_c+2], format='(A12,1X,F9.4)'
      print, 'Opacity',solution[num_c+2],             format='(A12,1X,F9.4)'
      print, 'Porosity',FIRSTsolution[num_c+3],        format='(A12,1X,F9.4)'
      print, 'efPor =',solution[num_c+3],             format='(A12,1X,F9.4)'
      print, 'normAt=',solution[num_c+4],             format='(A12,1X,F9.4)'
      print, ''

   ;ENDIF

   ;----------------------------------------------
   ; Write To File
   ;----------------------------------------------
   ;the  = strsplit(specfile, '/', /extract)
   ;name = the[n_elements(the)-1]
   olper = fracAdjust

   ;; openw, lun50, 'data/normSols'+strtrim(normVal,2)+'/'+name+'.data', /get_lun
   ;openw, lun50, 'data/'+name+'_data.txt', /get_lun, /append
   ;printf, lun50, name, solution[0], solution[1], solution[2], solution[3], olper, $
   ;               solution[4], solution[5], FIRSTsolution[7], avgDev, solProbs.typeProbs[0], $
   ;               solProbs.typeProbs[1], solProbs.typeProbs[2], format='(A20,1x,F8.4,1x,F8.4,1x,F8.4,1x,F8.4,1x,F8.4,1x,F8.4,1x,F8.4,1x,F8.4,1x,F8.4,1X,F6.2,1X,F6.2,1X,F6.2)'

   ;free_lun, lun50

   ;----------------------------------------------
   ; Plot the data and fit
   ;----------------------------------------------
   ; Get the solution and the fit
   allX = findgen(2000.)/951. + 0.4
   longSol  = FIT_FLUX(allX,solution)
   shortSol = FIT_FLUX(X,solution)
   residual = (shortSol(*,0)-Y)
   residual = (0.25-mean(residual))+residual

   ; decide how to label the figure
   IF strmatch(cmpFile,'*7*') THEN num='70'
   IF strmatch(cmpFile,'*6*') THEN num='60'

   ; some plot set up
   loadct,39,/silent            ;
   SET_PLOT, 'WIN'
   ;device, filename='plot/'+name+'_'+num+'.eps'
   ;device, /bold, /color, /encapsulated

   ; define symbols;
   B = findgen(16)*(!PI*2/16.)
   usersym, cos(B), sin(B), /FILL
   pm = '+-'

   ; adjust for when olper > or < 1 or 0
   IF olper GT 1.0 THEN olper=1.0
   IF olper LT 0.0 THEN olper=0.0

   ; find the y range of the plot
   ytop = 2.
   IF max(Y) GT 2. THEN BEGIN
      WHILE max(Y) GT ytop DO BEGIN
         ytop = ytop + 0.5
      ENDWHILE
   ENDIF

   ; Set up plot window
   plot, [0.,0.], [0.,0.], xrange=[0.,2.5], yrange=[0,ytop], title=specfile, $
                           ytitle='Relative Reflectance', xtitle='Wavelength',$
	                   charthick=1, xthick=2, ythick=2, /nodata ;, /noerase

   ; plot
   ; oploterr, X,Y,weighting
   oplot, X, Y, psym=8, symsize=0.75                      ; data
   oplot, allX, longSol, thick=3, color=250 ; the fit
   oplot, X, residual                       ; residual


   ;----------------------------------------------
   ; Print fractional infomration to plot
   ;----------------------------------------------
   ;xyouts,0.95,0.4,string('GS: ',solution(4),            format='(A-4,1x,F6.2)'),charthick=2
   ;xyouts,0.95,0.3,string('POR: ',FIRSTsolution(7),      format='(A-4,1x,F6.2)'),charthick=2
   ;xyouts,0.95,0.2,string('OPC:',solution(6),            format='(A-4,1x,F6.2)'),charthick=2
   ;xyouts,0.95,ytop-0.22,string('AVGDEV:',avgDev,        format='(A7,1x,F5.2)'),charthick=1
   ;xyouts,1.5,0.4,string('eGS: ',solution(4)*solution(6),format='(A-5,1x,F6.2)'),charthick=2
   ;xyouts,1.5,0.3,string('ePOR:',solution(7),            format='(A-5,1x,F6.2)'),charthick=2
   ;xyouts,1.5,0.2,string('CS:  ',solution(5),            format='(A-5,1x,F6.2)'),charthick=2
   ;xyouts,1.5,ytop-0.22,string('OL/(OL+PX):',olper,       format='(A11,F5.2)'),charthick=3

   ; Legend for fit and data
   oplot, [0.1,0.1], [ytop-0.18,ytop-0.18], psym=8, symsize=0.75
   oplot, [0.08,0.12], [ytop-0.23, ytop-0.23], thick=2, color=250
   oplot, [0.08,0.12], [ytop-0.28, ytop-0.28]
   xyouts, 0.15, ytop-0.2, 'Data', charthick=1
   xyouts, 0.15, ytop-0.25, 'Fit', charthick=1
   xyouts, 0.15, ytop-0.30, 'Residual', charthick=1
   xyouts,0.15,ytop-0.35,string('Grain Size: ',solution(num_c),            format='(A-15,1x,F6.2)'),charthick=1
   xyouts,0.15,ytop-0.40,string('Porosity: ',FIRSTsolution(num_c+3),            format='(A-15,1x,F6.2)'),charthick=1   
   xyouts,0.15,ytop-0.45,string('Opacity: ',solution(num_c+2),            format='(A-15,1x,F6.2)'),charthick=1
   FOR I = 0, num_c - 1 DO BEGIN
     xyouts,0.15,ytop-0.50-0.05*I,string(files(I)+': %',solution(I)*100, format='(A-15,F6.2,A2,F4.2)'),charthick=1
   ENDFOR

   ; peace y'all
   ;device,/close
   ;set_plot,'X'
   set_plot,'WIN'
   
END
