;+
; Description: Generates the wavelengths for genshkuratov
; Parameters:
;   c=c_list, list of concentrations for each endmember
;   q=q, the volume fraction filled with particles
;   s=s_list, list of sizes for each endmember
;   wave1-wave0, wavelengths for each of up to ten endmembers
;   real1-real0, real indices of refraction for each of up to ten endmembers
;   im1-im0, imaginary indices of refraction for each of up to ten endmembers
; Author: Szilard Gyalay
;-
function genshkuratovwaves, c=c_list, q=q, s=s_list, wave1=wave1_1, wave2=wave2_1, $
  wave3=wave3_1, wave4=wave4_1, wave5=wave5_1, wave6=wave6_1, wave7=wave7_1, $
  wave8=wave8_1, wave9=wave9_1, wave0=wave10_1, real1=real1_1, real2=real2_1, $
  real3=real3_1, real4=real4_1, real5=real5_1, real6=real6_1, real7=real7_1, $
  real8=real8_1, real9=real9_1, real0=real10_1, im1=im1_1, im2=im2_1, im3=im3_1, $
  im4=im4_1, im5=im5_1, im6=im6_1, im7=im7_1, im8=im8_1, im9=im9_1, im0=im10_1

  ;) To obtain the wavelengths for Shkuratov coarse mixture modeling
  ;) then can use genshkuratov to find the albedos
  ;) by finding n,k for several endmembers using these wavelengths
  
  ;) Parameters:
  ;) c_list, list of concentrations for each endmember
  ;) q, the volume fraction filled by particles
  ;) s_list, list of particle sizes for each endmember
  ;) wave#, array of wavelengths for each endmember optical constant dataset
  ;) real#, array of real indices of refraction for each wavelength
  ;) im#, as real but imaginary indices of refraction

  ;) At this moment, only supports ten endmembers

  ;) For simplicity, number of endmembers
  num=n_elements(c_list)

  ;) This is just to make sure that people stick parameters in
  if num eq 0 then print, 'Error: Concentrations not specified!'
  if n_elements(q) eq 0 then begin
    print, 'Error: Volume Fraction filled by particles not specified!'
  endif
  if n_elements(s_list) eq 0 then print, 'Error: Sizes not specified!'

  ;) Then make sure that c_list and s_list match up in size
  if num ne n_elements(s_list) then begin
    print, 'Error: Endmember numbers mismatch!'
  endif

  ;) Need to make sure the wavelength range accounts for all relevant datasets
  ;) Do this by finding the max first wavelength and min last wavelength
  ;) Then constrain the wavelength, real ind, and im ind lists to this range

  ;) put all the first wavelengths into an array
  wave_begins=dblarr(num)
  for i=1,num do begin
    if i eq 1 then wave_begins[0]=wave1_1[0]
    if i eq 2 then wave_begins[1]=wave2_1[0]
    if i eq 3 then wave_begins[2]=wave3_1[0]
    if i eq 4 then wave_begins[3]=wave4_1[0]
    if i eq 5 then wave_begins[4]=wave5_1[0]
    if i eq 6 then wave_begins[5]=wave6_1[0]
    if i eq 7 then wave_begins[6]=wave7_1[0]
    if i eq 8 then wave_begins[7]=wave8_1[0]
    if i eq 9 then wave_begins[8]=wave9_1[0]
    if i eq 10 then wave_begins[9]=wave10_1[0]
  endfor

  ;) Then the same for the last wavelengths
  wave_ends=dblarr(num)
  for i=1,num do begin
    if i eq 1 then wave_ends[0]=wave1_1[(n_elements(wave1_1)-1)]
    if i eq 2 then wave_ends[1]=wave2_1[(n_elements(wave2_1)-1)]
    if i eq 3 then wave_ends[2]=wave3_1[(n_elements(wave3_1)-1)]
    if i eq 4 then wave_ends[3]=wave4_1[(n_elements(wave4_1)-1)]
    if i eq 5 then wave_ends[4]=wave5_1[(n_elements(wave5_1)-1)]
    if i eq 6 then wave_ends[5]=wave6_1[(n_elements(wave6_1)-1)]
    if i eq 7 then wave_ends[6]=wave7_1[(n_elements(wave7_1)-1)]
    if i eq 8 then wave_ends[7]=wave8_1[(n_elements(wave8_1)-1)]
    if i eq 9 then wave_ends[8]=wave9_1[(n_elements(wave9_1)-1)]
    if i eq 10 then wave_ends[9]=wave10_1[(n_elements(wave10_1)-1)]
  endfor

  ;) Find max value of beginnings
  largest_begin=max(wave_begins)

  ;) Find min value of ends
  smallest_end=min(wave_ends)

  ;) Constrain each list to that amount
  for i=1,num do begin
    if i eq 1 then begin
      maxi=where(wave1_1 GE largest_begin)
      wave1_2=wave1_1[maxi]
      mini=where(wave1_2 LE smallest_end)
      wave1_3=wave1_2[mini]
    endif
    if i eq 2 then begin
      maxi=where(wave2_1 GE largest_begin)
      wave2_2=wave2_1[maxi]
      mini=where(wave2_2 LE smallest_end)
      wave2_3=wave2_2[mini]
    endif
    if i eq 3 then begin
      maxi=where(wave3_1 GE largest_begin)
      wave3_2=wave3_1[maxi]
      mini=where(wave3_2 LE smallest_end)
      wave3_3=wave3_2[mini]
    endif
    if i eq 4 then begin
      maxi=where(wave4_1 GE largest_begin)
      wave4_2=wave4_1[maxi]
      mini=where(wave4_2 LE smallest_end)
      wave4_3=wave4_2[mini]
    endif
    if i eq 5 then begin
      maxi=where(wave5_1 GE largest_begin)
      wave5_2=wave5_1[maxi]
      mini=where(wave5_2 LE smallest_end)
      wave5_3=wave5_2[mini]
    endif
    if i eq 6 then begin
      maxi=where(wave6_1 GE largest_begin)
      wave6_2=wave6_1[maxi]
      mini=where(wave6_2 LE smallest_end)
      wave6_3=wave6_2[mini]
    endif
    if i eq 7 then begin
      maxi=where(wave7_1 GE largest_begin)
      wave7_2=wave7_1[maxi]
      mini=where(wave7_2 LE smallest_end)
      wave7_3=wave7_2[mini]
    endif
    if i eq 8 then begin
      maxi=where(wave8_1 GE largest_begin)
      wave8_2=wave8_1[maxi]
      mini=where(wave8_2 LE smallest_end)
      wave8_3=wave8_2[mini]
    endif
    if i eq 9 then begin
      maxi=where(wave9_1 GE largest_begin)
      wave9_2=wave9_1[maxi]
      mini=where(wave9_2 LE smallest_end)
      wave9_3=wave9_2[mini]
    endif
    if i eq 10 then begin
      maxi=where(wave10_1 GE largest_begin)
      wave10_2=wave10_1[maxi]
      mini=where(wave10_2 LE smallest_end)
      wave10_3=wave10_2[mini]
    endif
  endfor

  ;) wavex_3 is only used to find the list with the best spectral resolution
  ;) Else, we'd cut off the ends needed later for finding the approximate value
  ;) of n,k at the wavelengths

  ;) Finding the highest spectral resolution (most wavelengths within constraint)
  res_list=dblarr(num)
  for i=1,num do begin
    if i eq 1 then res_list[0]=n_elements(wave1_3)
    if i eq 2 then res_list[1]=n_elements(wave2_3)
    if i eq 3 then res_list[2]=n_elements(wave3_3)
    if i eq 4 then res_list[3]=n_elements(wave4_3)
    if i eq 5 then res_list[4]=n_elements(wave5_3)
    if i eq 6 then res_list[5]=n_elements(wave6_3)
    if i eq 7 then res_list[6]=n_elements(wave7_3)
    if i eq 8 then res_list[7]=n_elements(wave8_3)
    if i eq 9 then res_list[8]=n_elements(wave9_3)
    if i eq 10 then res_list[9]=n_elements(wave10_3)
  endfor
  max_res=max(res_list)

  ;) Then, to get our final wavelength list
  if n_elements(wave1_3) eq max_res then lambda=wave1_3 else $
    if n_elements(wave2_3) eq max_res then lambda=wave2_3 else $
    if n_elements(wave3_3) eq max_res then lambda=wave3_3 else $
    if n_elements(wave4_3) eq max_res then lambda=wave4_3 else $
    if n_elements(wave5_3) eq max_res then lambda=wave5_3 else $
    if n_elements(wave6_3) eq max_res then lambda=wave6_3 else $
    if n_elements(wave7_3) eq max_res then lambda=wave7_3 else $
    if n_elements(wave8_3) eq max_res then lambda=wave8_3 else $
    if n_elements(wave9_3) eq max_res then lambda=wave9_3 else $
    if n_elements(wave10_3) eq max_res then lambda=wave10_3

  ;) now lambda is the list of wavelengths we'll have the rest of the data to
  
  return, lambda
  
end