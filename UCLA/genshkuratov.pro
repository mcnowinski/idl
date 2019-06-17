;+
; Description: Calculates the albedo of a coarse mixture using Shkuratov et al. 1999
; Parameters:
;   c=c_list, list of concentrations for each endmember
;   q=q, the volume fraction filled with particles
;   s=s_list, list of sizes for each endmember
;   wave1-wave0, wavelengths for each of up to ten endmembers
;   real1-real0, real indices of refraction for each of up to ten endmembers
;   im1-im0, imaginary indices of refraction for each of up to ten endmembers
; Author: Szilard Gyalay
;-
function genshkuratov, c=c_list, q=q, s=s_list, wave1=wave1_1, wave2=wave2_1, $
  wave3=wave3_1, wave4=wave4_1, wave5=wave5_1, wave6=wave6_1, wave7=wave7_1, $
  wave8=wave8_1, wave9=wave9_1, wave0=wave10_1, real1=real1_1, real2=real2_1, $
  real3=real3_1, real4=real4_1, real5=real5_1, real6=real6_1, real7=real7_1, $
  real8=real8_1, real9=real9_1, real0=real10_1, im1=im1_1, im2=im2_1, im3=im3_1, $
  im4=im4_1, im5=im5_1, im6=im6_1, im7=im7_1, im8=im8_1, im9=im9_1, im0=im10_1
  
  ;) This is an attempt to make coarse intramixtures 
  ;) using Shkuratov et al. 1999 (S99)
  
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
  
  ;) Using another function, we can find the wavelengths for n,k to correspond to
  ;) Need two functions to return both the wavelengths and albedos for plotting
  lambda=genshkuratovwaves(c=c_list, q=q, s=s_list, wave1=wave1_1, wave2=wave2_1, $
    wave3=wave3_1, wave4=wave4_1, wave5=wave5_1, wave6=wave6_1, wave7=wave7_1, $
    wave8=wave8_1, wave9=wave9_1, wave0=wave10_1, real1=real1_1, real2=real2_1, $
    real3=real3_1, real4=real4_1, real5=real5_1, real6=real6_1, real7=real7_1, $
    real8=real8_1, real9=real9_1, real0=real10_1, im1=im1_1, im2=im2_1, im3=im3_1, $
    im4=im4_1, im5=im5_1, im6=im6_1, im7=im7_1, im8=im8_1, im9=im9_1, im0=im10_1)
  
  ;) Here we find the values of n,k for the wavelengths in lambda
  ;) Basically, for each endmember, we cycle through lambda.
  ;) We then find the wavelengths for each endmember before and after the lambda
  ;) Depending on how far lambda is from the wavelength before compared to after
  ;) we can find how much n,k should change. 
  for i=1,num do begin
    if i eq 1 then begin
      real1_4=dblarr(n_elements(lambda))
      im1_4=dblarr(n_elements(lambda))
      for j=0,(n_elements(lambda)-1) do begin
        wavelen=lambda[j]
        findnum=where(wave1_1 GE wavelen)
        if wavelen eq wave1_1[findnum[0]] then begin
          real1_4[j]=real1_1[findnum[0]]
          im1_4[j]=im1_1[findnum[0]]
        endif else begin
        wavebefore=wave1_1[findnum[0]-1]
        waveafter=wave1_1[findnum[0]]
        realbefore=real1_1[findnum[0]-1]
        realafter=real1_1[findnum[0]]
        imbefore=im1_1[findnum[0]-1]
        imafter=im1_1[findnum[0]]
        real1_4[j]=realbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(realafter-realbefore)
        im1_4[j]=imbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(imafter-imbefore)
        endelse
      endfor
    endif
    if i eq 2 then begin
      real2_4=dblarr(n_elements(lambda))
      im2_4=dblarr(n_elements(lambda))
      for j=0,(n_elements(lambda)-1) do begin
        wavelen=lambda[j]
        findnum=where(wave2_1 GE wavelen)
        if wavelen eq wave2_1[findnum[0]] then begin
          real2_4[j]=real2_1[findnum[0]]
          im2_4[j]=im2_1[findnum[0]]
        endif else begin
        wavebefore=wave2_1[findnum[0]-1]
        waveafter=wave2_1[findnum[0]]
        realbefore=real2_1[findnum[0]-1]
        realafter=real2_1[findnum[0]]
        imbefore=im2_1[findnum[0]-1]
        imafter=im2_1[findnum[0]]
        real2_4[j]=realbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(realafter-realbefore)
        im2_4[j]=imbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(imafter-imbefore)
        endelse
      endfor
    endif
    if i eq 3 then begin
      real3_4=dblarr(n_elements(lambda))
      im3_4=dblarr(n_elements(lambda))
      for j=0,(n_elements(lambda)-1) do begin
        wavelen=lambda[j]
        findnum=where(wave3_1 GE wavelen)
        if wavelen eq wave3_1[findnum[0]] then begin
          real3_4[j]=real3_1[findnum[0]]
          im3_4[j]=im3_1[findnum[0]]
        endif else begin
        wavebefore=wave3_1[findnum[0]-1]
        waveafter=wave3_1[findnum[0]]
        realbefore=real3_1[findnum[0]-1]
        realafter=real3_1[findnum[0]]
        imbefore=im3_1[findnum[0]-1]
        imafter=im3_1[findnum[0]]
        real3_4[j]=realbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(realafter-realbefore)
        im3_4[j]=imbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(imafter-imbefore)
        endelse
      endfor
    endif
    if i eq 4 then begin
      real4_4=dblarr(n_elements(lambda))
      im4_4=dblarr(n_elements(lambda))
      for j=0,(n_elements(lambda)-1) do begin
        wavelen=lambda[j]
        findnum=where(wave4_1 GE wavelen)
        if wavelen eq wave4_1[findnum[0]] then begin
          real4_4[j]=real4_1[findnum[0]]
          im4_4[j]=im4_1[findnum[0]]
        endif else begin
        wavebefore=wave4_1[findnum[0]-1]
        waveafter=wave4_1[findnum[0]]
        realbefore=real4_1[findnum[0]-1]
        realafter=real4_1[findnum[0]]
        imbefore=im4_1[findnum[0]-1]
        imafter=im4_1[findnum[0]]
        real4_4[j]=realbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(realafter-realbefore)
        im4_4[j]=imbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(imafter-imbefore)
        endelse
      endfor
    endif
    if i eq 5 then begin
      real5_4=dblarr(n_elements(lambda))
      im5_4=dblarr(n_elements(lambda))
      for j=0,(n_elements(lambda)-1) do begin
        wavelen=lambda[j]
        findnum=where(wave5_1 GE wavelen)
        if wavelen eq wave5_1[findnum[0]] then begin
          real5_4[j]=real5_1[findnum[0]]
          im5_4[j]=im5_1[findnum[0]]
        endif else begin
        wavebefore=wave5_1[findnum[0]-1]
        waveafter=wave5_1[findnum[0]]
        realbefore=real5_1[findnum[0]-1]
        realafter=real5_1[findnum[0]]
        imbefore=im5_1[findnum[0]-1]
        imafter=im5_1[findnum[0]]
        real5_4[j]=realbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(realafter-realbefore)
        im5_4[j]=imbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(imafter-imbefore)
        endelse
      endfor
    endif
    if i eq 6 then begin
      real6_4=dblarr(n_elements(lambda))
      im6_4=dblarr(n_elements(lambda))
      for j=0,(n_elements(lambda)-1) do begin
        wavelen=lambda[j]
        findnum=where(wave6_1 GE wavelen)
        if wavelen eq wave6_1[findnum[0]] then begin
          real6_4[j]=real6_1[findnum[0]]
          im6_4[j]=im6_1[findnum[0]]
        endif else begin
        wavebefore=wave6_1[findnum[0]-1]
        waveafter=wave6_1[findnum[0]]
        realbefore=real6_1[findnum[0]-1]
        realafter=real6_1[findnum[0]]
        imbefore=im6_1[findnum[0]-1]
        imafter=im6_1[findnum[0]]
        real6_4[j]=realbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(realafter-realbefore)
        im6_4[j]=imbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(imafter-imbefore)
        endelse
      endfor
    endif
    if i eq 7 then begin
      real7_4=dblarr(n_elements(lambda))
      im7_4=dblarr(n_elements(lambda))
      for j=0,(n_elements(lambda)-1) do begin
        wavelen=lambda[j]
        findnum=where(wave7_1 GE wavelen)
        if wavelen eq wave7_1[findnum[0]] then begin
          real7_4[j]=real7_1[findnum[0]]
          im7_4[j]=im7_1[findnum[0]]
        endif else begin
        wavebefore=wave7_1[findnum[0]-1]
        waveafter=wave7_1[findnum[0]]
        realbefore=real7_1[findnum[0]-1]
        realafter=real7_1[findnum[0]]
        imbefore=im7_1[findnum[0]-1]
        imafter=im7_1[findnum[0]]
        real7_4[j]=realbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(realafter-realbefore)
        im7_4[j]=imbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(imafter-imbefore)
        endelse
      endfor
    endif
    if i eq 8 then begin
      real8_4=dblarr(n_elements(lambda))
      im8_4=dblarr(n_elements(lambda))
      for j=0,(n_elements(lambda)-1) do begin
        wavelen=lambda[j]
        findnum=where(wave8_1 GE wavelen)
        if wavelen eq wave8_1[findnum[0]] then begin
          real8_4[j]=real8_1[findnum[0]]
          im8_4[j]=im8_1[findnum[0]]
        endif else begin
        wavebefore=wave8_1[findnum[0]-1]
        waveafter=wave8_1[findnum[0]]
        realbefore=real8_1[findnum[0]-1]
        realafter=real8_1[findnum[0]]
        imbefore=im8_1[findnum[0]-1]
        imafter=im8_1[findnum[0]]
        real8_4[j]=realbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(realafter-realbefore)
        im8_4[j]=imbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(imafter-imbefore)
        endelse
      endfor
    endif
    if i eq 9 then begin
      real9_4=dblarr(n_elements(lambda))
      im9_4=dblarr(n_elements(lambda))
      for j=0,(n_elements(lambda)-1) do begin
        wavelen=lambda[j]
        findnum=where(wave9_1 GE wavelen)
        if wavelen eq wave9_1[findnum[0]] then begin
          real9_4[j]=real9_1[findnum[0]]
          im9_4[j]=im9_1[findnum[0]]
        endif else begin
        wavebefore=wave9_1[findnum[0]-1]
        waveafter=wave9_1[findnum[0]]
        realbefore=real9_1[findnum[0]-1]
        realafter=real9_1[findnum[0]]
        imbefore=im9_1[findnum[0]-1]
        imafter=im9_1[findnum[0]]
        real9_4[j]=realbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(realafter-realbefore)
        im9_4[j]=imbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(imafter-imbefore)
        endelse
      endfor
    endif
    if i eq 10 then begin
      real10_4=dblarr(n_elements(lambda))
      im10_4=dblarr(n_elements(lambda))
      for j=0,(n_elements(lambda)-1) do begin
        wavelen=lambda[j]
        findnum=where(wave10_1 GE wavelen)
        if wavelen eq wave10_1[findnum[0]] then begin
          real10_4[j]=real10_1[findnum[0]]
          im10_4[j]=im10_1[findnum[0]]
        endif else begin
        wavebefore=wave10_1[findnum[0]-1]
        waveafter=wave10_1[findnum[0]]
        realbefore=real10_1[findnum[0]-1]
        realafter=real10_1[findnum[0]]
        imbefore=im10_1[findnum[0]-1]
        imafter=im10_1[findnum[0]]
        real10_4[j]=realbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(realafter-realbefore)
        im10_4[j]=imbefore+ $
          ((wavelen-wavebefore)/(waveafter-wavebefore))*(imafter-imbefore)
        endelse
      endfor
    endif
  endfor
  
  ;) Now we have n,k for all endmembers at the same wavelengths
  ;) So now we can start the program in ernest.
  
  ;) First, we use S99 eqs 7 to find reflection and transmission coefficients
  ;) These equations work for n=1.2-1.7
  ;) Also fractions of light flux scattered by a particle into the
  ;) backward (r_b) and forward (r_f) hemispheres, S99 eqs 9
  For i=1,num do begin
    if i eq 1 then begin
      r_0_1=((real1_4-1.)^2)/((real1_4+1.)^2)
      rbig_e_1=r_0_1+0.05
      rbig_b_1=(0.28*real1_4-.2)*rbig_e_1
      rbig_f_1=rbig_e_1-rbig_b_1
      rbig_i_1=1.04-1./(real1_4^2)
      T_e_1=1.-rbig_e_1
      T_i_1=1.-rbig_i_1
      tau_1=4.*!pi*im1_4*S_list[0]/lambda
      r_b_1=rbig_b_1+0.5*T_e_1*T_i_1*rbig_i_1*exp(-2.*tau_1)/ $
        (1.-rbig_i_1*exp(-tau_1))
      r_f_1=rbig_f_1+T_e_1*T_i_1*exp(-tau_1)+ $
        0.5*T_e_1*T_i_1*rbig_i_1*exp(-2.*tau_1)/(1.-rbig_i_1*exp(-tau_1))
    endif
    if i eq 2 then begin
      r_0_2=((real2_4-1.)^2)/((real2_4+1.)^2)
      rbig_e_2=r_0_2+0.05
      rbig_b_2=(0.28*real2_4-.2)*rbig_e_2
      rbig_f_2=rbig_e_2-rbig_b_2
      rbig_i_2=1.04-1./(real2_4^2)
      T_e_2=1.-rbig_e_2
      T_i_2=1.-rbig_i_2
      tau_2=4.*!pi*im2_4*S_list[1]/lambda
      r_b_2=rbig_b_2+0.5*T_e_2*T_i_2*rbig_i_2*exp(-2.*tau_2)/ $
        (1.-rbig_i_2*exp(-tau_2))
      r_f_2=rbig_f_2+T_e_2*T_i_2*exp(-tau_2)+ $
        0.5*T_e_2*T_i_2*rbig_i_2*exp(-2.*tau_2)/(1.-rbig_i_2*exp(-tau_2))
    endif  
    if i eq 3 then begin
      r_0_3=((real3_4-1.)^2)/((real3_4+1.)^2)
      rbig_e_3=r_0_3+0.05
      rbig_b_3=(0.28*real3_4-.2)*rbig_e_3
      rbig_f_3=rbig_e_3-rbig_b_3
      rbig_i_3=1.04-1./(real3_4^2)
      T_e_3=1.-rbig_e_3
      T_i_3=1.-rbig_i_3
      tau_3=4.*!pi*im3_4*S_list[2]/lambda
      r_b_3=rbig_b_3+0.5*T_e_3*T_i_3*rbig_i_3*exp(-2.*tau_3)/ $
        (1.-rbig_i_3*exp(-tau_3))
      r_f_3=rbig_f_3+T_e_3*T_i_3*exp(-tau_3)+ $
        0.5*T_e_3*T_i_3*rbig_i_3*exp(-2.*tau_3)/(1.-rbig_i_3*exp(-tau_3))
    endif
    if i eq 4 then begin
      r_0_4=((real4_4-1.)^2)/((real4_4+1.)^2)
      rbig_e_4=r_0_4+0.05
      rbig_b_4=(0.28*real4_4-.2)*rbig_e_4
      rbig_f_4=rbig_e_4-rbig_b_4
      rbig_i_4=1.04-1./(real4_4^2)
      T_e_4=1.-rbig_e_4
      T_i_4=1.-rbig_i_4
      tau_4=4.*!pi*im4_4*S_list[3]/lambda
      r_b_4=rbig_b_4+0.5*T_e_4*T_i_4*rbig_i_4*exp(-2.*tau_4)/ $
        (1.-rbig_i_4*exp(-tau_4))
      r_f_4=rbig_f_4+T_e_4*T_i_4*exp(-tau_4)+ $
        0.5*T_e_4*T_i_4*rbig_i_4*exp(-2.*tau_4)/(1.-rbig_i_4*exp(-tau_4))
    endif
    if i eq 5 then begin
      r_0_5=((real5_4-1.)^2)/((real5_4+1.)^2)
      rbig_e_5=r_0_5+0.05
      rbig_b_5=(0.28*real5_4-.2)*rbig_e_5
      rbig_f_5=rbig_e_5-rbig_b_5
      rbig_i_5=1.04-1./(real5_4^2)
      T_e_5=1.-rbig_e_5
      T_i_5=1.-rbig_i_5
      tau_5=4.*!pi*im5_4*S_list[4]/lambda
      r_b_5=rbig_b_5+0.5*T_e_5*T_i_5*rbig_i_5*exp(-2.*tau_5)/ $
        (1.-rbig_i_5*exp(-tau_5))
      r_f_5=rbig_f_5+T_e_5*T_i_5*exp(-tau_5)+ $
        0.5*T_e_5*T_i_5*rbig_i_5*exp(-2.*tau_5)/(1.-rbig_i_5*exp(-tau_5))
    endif
    if i eq 6 then begin
      r_0_6=((real6_4-1.)^2)/((real6_4+1.)^2)
      rbig_e_6=r_0_6+0.05
      rbig_b_6=(0.28*real6_4-.2)*rbig_e_6
      rbig_f_6=rbig_e_6-rbig_b_6
      rbig_i_6=1.04-1./(real6_4^2)
      T_e_6=1.-rbig_e_6
      T_i_6=1.-rbig_i_6
      tau_6=4.*!pi*im6_4*S_list[5]/lambda
      r_b_6=rbig_b_6+0.5*T_e_6*T_i_6*rbig_i_6*exp(-2.*tau_6)/ $
        (1.-rbig_i_6*exp(-tau_6))
      r_f_6=rbig_f_6+T_e_6*T_i_6*exp(-tau_6)+ $
        0.5*T_e_6*T_i_6*rbig_i_6*exp(-2.*tau_6)/(1.-rbig_i_6*exp(-tau_6))
    endif
    if i eq 7 then begin
      r_0_7=((real7_4-1.)^2)/((real7_4+1.)^2)
      rbig_e_7=r_0_7+0.05
      rbig_b_7=(0.28*real7_4-.2)*rbig_e_7
      rbig_f_7=rbig_e_7-rbig_b_7
      rbig_i_7=1.04-1./(real7_4^2)
      T_e_7=1.-rbig_e_7
      T_i_7=1.-rbig_i_7
      tau_7=4.*!pi*im7_4*S_list[6]/lambda
      r_b_7=rbig_b_7+0.5*T_e_7*T_i_7*rbig_i_7*exp(-2.*tau_7)/ $
        (1.-rbig_i_7*exp(-tau_7))
      r_f_7=rbig_f_7+T_e_7*T_i_7*exp(-tau_7)+ $
        0.5*T_e_7*T_i_7*rbig_i_7*exp(-2.*tau_7)/(1.-rbig_i_7*exp(-tau_7))
    endif
    if i eq 8 then begin
      r_0_8=((real8_4-1.)^2)/((real8_4+1.)^2)
      rbig_e_8=r_0_8+0.05
      rbig_b_8=(0.28*real8_4-.2)*rbig_e_8
      rbig_f_8=rbig_e_8-rbig_b_8
      rbig_i_8=1.04-1./(real8_4^2)
      T_e_8=1.-rbig_e_8
      T_i_8=1.-rbig_i_8
      tau_8=4.*!pi*im8_4*S_list[7]/lambda
      r_b_8=rbig_b_8+0.5*T_e_8*T_i_8*rbig_i_8*exp(-2.*tau_8)/ $
        (1.-rbig_i_8*exp(-tau_8))
      r_f_8=rbig_f_8+T_e_8*T_i_8*exp(-tau_8)+ $
        0.5*T_e_8*T_i_8*rbig_i_8*exp(-2.*tau_8)/(1.-rbig_i_8*exp(-tau_8))
    endif
    if i eq 9 then begin
      r_0_9=((real9_4-1.)^2)/((real9_4+1.)^2)
      rbig_e_9=r_0_9+0.05
      rbig_b_9=(0.28*real9_4-.2)*rbig_e_9
      rbig_f_9=rbig_e_9-rbig_b_9
      rbig_i_9=1.04-1./(real9_4^2)
      T_e_9=1.-rbig_e_9
      T_i_9=1.-rbig_i_9
      tau_9=4.*!pi*im9_4*S_list[8]/lambda
      r_b_9=rbig_b_9+0.5*T_e_9*T_i_9*rbig_i_9*exp(-2.*tau_9)/ $
        (1.-rbig_i_9*exp(-tau_9))
      r_f_9=rbig_f_9+T_e_9*T_i_9*exp(-tau_9)+ $
        0.5*T_e_9*T_i_9*rbig_i_9*exp(-2.*tau_9)/(1.-rbig_i_9*exp(-tau_9))
    endif
    if i eq 10 then begin
      r_0_10=((real10_4-1.)^2)/((real10_4+1.)^2)
      rbig_e_10=r_0_10+0.05
      rbig_b_10=(0.28*real10_4-.2)*rbig_e_10
      rbig_f_10=rbig_e_10-rbig_b_10
      rbig_i_10=1.04-1./(real10_4^2)
      T_e_10=1.-rbig_e_10
      T_i_10=1.-rbig_i_10
      tau_10=4.*!pi*im10_4*S_list[9]/lambda
      r_b_10=rbig_b_10+0.5*T_e_10*T_i_10*rbig_i_10*exp(-2.*tau_10)/ $
        (1.-rbig_i_10*exp(-tau_10))
      r_f_10=rbig_f_10+T_e_10*T_i_10*exp(-tau_10)+ $
        0.5*T_e_10*T_i_10*rbig_i_10*exp(-2.*tau_10)/(1.-rbig_i_10*exp(-tau_10))
    endif
  endfor
  
  ;) S99 eqs 14, 12
  case num of
    1: begin
        rho_b=q*(c_list[0]*r_b_1)
        rho_f=q*(c_list[0]*r_f_1)+1.-q
        thing=(1.+rho_b^2-rho_f^2)/(2*rho_b)
        albedo=thing-sqrt(thing^2-1.)
      end
    2: begin
        rho_b=q*(c_list[0]*r_b_1+c_list[1]*r_b_2)
        rho_f=q*(c_list[0]*r_f_1+c_list[1]*r_f_2)+1.-q
        thing=(1.+rho_b^2-rho_f^2)/(2*rho_b)
        albedo=thing-sqrt(thing^2-1.)
      end
    3: begin
        rho_b=q*(c_list[0]*r_b_1+c_list[1]*r_b_2+c_list[2]*r_b_3)
        rho_f=q*(c_list[0]*r_f_1+c_list[1]*r_f_2+c_list[2]*r_f_3)+1.-q
        thing=(1.+rho_b^2-rho_f^2)/(2*rho_b)
        albedo=thing-sqrt(thing^2-1.)
      end
    4: begin
        rho_b=q*(c_list[0]*r_b_1+c_list[1]*r_b_2+c_list[2]*r_b_3+ $
          c_list[3]*r_b_4)
        rho_f=q*(c_list[0]*r_f_1+c_list[1]*r_f_2+c_list[2]*r_f_3+ $
          c_list[3]*r_f_4)+1.-q
        thing=(1.+rho_b^2-rho_f^2)/(2*rho_b)
        albedo=thing-sqrt(thing^2-1.)
      end
    5: begin
        rho_b=q*(c_list[0]*r_b_1+c_list[1]*r_b_2+c_list[2]*r_b_3+ $
          c_list[3]*r_b_4+c_list[4]*r_b_5)
        rho_f=q*(c_list[0]*r_f_1+c_list[1]*r_f_2+c_list[2]*r_f_3+ $
          c_list[3]*r_f_4+c_list[4]*r_f_5)+1.-q
        thing=(1.+rho_b^2-rho_f^2)/(2*rho_b)
        albedo=thing-sqrt(thing^2-1.)
      end
    6: begin
        rho_b=q*(c_list[0]*r_b_1+c_list[1]*r_b_2+c_list[2]*r_b_3+ $
          c_list[3]*r_b_4+c_list[4]*r_b_5+c_list[5]*r_b_6)
        rho_f=q*(c_list[0]*r_f_1+c_list[1]*r_f_2+c_list[2]*r_f_3+ $
          c_list[3]*r_f_4+c_list[4]*r_f_5+c_list[5]*r_f_6)+1.-q
        thing=(1.+rho_b^2-rho_f^2)/(2*rho_b)
        albedo=thing-sqrt(thing^2-1.)
      end
    7: begin
        rho_b=q*(c_list[0]*r_b_1+c_list[1]*r_b_2+c_list[2]*r_b_3+ $
          c_list[3]*r_b_4+c_list[4]*r_b_5+c_list[5]*r_b_6+c_list[6]*r_b_7)
        rho_f=q*(c_list[0]*r_f_1+c_list[1]*r_f_2+c_list[2]*r_f_3+ $
          c_list[3]*r_f_4+c_list[4]*r_f_5+c_list[5]*r_f_6+c_list[6]*r_f_7+ $
          0)+1.-q
        thing=(1.+rho_b^2-rho_f^2)/(2*rho_b)
        albedo=thing-sqrt(thing^2-1.)
      end
    8: begin
        rho_b=q*(c_list[0]*r_b_1+c_list[1]*r_b_2+c_list[2]*r_b_3+ $
          c_list[3]*r_b_4+c_list[4]*r_b_5+c_list[5]*r_b_6+c_list[6]*r_b_7+ $
          c_list[7]*r_b_8)
        rho_f=q*(c_list[0]*r_f_1+c_list[1]*r_f_2+c_list[2]*r_f_3+ $
          c_list[3]*r_f_4+c_list[4]*r_f_5+c_list[5]*r_f_6+c_list[6]*r_f_7+ $
          c_list[7]*r_f_8)+1.-q
        thing=(1.+rho_b^2-rho_f^2)/(2*rho_b)
        albedo=thing-sqrt(thing^2-1.)
      end
    9: begin
        rho_b=q*(c_list[0]*r_b_1+c_list[1]*r_b_2+c_list[2]*r_b_3+ $
          c_list[3]*r_b_4+c_list[4]*r_b_5+c_list[5]*r_b_6+c_list[6]*r_b_7+ $
          c_list[7]*r_b_8+c_list[8]*r_b_9)
        rho_f=q*(c_list[0]*r_f_1+c_list[1]*r_f_2+c_list[2]*r_f_3+ $
          c_list[3]*r_f_4+c_list[4]*r_f_5+c_list[5]*r_f_6+c_list[6]*r_f_7+ $
          c_list[7]*r_f_8+c_list[8]*r_f_9)+1.-q
        thing=(1.+rho_b^2-rho_f^2)/(2*rho_b)
        albedo=thing-sqrt(thing^2-1.)
      end
    10: begin
        rho_b=q*(c_list[0]*r_b_1+c_list[1]*r_b_2+c_list[2]*r_b_3+ $
          c_list[3]*r_b_4+c_list[4]*r_b_5+c_list[5]*r_b_6+c_list[6]*r_b_7+ $
          c_list[7]*r_b_8+c_list[8]*r_b_9+c_list[9]*r_b_10)
        rho_f=q*(c_list[0]*r_f_1+c_list[1]*r_f_2+c_list[2]*r_f_3+ $
          c_list[3]*r_f_4+c_list[4]*r_f_5+c_list[5]*r_f_6+c_list[6]*r_f_7+ $
          c_list[7]*r_f_8+c_list[8]*r_f_9+c_list[9]*r_f_10)+1.-q
        thing=(1.+rho_b^2-rho_f^2)/(2*rho_b)
        albedo=thing-sqrt(thing^2-1.)
      end
  endcase
  
  return, albedo
end
  