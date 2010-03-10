function hi_photox_quick, ryd

  sigma=ryd
  eps=ryd
  num=ryd
  den=ryd

  SIGMA0_HI = 6.30d-18

  n = n_elements(ryd)
  iinvalid = where( ryd lt 1, badcnt )
  
  if badcnt ne 0 then begin
     print, 'energies below 1 Rydberg not allowed in calls to hi_photox'
     sigma=-1
     return, sigma
  endif else begin
     sigma = SIGMA0_HI * ryd^(-3.d0)
     return, sigma
  endelse






end





