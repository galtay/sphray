function hi_photox_analytic, ryd

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
     ieq1 = where(ryd eq 1, ieq1cnt, complement=igt1)
     if ieq1cnt ne 0 then sigma[ieq1] = SIGMA0_HI
     eps[igt1] = sqrt( ryd[igt1] - 1.0d0 )
     den[igt1] = 1.0d0 - exp( -2*!DPI / eps[igt1] )
     num[igt1] = exp( 4.d0 - 4.d0 * atan(eps[igt1]) / eps[igt1] )
     sigma[igt1] = SIGMA0_HI * ryd[igt1]^(-4) * num[igt1] / den[igt1]
     return, sigma
  endelse






end





