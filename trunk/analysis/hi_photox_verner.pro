function hi_photox_verner, eV

  sigma=eV

  Eth = 13.6d0
  Emax = 5.0d4
  E0 = 4.298d-1
  sig0 = 5.475d4
  ya = 3.288d1
  P = 2.963

  x = eV / E0
  y = x
  
  sigma = sig0 * (x-1)^2 * y^(0.5d0 * P - 5.5d0) * (1 + sqrt(y/ya))^(-P)
  sigma = sigma * 1.0d-18

  return, sigma
  






end
