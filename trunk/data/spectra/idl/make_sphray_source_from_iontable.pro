; This idl procedure reads an Ionization Table in HDF5 format,
; extracts the background UV spectra and creates .cdf, .pdf, and 
; .src files for SPHRAY.  
;




pro make_sphray_source_from_iontable, z, numin, numax, ngamma, my_gammaHI, $
                                      doplot=doplot, dofile=dofile

on_error,2

if n_params() lt 4 then begin
      print
      print,'Usage: '
      print,'make_sphray_source_from_iontable, z, numin, numax, outbase'
      print
      print,' Input:'
      print,'       z      : red shift'
      print,'       numin  : min frequency for spectra [Rydbergs]'
      print,'       numax  : max frequency for spectra [Rydbergs]'
      print
      print,' Output:'
      print,'       creates the files *.cdf, *.pdf, and *.src '
      print,'       ngamma     : number density of background photons'
      print,'       my_gammaHI : photoionization rate from integrating' 
      return
  endif

if n_elements(doplot) eq 0 then doplot=0
if n_elements(dofile) eq 0 then dofile=0

@physical_constants


outbase="hm01qg_z" + string(z, format='(F4.2)')

input = create_struct( 'z', z, 'numin', numin, 'numax', numax, $
                       'lognumin', 0.d0, 'lognumax', 0.0d0, $
                       'outbase', outbase, 'zstr', string(z, format='(F4.2)') )


; get logs of input frequencies
;-------------------------------
input.lognumin = alog10(input.numin)
input.lognumax = alog10(input.numax)


; set constants 
;---------------------
HIth = 13.606
HeIth = 24.587
HeIIth = 54.416

HIth_ry = 1.0
HeIth_ry = HeIth / HIth
HeIIth_ry = HeIIth / HIth

Ryd2erg = 2.17896d-11
HIth_nu = 3.2899d15  ; Hz
SIGMA0 = 6.3d-18





; do input
;-------------------
ion_table_file   = "../../ionization_tables/h1.hdf5"
ion_table_data   = h5_parse(ion_table_file,/read)
ion_table_header = ion_table_data.header
bckgnd_spectrum  = ion_table_header.spectrum
gammahi          = bckgnd_spectrum.gammahi._data
zfile            = bckgnd_spectrum.redshift._data
logflux          = bckgnd_spectrum.logflux._data ; [ergs/s/Hz/cm^2/sr]
logryd           = bckgnd_spectrum.logenergy_ryd._data

fulltable = create_struct( 'GHI', gammahi, 'z', zfile, 'logflux', logflux, $
                           'logryd', logryd )


; find the redshift closest to the redshift we want
;------------------------------------------------------
zdif    = abs( fulltable.z - input.z )
zdifmin = min(zdif, iz)


; use numin and numax to find the indices we want to integrate over
;--------------------------------------------------------------------
nudif    = abs( fulltable.logryd - input.lognumin )
nudifmin = min( nudif, inu )

nudif    = abs( fulltable.logryd - input.lognumax )
nudifmin = min( nudif, fnu ) 
 
if ( fulltable.logryd[inu] gt input.lognumin ) then inu = inu-1
if ( fulltable.logryd[fnu] lt input.lognumax ) then fnu = fnu+1


npts = fnu - inu + 1

table = create_struct( 'GHI', fulltable.GHI[iz], $
                       'z', fulltable.z[iz], $
                       'logryd', fulltable.logryd[inu:fnu], $
                       'logflux', fltarr(npts), $
                       'ryd', fltarr(npts), $
                       'flux', fltarr(npts) )

table.logflux = fulltable.logflux[iz,inu:fnu]
table.ryd = 10^table.logryd
table.flux = 10^table.logflux





; now logryd[inu] and logryd[fnu] bracket the requested frequencies
;--------------------------------------------------------------------

print 
print, "INPUT"
print, "z                   = ", input.z
print, "nu[ryd] min/max     = ", input.numin, input.numax
print, "log nu[ryd] min/max = ", input.lognumin, input.lognumax
print, "out base            = ", outbase
print
print, "CALCULATED"
print, "closest z in file     : ", table.z
print, "log nu bracket values : ", minmax(table.logryd)


; find discontinuities in spectrum
;-------------------------------------
print
print, 'DISCONTINUITIES'
for i = 1, npts-1 do begin
   dx = table.logryd[i] - table.logryd[i-1]
   if dx lt 0.02/10. then begin
      print, 'i,logryd,ryd,dlogryd: ', i, table.logryd[i], table.ryd[i],dx
   endif
endfor
print



; interpolate the flux between the input energies to smooth it slightly.
; below we do an interpolation in points logarithmically spaced
; in energy.
;---------------------------------------------------------------------------
nintrp = npts*2

intrp = create_struct( 'ryd',        dblarr(nintrp),   $
                       'logryd',     dblarr(nintrp),   $
                       'nu',         dblarr(nintrp),   $
                       'lognu',      dblarr(nintrp),   $
                       'flux',       dblarr(nintrp),   $
                       'logflux',    dblarr(nintrp),   $
                       'fluxnu',     dblarr(nintrp),   $
                       'logfluxnu',  dblarr(nintrp),   $
                       'pflux',      dblarr(nintrp),   $
                       'logpflux',   dblarr(nintrp),   $
                       'sigmaHI',    dblarr(nintrp),   $
                       'pdf',        dblarr(nintrp),   $
                       'cdf',        dblarr(nintrp)    ) 


dxintrp         = (input.lognumax - input.lognumin) / (nintrp-1)

intrp.logryd    = input.lognumin + indgen(nintrp) * dxintrp
intrp.ryd       = 10^intrp.logryd

intrp.nu        = intrp.ryd * HIth_nu
intrp.lognu     = alog10( intrp.nu )

intrp.logflux   = interpol( table.logflux, table.logryd, intrp.logryd )
intrp.flux      = 10^intrp.logflux

intrp.fluxnu    = intrp.flux * intrp.ryd * HIth_nu
intrp.logfluxnu = alog10(intrp.fluxnu)

intrp.pflux     = intrp.flux / (intrp.nu * PLANCK)
intrp.logpflux  = alog10( intrp.pflux )

intrp.sigmaHI   = hi_photox_verner( intrp.ryd * 13.6d0 ) 
 


; check if integrated spectrum is equal to table gammaHI
;--------------------------------------------------------------------
xx = intrp.logryd
yy = 4.0d0 * !DPI * alog(10.d0) * intrp.pflux * intrp.nu * intrp.sigmaHI 
my_gammaHI = int_tabulated(xx,yy, /double)

print
print, 'DUMMY CHECK - GAMMAHI'
print, '(log) GHI from table:       ', alog10(table.GHI), table.GHI
print, '(log) GHI from integration: ', alog10(my_gammaHI), my_gammaHI
print, 'table/int:                  ', table.GHI / my_gammaHI
print



; now we integrate the interpolated function to calculate a 
; probability distribution function (PDF) 
;---------------------------------------------------------------
sum = int_tabulated( intrp.nu, intrp.pflux, /double )
intrp.pdf = intrp.pflux/sum
print, 'PDF check: ', int_tabulated( intrp.nu, intrp.pdf, /double )

sum = int_tabulated( intrp.ryd, intrp.pflux * HIth_nu, /double )
intrp.pdf = intrp.pflux * HIth_nu / sum
print, 'PDF check: ', int_tabulated( intrp.ryd, intrp.pdf, /double )




; now we integrate the PDF to calculate a cumulative distribution 
; function (CDF). 
;---------------------------------------------------------------
intrp.cdf[*]= 0.0d0
cdfsum = 0.0d0
for i = 1, nintrp-1 do begin
   a = intrp.ryd[i-1]
   b = intrp.ryd[i]
   ab2 = (a+b)/2.d0
   fy2 = (intrp.pdf[i-1] + intrp.pdf[i])/2.d0
   tmp = (b-a)/6.d0 * ( intrp.pdf[i-1] + 4*fy2 + intrp.pdf[i] )
   cdfsum = cdfsum + tmp
   intrp.cdf[i] = cdfsum
   if intrp.cdf[i] lt intrp.cdf[i-1] then begin
      print, 'cdf inverted'
      stop
   endif
endfor
print, "last entry in cdf array = ", intrp.cdf[nintrp-1]
intrp.cdf = intrp.cdf / intrp.cdf[nintrp-1]






; plot interpolated values over table values
;-------------------------------------------------------------------
if doplot then begin
   altay_set_x, xsize=1200, ysize=400
   !P.multi=[0,3,1]

   ytitle= TexToIdl("Log J \nu [erg cm^{-2} s^{-1} sr^{-1}]")
   xtitle= TexToIdl("Log Energy [Rydbergs]")
   ybigrange = [-50,50]
   HeIplt  = alog10(HeIth_ry)
   HeIIplt = alog10(HeIIth_ry)

   altay_plot, table.logryd, alog10(table.flux * table.ryd * HIth_nu), $
               xtitle=xtitle, ytitle=ytitle, $
               thick=8, charsize=3, charthick=2, $
               yrange=[-10,-5]   
   altay_oplot, intrp.logryd, intrp.logfluxnu, color=250   
   altay_oplot, [HeIplt,  HeIplt],  ybigrange, color=150, thick=1
   altay_oplot, [HeIIplt, HeIIplt], ybigrange, color=150, thick=1
   
   
   ytitle= TexToIdl("Log J [erg cm^{-2} Hz^{-1} s^{-1} sr^{-1}]")
   altay_plot, table.logryd, table.logflux, $
               xtitle=xtitle, ytitle=ytitle, $
               thick=8, charsize=3, charthick=2, $
               yrange=[-26,-20.5]
   altay_oplot, intrp.logryd, intrp.logflux, color=250
   altay_oplot, [HeIplt,  HeIplt],  ybigrange, color=150, thick=1
   altay_oplot, [HeIIplt, HeIIplt], ybigrange, color=150, thick=1
   
   
   xyouts, 0.3, -21.0, "z = " + input.zstr, $
           color=0, charthick=2, charsize=2


altay_plot, intrp.logryd, intrp.pdf, thick=2, $
            xtitle="Log Energy [Rydbergs]", $
            ytitle=TexToIdl("PDF and CDF"), $
            charsize=3, charthick=2, color=50, $
            yrange=[-0.05,1.7]

altay_oplot, intrp.logryd, intrp.cdf, thick=2, color=250


!P.multi=0


endif



; now we integrate to get the number density of photons
; note n = (photons cm-2  s-1) / c
;-------------------------------------------------------
xx = intrp.ryd
yy = (4 * !DPI / LIGHT) * intrp.pflux * HIth_nu
ngamma = int_tabulated(xx, yy, /double) 
plane_flux = ngamma * LIGHT

print
print, " n [photons / cm^3]     = ", ngamma
print, " F [photons / cm^2 / s] = ", plane_flux
print


; calculate gammaHI from number density.  this is the 
; fairest number to compare to the code output.
;--------------------------------------------------------------------
xx = intrp.ryd
yy = LIGHT * ngamma * intrp.pdf * intrp.sigmaHI
my_gammaHI = int_tabulated(xx,yy, /double)


print
print, 'GAMMAHI from ngamma - Compare to SPHRAY output'
print, 'GHI from ngamma: ', alog10(my_gammaHI)
print



; write out to .cdf, .pdf, and .src files
;-----------------------------------------
if dofile then begin
   cdffile='../hm01/'+outbase+".cdf"
   spcfile='../hm01/'+outbase+".pdf"
   srcfile='../../source_templates/'+outbase+".src"

   openw, luncdf, cdffile, /get_lun
   openw, lunspc, spcfile, /get_lun
   openw, lunsrc, srcfile, /get_lun

   printf, luncdf, 1
   printf, luncdf, nintrp, numin, numax

   for i = 0,nintrp-1 do begin
      printf, luncdf, intrp.ryd[i], intrp.cdf[i], format='(2E17.7)'
      printf, lunspc, intrp.ryd[i], intrp.pdf[i], format='(2E17.7)'
   endfor

   free_lun, luncdf
   free_lun, lunspc


   printf, lunsrc, '1       ! number of sources in this file'
   printf, lunsrc, '1       ! number of sources in this source snapshot'
   printf, lunsrc, '1       ! number of files in the source snapshot'
   printf, lunsrc, '1.0d7   ! number of rays to trace [depends on config file]'
   printf, lunsrc, '1.0d50  ! luminosity unit in photons/s'
   printf, lunsrc, ''
   printf, lunsrc, '! pos     | vel    | lum    |   spec    |  emis'
   printf, lunsrc, '0.0 0.0 0.0    0.0 0.0 0.0  ', ngamma, '    1    -3'
   free_lun, lunsrc
endif




end
