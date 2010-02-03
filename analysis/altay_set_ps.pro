pro altay_set_ps, psfile, xsize=xsize, ysize=ysize

on_error, 2

if n_elements(psfile) eq 0 or n_params() ne 1 then begin
    print
    print,'Usage: altay_set_ps, psfile'
    print
    print,' Input:'
    print,'       psfile  : name of postscript file for output'
    print
    print,' Keywords:'
    print,'       xsize   : x size in cm (default=16.0)'
    print,'       ysize   : y size in cm (default=16.0)'
    print
    print,' Output:'
    print
    return
endif


if n_elements(xsize) eq 0 then xsize=16.
if n_elements(ysize) eq 0 then ysize=16.

set_plot, 'PS', /copy, /interpolate
device, encapsulated=1, filename=psfile, $
  xsize=xsize, ysize=ysize, /color, bits_per_pixel=8

return

end
