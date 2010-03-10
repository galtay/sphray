pro altay_set_x, xsize=xsize, ysize=ysize, windownum=windownum

on_error, 2

if n_params() ne 0 then begin
    print
    print,'Usage: altay_set_x'
    print
    print,' Input:'
    print
    print,' Keywords:'
    print,'       xsize   : x size in pixels (default=700)'
    print,'       ysize   : y size in pixels (default=700)'
    print
    print,' Output:'
    print
    return
endif


if n_elements(xsize) eq 0 then xsize=700
if n_elements(ysize) eq 0 then ysize=700
if n_elements(windownum) eq 0 then windownum=0

set_plot, 'X', /copy
device, decomposed=0, /true_color, retain=2 
window, windownum, xsize=xsize, ysize=ysize

return

end
