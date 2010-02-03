pro altay_plot, x, y, ct=ct, axiscolor=axiscolor, axisct=axisct, _extra=extra

on_error,2


if n_elements(ct) eq 0 then ct=39
if n_elements(axiscolor) eq 0 then axiscolor=0
if n_elements(axisct) eq 0 then axisct=39

background=255
color=0

if (!D.name eq "PS") then begin

    charsize=1.5
    charthick=6.0
    xthick=6.0
    ythick=6.0
    thick=6.0

endif else if (!D.name eq "X") then begin

    charsize=2.0
    charthick=2.0
    xthick=2.0
    ythick=2.0
    thick=2.0

endif else begin

    message, "  Plotting device not recognized"

endelse



; dont want extra.color to overide axiscolor
;--------------------------------------------
if n_elements(extra) ne 0 then begin
    tags = tag_names(extra)
    tagmatch = where( tags eq "COLOR", count )
    if count ne 0 then begin
        origcolor = extra.color
        extra.color=axiscolor
    endif
endif


; plot axes
;------------
xmin=min(x)
xmax=max(x)
ymin=min(y)
ymax=max(y)

dx = xmax-xmin
dy = ymax-ymin

xmin = xmin - dx * 0.05
xmax = xmax + dx * 0.05
ymin = ymin - dy * 0.05
ymax = ymax + dy * 0.05


loadct, axisct
plot, x, y, $
  /xstyle, xrange=[xmin,xmax], xthick=xthick, $
  /ystyle, yrange=[ymin,ymax], ythick=ythick, $
  charsize=charsize, charthick=charthick, $
  background=background, color=axiscolor, /nodata, _strict_extra=extra


; if we need to, switch extra.color back for plotting the data
;--------------------------------------------------------------
if n_elements(extra) ne 0 then begin
    if count ne 0 then extra.color = origcolor
endif

; plot data
;-----------
loadct, ct
oplot, x, y, color=color, thick=thick, _extra=extra



end
