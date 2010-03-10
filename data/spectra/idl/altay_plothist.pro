pro altay_plothist, x, ct=ct, axiscolor=axiscolor, _extra=extra

on_error,2


if n_elements(ct) eq 0 then ct=39
if n_elements(axiscolor) eq 0 then axiscolor=0


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


plothist, x, xhist, yhist, /noplot, /overplot, _strict_extra=extra


; plot axes
;------------
xmin=min(xhist)
xmax=max(xhist)
ymin=min(yhist)
ymax=max(yhist)

dx = xmax-xmin
xmin = xmin - dx * 0.05
xmax = xmax + dx * 0.05


; account for /ylog choice in range selection
;-----------------------------------------------
if n_elements(extra) ne 0 then begin

    tags = tag_names(extra)
    tagmatch = where( tags eq "YLOG", ylogcount )
    if ylogcount ne 0 then begin

        dy = alog10(ymax) - alog10(ymin)
        ymin = alog10(ymin) - dy * 0.05
        ymax = alog10(ymax) + dy * 0.05

        ymin = 10^ymin
        ymax = 10^ymax

    endif else begin

        dy = ymax-ymin
        ymin = ymin - dy * 0.05
        ymax = ymax + dy * 0.05
        
    endelse

endif else begin

    dy = ymax-ymin
    ymin = ymin - dy * 0.05
    ymax = ymax + dy * 0.05

endelse




loadct, ct
plothist, x, $
  /xstyle, xrange=[xmin,xmax], xthick=xthick, $
  /ystyle, yrange=[ymin,ymax], ythick=ythick, $
  thick=thick, charsize=charsize, charthick=charthick, $
  background=background, color=color, axiscolor=axiscolor, $
  _strict_extra=extra
 






end
