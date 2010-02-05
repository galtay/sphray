pro altay_plot, x, y, _extra=extra

on_error,2

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



plot, x, y, $
  /xstyle, xrange=[xmin,xmax], xthick=xthick, $
  /ystyle, yrange=[ymin,ymax], ythick=ythick, $
  charsize=charsize, charthick=charthick, $
  background=background, color=color, _strict_extra=extra



end
