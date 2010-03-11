files=["../rate_comp_D1n5.txt","../rate_comp_D1n1.txt"]

logT = fltarr(100,2)
x_h  = fltarr(100,2)
x_c  = fltarr(100,2)
x_i  = fltarr(100,2)
x_hB = fltarr(100,2)

for i = 0, n_elements(files)-1 do begin
    readcol, files[i], a, b, c, d, e
    logT[*,i] = a
    x_h[*,i] = b
    x_c[*,i] = c
    x_i[*,i] = d
    x_hB[*,i] = e
endfor


altay_set_x

altay_plot, [0], [0], xrange=[0.9,6.1], yrange=[-8.0,0.0], $
  xtitle="Log T", ytitle=TexToIdl("Equilibrium Log x_{HI}"), charsize=2

;altay_oplot, logT[*,0], alog10(x_c[*,0]), color=50
altay_oplot, logT[*,0], alog10(x_h[*,0]), color=50
altay_oplot, logT[*,0], alog10(x_hB[*,0]), color=100
altay_oplot, logT[*,0], alog10(x_i[*,0]), color=254



;altay_oplot, logT[*,1], alog10(x_c[*,1]), color=50
altay_oplot, logT[*,1], alog10(x_h[*,1]), color=50
altay_oplot, logT[*,1], alog10(x_hB[*,1]), color=100
altay_oplot, logT[*,1], alog10(x_i[*,1]), color=254


xyouts, 3.5, -7., "nH = 1.0e-5", color=1, charsize=2, charthick=2
xyouts, 4.0, -1., "nH = 1.0e-1", color=1, charsize=2, charthick=2

legend, ["Hui97-A", "Hui97-B", "Cloudy"], $
  linestyle=[0,0,0], color=[50,100,254], textcolor=[0,0,0], $
  charsize=2, charthick=2, thick=2, position=[1.2, -2.], box=0

screen_to_png, "rate_compare.png"

end
