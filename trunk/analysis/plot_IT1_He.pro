; SPHRAY readin and plot for IT1
;========================================

;====================================================================
; set plot to screen/file option and input files
;====================================================================
ps=0       ; ps=0 directs output to screen, ps=1 directs output to psfile
makepng=0  ; if ps=0 and makepng=1 then tries a screen capture to png

; SPHRAY file IO
;----------------
snapdir  = "../../sphray_output/IT1_He"
snapbase = "snap"
snapnum = 1
snapnumstr = string(snapnum, format="(I3.3)")


sfile  = snapdir + "/" + snapbase + "_" + snapnumstr
psfile = "T1x_He" + snapnumstr + ".eps"
pngfile = "T1x_He" + snapnumstr + ".png"

print
print, "sphray output file:", sfile
print, "eps file if doing post script output:", psfile
print



; define some physical constants
;----------------------------------

PROTONMASS = 1.6726d-24 ; proton mass [g]
GLEN  = 3.085678d21   ; [cm h^-1]
GMASS = 1.989d43      ; [g h^-1]
GVEL  = 1.0d5         ; [cm s^-1]

GTIME = GLEN / GVEL               ; [s h^-1]
GRHO  = GMASS / GLEN^3            ; [g cm^-3 h^2]
GPRS  = GMASS / GLEN / GTIME^2    ; [g cm^-1 s^-2 h^2]
GENRG = GMASS * GLEN^2 / GTIME^2  ; [erg h^-1]
GLUM  = GENRG / GTIME             ; [erg/s]



;====================================================================
; read sphray snapshot
;====================================================================
data = read_sphray(sfile) 


;====================================================================
; sort the particles into radial bins. 
;====================================================================
data.pos = data.pos - 6.6
cen = [0.0,0.0,0.0]
r = sqrt( (data.pos[0,*]-cen[0])^2 + $
          (data.pos[1,*]-cen[1])^2 + $
          (data.pos[2,*]-cen[2])^2   )
  
bins=50L
mlocs = fltarr(bins)
minr=0.0
maxr=1.0
bs=(maxr-minr)/bins
 
print, "SPHRAY radial bins = ", bins

hist=histogram(r/6.6,min=0.0,max=1.0-bs,nbins=bins, $
                 reverse_indices=ri)

y1 = mlocs
y2 = mlocs

w1 = mlocs
w2 = mlocs
w3 = mlocs



for i = 0,bins-1 do begin
    mlocs[i] = (i + 0.5) * bs
    indx = ri[ ri[i] : ri[i+1]-1 ]

    y1[i] = mean(data.xHI[indx], /double)
    y2[i] = mean(1.0d0 - data.xHI[indx], /double)

    w1[i] = mean(data.xHeI[indx], /double)
    w2[i] = mean(data.xHeII[indx], /double)
    w3[i] = 1.0d0 - w1[i] - w2[i]
endfor



;====================================================================
; create plot
;====================================================================
if (ps EQ 1) then begin
   altay_set_ps, psfile
    mythick=7.0
    charsize = 2.0   ; characters
    charthick = 4.0
endif else begin
   altay_set_x
    mythick=3.0
    charsize = 2.0   ; characters
    charthick = 2.0
endelse


 

Hcolor=1     ; black
HeIcolor=254
HeIIcolor=150
HeIIIcolor=70

Hline=2   
Heline=0

               
loadct, 39
plot,      mlocs, y1,  $
           xstyle=1, xrange=[0.0,0.99], xtitle="r/L!lbox", $
           ystyle=1, yrange=[1.0E-5,1.1], /ylog, $
           ytitle="H and He Ionization Fractions", $
           position=[0.18,0.15,0.95,0.95], /nodata, color=0, $
           background=255, charsize=charsize, charthick=charthick, $
           xthick=mythick, ythick=mythick


oplot,mlocs,y1, linestyle=Hline, color=sphraycolor, thick=mythick
oplot,mlocs,y2, linestyle=Hline, color=sphraycolor, thick=mythick

oplot,mlocs,w1, linestyle=Heline, color=HeIcolor, thick=mythick
oplot,mlocs,w2, linestyle=Heline, color=HeIIcolor, thick=mythick
oplot,mlocs,w3, linestyle=Heline, color=HeIIIcolor, thick=mythick

legend, ["HI,HII", "HeI", "HeII", "HeIII"], $
        linestyle=[Hline,Heline,Heline,Heline], $
        color=[Hcolor,HeIcolor,HeIIcolor,HeIIIcolor], textcolor=0, $
        charsize=2, charthick=2, thick=2, position=[0.1, 2.0e-4], $
        box=1, outline_color=0, number=0.01, /clear


if ps then begin
   device, /close
endif else begin
   if makepng then screen_to_png, pngfile
endelse
set_plot, "x"
  

end
