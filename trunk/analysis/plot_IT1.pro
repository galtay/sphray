; SPHRAY readin and plot for IT1
;========================================

;====================================================================
; set plot to screen/file option and input files
;====================================================================
ps=0       ; ps=0 directs output to screen, ps=1 directs output to psfile
makepng=0  ; if ps=0 and makepng=1 then tries a screen capture to png

; SPHRAY file IO
;----------------
snapdir  = "../../sphray_output/IT1/r6"
snapbase = "snap"
snapnum = 1
snapnumstr = string(snapnum, format="(I3.3)")


sfile  = snapdir + "/" + snapbase + "_" + snapnumstr
psfile = "T1x" + snapnumstr + ".eps"
pngfile = "T1x" + snapnumstr + ".png"
cmpfile= "CmpData/CmpT1_" + snapnumstr + "x.txt"

print
print, "comparison project file:", cmpfile
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
  
bins=60L
mlocs = fltarr(bins)
minr=0.0
maxr=1.0
bs=(maxr-minr)/bins
 
print, "SPHRAY radial bins = ", bins

hist=histogram(r/6.6,min=0.0,max=1.0-bs,nbins=bins, $
                 reverse_indices=ri)

y1 = mlocs
y2 = mlocs

for i = 0,bins-1 do begin
    mlocs[i] = (i + 0.5) * bs
    indx = ri[ ri[i] : ri[i+1]-1 ]
    y1[i] = 1 - mean(data.xHI[indx], /double)
    y2[i] = mean(data.xHI[indx], /double)
endfor

;====================================================================
; define arrays for comparison project data
;====================================================================
Ncmpbins=100
locs = fltarr(Ncmpbins)

xc2ray1   = fltarr(Ncmpbins)
xotvet1   = fltarr(Ncmpbins)
xcrash1   = fltarr(Ncmpbins)
xrsph1    = fltarr(Ncmpbins)
xart1     = fltarr(Ncmpbins)
xftte1    = fltarr(Ncmpbins)
xsimplex1 = fltarr(Ncmpbins)
xzeus1    = fltarr(Ncmpbins)
xflash1   = fltarr(Ncmpbins)
xift1     = fltarr(Ncmpbins)

;====================================================================
; read in comparison project data
;====================================================================
openr, lun, cmpfile, /get_lun
for i = 0,Ncmpbins-1 do begin
    readf, lun, xx
    locs[i] = xx
endfor  

; c2ray data
for i = 0,Ncmpbins-1 do begin
    readf, lun, xx
    xc2ray1[i]=xx
endfor 

; otvet data
for i = 0,Ncmpbins-1 do begin
    readf, lun, xx
    xotvet1[i]=xx
endfor
  
; crash data
for i = 0,Ncmpbins-1 do begin
    readf, lun, xx
    xcrash1[i]=xx
endfor  
  
; rsph data
for i = 0,Ncmpbins-1 do begin
    readf,lun,xx
    xrsph1[i]=xx
endfor  

; art data
for i = 0,Ncmpbins-1 do begin
    readf, lun, xx
    xart1[i]=xx
endfor  
  
; ftte data
for i = 0,Ncmpbins-1 do begin
    readf,lun,xx
    xftte1[i]=xx
endfor  

; simplex data
for i = 0,Ncmpbins-1 do begin
    readf,lun,xx
    xsimplex1[i]=xx
endfor  

; zeus data
for i = 0,Ncmpbins-1 do begin
    readf,lun,xx
    xzeus1[i]=xx
endfor  

; flash data
for i = 0,Ncmpbins-1 do begin
    readf,lun,xx
    xflash1[i]=xx
endfor  

; ift data
for i = 0,Ncmpbins-1 do begin
    readf,lun,xx
    xift1[i]=xx
endfor   

free_lun, lun

;====================================================================
; create plot
;====================================================================
if (ps EQ 1) then begin
   altay_set_ps, psfile
    cmpthick=5.0     ; line thickness
    mythick=7.0
    charsize = 2.0   ; characters
    charthick = 4.0
endif else begin
   altay_set_x
    cmpthick=2.0     ; line thickness
    mythick=3.0
    charsize = 2.0   ; characters
    charthick = 2.0
endelse


 
c2color=254       ; red     ct 39
otvetcolor=50     ; blue    ct 39
crashcolor=150    ; green   ct 39
rsphcolor=1       ; black   ct 39
artcolor=120      ; teal    ct 39
fttecolor=130     ; pink    ct 12
simplexcolor=150  ; brown ct 24 
zeuscolor=100     ; other green ct 24
flashcolor=206    ; flesh ct 28
iftcolor=95       ; blue/green ct 39
coralcolor=254    ; red ct 39
sphraycolor=1     ; black

c2line=0       ; solid line
otvetline=5    ; long dashes
crashline=2    ; short dashes
rsphline=1     ; dots
artline=3      ; dash dot
ftteline=0     ; solid line
simplexline=3  ; dashed dot
zeusline=1     ; dots
flashline=2    ; short dashes
iftline=3      ; dash dot
sphrayline=0   ; solid line

               
loadct, 39
plot,      mlocs, y1,  $
           xstyle=1, xrange=[0.0,0.99], xtitle="r/L!lbox", $
           ystyle=1, yrange=[1.0E-5,1.1], /ylog, ytitle="x, 1-x", $
           position=[0.18,0.15,0.95,0.95], /nodata, color=0, $
           background=255, charsize=charsize, charthick=charthick, $
           xthick=mythick, ythick=mythick

oplot,locs,xc2ray1, linestyle=c2line, color=c2color, thick=cmpthick
oplot,locs,1-xc2ray1, linestyle=c2line, color=c2color, thick=cmpthick
oplot,locs,xotvet1, linestyle=otvetline, color=otvetcolor, thick=cmpthick
oplot,locs,1-xotvet1, linestyle=otvetline, color=otvetcolor, thick=cmpthick

oplot,locs,xcrash1,linestyle=crashline, color=crashcolor, thick=cmpthick
oplot,locs,1-xcrash1,linestyle=crashline, color=crashcolor, thick=cmpthick
oplot,locs,xrsph1,linestyle=rsphline, color=rsphcolor, thick=cmpthick
oplot,locs,1-xrsph1,linestyle=rsphline, color=rsphcolor, thick=cmpthick

oplot,locs,xart1,linestyle=artline, color=artcolor, thick=cmpthick
oplot,locs,1-xart1,linestyle=artline, color=artcolor, thick=cmpthick

loadct, 12
oplot,locs,xftte1,linestyle=ftteline, color=fttecolor, thick=cmpthick
oplot,locs,1-xftte1,linestyle=ftteline, color=fttecolor, thick=cmpthick

loadct, 24
oplot,locs,xsimplex1,linestyle=simplexline, color=simplexcolor, thick=cmpthick
oplot,locs,1-xsimplex1,linestyle=simplexline, color=simplexcolor, thick=cmpthick
oplot,locs,xzeus1,linestyle=zeusline, color=zeuscolor, thick=cmpthick
oplot,locs,1-xzeus1,linestyle=zeusline, color=zeuscolor, thick=cmpthick

loadct, 28
oplot,locs,xflash1,linestyle=flashline, color=flashcolor, thick=cmpthick
oplot,locs,1-xflash1,linestyle=flashline, color=flashcolor, thick=cmpthick

loadct, 39
oplot,locs,xift1,linestyle=iftline, color=iftcolor, thick=cmpthick
oplot,locs,1-xift1,linestyle=iftline, color=iftcolor, thick=cmpthick

oplot,mlocs,y1, linestyle=sphrayline, color=sphraycolor, thick=mythick
oplot,mlocs,y2, linestyle=sphrayline, color=sphraycolor, thick=mythick


if ps then begin
   device, /close
endif else begin
   if makepng then screen_to_png, pngfile
endelse
set_plot, "x"
  

end
