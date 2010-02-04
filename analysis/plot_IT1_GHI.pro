; SPHRAY readin and plot for IT1 photoionization rates
;========================================================

;====================================================================
; set plot to screen/file option and input files
;====================================================================
ps=0       ; ps=0 directs output to screen, ps=1 directs output to psfile
makepng=1  ; if ps=0 and makepng=1 then tries a screen capture to png

; SPHRAY file IO
;----------------
snapdir  = "../../sphray_output/IT1"
snapbase = "iliev_test1"
snapnum = [1,3,5]
snapnumstr = string(snapnum, format="(I3.3)")


sfile  = snapdir + "/" + snapbase + "_" + snapnumstr
psfile = "T1gHI.eps"
pngfile = "T1gHI.png"


print
print, "sphray output file:", sfile
print, "eps file if doing post script output:", psfile
print, "png file if doing png output:", pngfile
print



; define some physical constants
;----------------------------------

PROTONMASS = 1.6726d-24 ; proton mass [g]
KPC2CM = 3.08568025d21  ; cm per kpc
SIGMA0 = 6.3d-18        ; HI photo-absorption cross-section

GLEN   = 3.085678d21   ; [cm h^-1]
GMASS  = 1.989d43      ; [g h^-1]
GVEL   = 1.0d5         ; [cm s^-1]

GTIME  = GLEN / GVEL               ; [s h^-1]
GRHO   = GMASS / GLEN^3            ; [g cm^-3 h^2]
GPRS   = GMASS / GLEN / GTIME^2    ; [g cm^-1 s^-2 h^2]
GENRG  = GMASS * GLEN^2 / GTIME^2  ; [erg h^-1]
GLUM   = GENRG / GTIME             ; [erg/s]



;====================================================================
; read sphray snapshots
;====================================================================
gHI = create_struct( "name", "Hydrogen Photoionization Rates" )
nHI = create_struct( "name", "Neutral Hydrogen Number Density" )
for isnap = 0, n_elements(snapnum)-1 do begin
    data = read_sphray(sfile[isnap])
    ngas = data.heads[0].npar_snap[0]

    gHI = create_struct( gHI, $
          "_" + string( snapnum[isnap], format="(I3.3)"), fltarr(ngas) )
    gHI.(isnap+1) = data.GHI

    nH = data.rho * GRHO / PROTONMASS
    nHI = create_struct( nHI, $
          "_" + string( snapnum[isnap], format="(I3.3)"), dblarr(ngas) )
    nHI.(isnap+1) = data.xHI * nH
endfor







;====================================================================
; sort the particles into radial bins. 
;====================================================================
data.pos = data.pos - 6.6
cen = [0.0,0.0,0.0]
r = sqrt( (data.pos[0,*]-cen[0])^2 + $
          (data.pos[1,*]-cen[1])^2 + $
          (data.pos[2,*]-cen[2])^2   )

; calculate analytic gammaHI assuming equilibrium
;---------------------------------------------------

nbins=60L
mlocs = fltarr(nbins)
minr=0.0
maxr=1.0
bs=(maxr-minr)/nbins

bdata = create_struct( "locs",     dblarr(nbins), $
                       "rin",      dblarr(nbins), $
                       "rcen",     dblarr(nbins), $
                       "rout",     dblarr(nbins), $
                       "dl",       dblarr(nbins), $
                       "flux",     dblarr(nbins), $
                       "nHI_001",  dblarr(nbins), $
                       "nHI_003",  dblarr(nbins), $
                       "nHI_005",  dblarr(nbins), $
                       "tau_001",  dblarr(nbins), $
                       "tau_003",  dblarr(nbins), $
                       "tau_005",  dblarr(nbins), $
                       "GHIs_001", dblarr(nbins), $
                       "GHIs_003", dblarr(nbins), $
                       "GHIs_005", dblarr(nbins), $
                       "GHIa_001", dblarr(nbins), $
                       "GHIa_003", dblarr(nbins), $
                       "GHIa_005", dblarr(nbins)  )                       

print, "SPHRAY radial bins = ", nbins

hist=histogram(r/6.6,min=0.0,max=1.0-bs,nbins=nbins, $
                 reverse_indices=ri)

bdata.tau_001 = 0.0d0
bdata.tau_003 = 0.0d0
bdata.tau_005 = 0.0d0

for i = 0,nbins-1 do begin

    bdata.locs[i] = (i + 0.5) * bs
    bdata.rin[i]  = (i      ) * bs * 6.6d0 * KPC2CM
    bdata.rcen[i] = (i + 0.5) * bs * 6.6d0 * KPC2CM
    bdata.rout[i] = (i + 1  ) * bs * 6.6d0 * KPC2CM

    bdata.dl[i]   = bdata.rout[i] - bdata.rin[i]
    
    indx = ri[ ri[i] : ri[i+1]-1 ]

    ; photoionization rates from SPHRAY
    bdata.GHIs_001[i] = mean(gHI.(1)[indx], /double)
    bdata.GHIs_003[i] = mean(gHI.(2)[indx], /double)
    bdata.GHIs_005[i] = mean(gHI.(3)[indx], /double)

    ; analytic photoionization rates using the SPHRAY xHI
    bdata.flux[i] = 5.0d48 / (4.0d0 * !DPI * bdata.rcen[i]^2)    

    bdata.nHI_001[i] = mean( nHI.(1)[indx], /double )
    bdata.nHI_003[i] = mean( nHI.(2)[indx], /double )
    bdata.nHI_005[i] = mean( nHI.(3)[indx], /double )
 
    bdata.tau_001[i] = bdata.dl[i] * bdata.nHI_001[i] * SIGMA0
    bdata.tau_003[i] = bdata.dl[i] * bdata.nHI_003[i] * SIGMA0
    bdata.tau_005[i] = bdata.dl[i] * bdata.nHI_005[i] * SIGMA0

    bdata.GHIa_001[i] = bdata.flux[i] * SIGMA0 * exp( -total(bdata.tau_001[0:i]) )
    bdata.GHIa_003[i] = bdata.flux[i] * SIGMA0 * exp( -total(bdata.tau_003[0:i]) )
    bdata.GHIa_005[i] = bdata.flux[i] * SIGMA0 * exp( -total(bdata.tau_005[0:i]) )

endfor

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

sphrayline=2
sphraycolor=254

anacolor=50
analine=0
 
      
loadct, 39
plot,      [0], [0],  $
           xstyle=1, xrange=[0.0,0.99], xtitle="r/L!lbox", $
           ystyle=1, yrange=[-18,-10], ytitle=TexToIdl("Log \Gamma_{HI}"), $
           position=[0.18,0.15,0.95,0.95], /nodata, color=0, $
           background=255, charsize=charsize, charthick=charthick, $
           xthick=mythick, ythick=mythick

oplot,bdata.locs, alog10(bdata.GHIs_001), linestyle=sphrayline, color=sphraycolor, thick=mythick
oplot,bdata.locs, alog10(bdata.GHIs_003), linestyle=sphrayline, color=sphraycolor, thick=mythick
oplot,bdata.locs, alog10(bdata.GHIs_005), linestyle=sphrayline, color=sphraycolor, thick=mythick
 
oplot,bdata.locs, alog10(bdata.GHIa_001), linestyle=analine, color=anacolor, thick=mythick
oplot,bdata.locs, alog10(bdata.GHIa_003), linestyle=analine, color=anacolor, thick=mythick
oplot,bdata.locs, alog10(bdata.GHIa_005), linestyle=analine, color=anacolor, thick=mythick


legend, ["Analytic", "SPHRAY"], $
  linestyle=[0,2], color=[50,254], textcolor=[0,0], $
  charsize=2, charthick=2, thick=2, position=[0.4, -10.5], box=0



if ps then begin
   device, /close
endif else begin
   if makepng then screen_to_png, pngfile
endelse
set_plot, "x"
  

end
