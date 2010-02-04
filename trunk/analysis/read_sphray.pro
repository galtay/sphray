function read_sphray, snapbase



gheads = gadget_header_read(snapbase, /flag_sphray)


ngas   = gheads[0].npar_snap[0]
nfiles = gheads[0].nfiles


data = create_struct( "heads", gheads,         $
                      "pos",   fltarr(3,ngas), $
                      "vel",   fltarr(3,ngas), $
                      "id",    lonarr(ngas),   $
                      "mass",  fltarr(ngas),   $
                      "u",     fltarr(ngas),   $
                      "rho",   fltarr(ngas),   $
                      "ye",    fltarr(ngas),   $
                      "xHI",   fltarr(ngas),   $
                      "hsml",  fltarr(ngas),   $
                      "T",     fltarr(ngas)    )

if gheads[0].flag_Hmf then data = create_struct( data, $
                                                 "Hmf", fltarr(ngas) )

if gheads[0].flag_Hemf then data = create_struct( data, $
                                                  "Hemf", fltarr(ngas) )

if gheads[0].flag_helium then data = create_struct( data,                 $
                                                    "xHeI", fltarr(ngas), $
                                                    "xHeII", fltarr(ngas) )

if gheads[0].flag_gammaHI then data = create_struct( data, $
                                                     "GHI", fltarr(ngas) )

if gheads[0].flag_cloudy then data = create_struct( data, $
                                                    "xHI_cloudy", fltarr(ngas) )

if gheads[0].flag_eos then data = create_struct( data, $
                                                 "eos", fltarr(ngas) )



dumhead = gadget_header_define( /flag_sphray )

nred=0LL
for ifile = 0, nfiles-1 do begin

   ngas1 = gheads[ifile].npar_file[0]
   if ngas1 eq 0 then continue

   if nfiles eq 1 then begin
      f = snapbase
   endif else begin
      f = snapbase + "." + strcompress(ifile, /remove_all)
   endelse

   openr, lun, f, /f77_unformatted, /get_lun
   readu, lun, dumhead


   ; positions
   tmp = fltarr(3,ngas1)
   readu, lun, tmp
   data.pos[*,nred:nred+ngas1-1] = tmp

   ;velocities
   tmp = fltarr(3,ngas1)
   readu, lun, tmp
   data.vel[*,nred:nred+ngas1-1] = tmp

   ; IDs
   tmp = lonarr(ngas1)
   readu, lun, tmp
   data.id[nred:nred+ngas1-1] = tmp

   ; mass
   tmp = fltarr(ngas1)
   readu, lun, tmp
   data.mass[nred:nred+ngas1-1] = tmp

   ; u
   tmp = fltarr(ngas1)
   readu, lun, tmp
   data.u[nred:nred+ngas1-1] = tmp

   ; rho
   tmp = fltarr(ngas1)
   readu, lun, tmp
   data.rho[nred:nred+ngas1-1] = tmp

   ; ye
   tmp = fltarr(ngas1)
   readu, lun, tmp
   data.ye[nred:nred+ngas1-1] = tmp

   ; xHI
   tmp = fltarr(ngas1)
   readu, lun, tmp
   data.xHI[nred:nred+ngas1-1] = tmp

   ; hsml
   tmp = fltarr(ngas1)
   readu, lun, tmp
   data.hsml[nred:nred+ngas1-1] = tmp
   
   ; T
   tmp = fltarr(ngas1)
   readu, lun, tmp
   data.T[nred:nred+ngas1-1] = tmp


   if gheads[0].flag_Hmf then begin
       tmp = fltarr(ngas1) 
       readu, lun, tmp
       data.Hmf[nred:nred+ngas1-1] = tmp
   endif

   if gheads[0].flag_Hemf then begin
       tmp = fltarr(ngas1) 
       readu, lun, tmp
       data.Hemf[nred:nred+ngas1-1] = tmp
   endif

   if gheads[0].flag_helium then begin
       tmp = fltarr(ngas1) 
       readu, lun, tmp
       data.xHeI[nred:nred+ngas1-1] = tmp
       readu, lun, tmp
       data.xHeII[nred:nred+ngas1-1] = tmp
   endif

   if gheads[0].flag_gammaHI then begin
       tmp = fltarr(ngas1) 
       readu, lun, tmp
       data.GHI[nred:nred+ngas1-1] = tmp
   endif

   if gheads[0].flag_cloudy then begin
       tmp = fltarr(ngas1) 
       readu, lun, tmp
       data.xHI_cloudy[nred:nred+ngas1-1] = tmp
   endif

   if gheads[0].flag_eos then begin
       tmp = fltarr(ngas1) 
       readu, lun, tmp
       data.eos[nred:nred+ngas1-1] = tmp
   endif

   nred = nred + ngas1
   tmp=0
   free_lun, lun

endfor



gheads=0





return, data

end
