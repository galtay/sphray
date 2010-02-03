pro gadget_header_print, ghead

;+
; NAME:
;    gadget_header_print
;
; PURPOSE:
;    Prints a gadget header to screen. 
;
; CATEGORY:
;    IDL Hands Library - Gadget Routines
;
; CALLING SEQUENCE:
;    gadget_header_print, ghead
;
; INPUTS:
;    ghead - Gadget header structure to print. 
;
; OUTPUTS:
;    Header information to screen.
;
; RESTRICTIONS:
;    Header structure must be scalar.
;
; EXAMPLE:
;    snapbase = "/home/galtay/Data/DMG_N64_L1/snapshot_002"
;    gheads = gadget_header_read(snapbase, /silent)
;    gadget_header_print, gheads[0]
;
; MODIFICATION HISTORY:
;    Version 1, Gabriel Altay, CMU, 14 Oct 2009.
;-


; check arguments
;-----------------

on_error, 2 ; return to calling routine on error

if (n_elements(ghead) eq 0) then begin
    message, "usage: gadget_header_print, <gadget_header_struct>"
endif

if (n_elements(ghead) ne 1) then begin
    message, "<gadget_header_struct> must be scalar"
endif




; formatted print the header
;-----------------------------

print, " gadget header ", format='( "* ", 30("-"), A, 30("-"), " *" )' 
print, format='( "* ", T78, " *" )'
  
print, "gas", "halo", "disk", "bulge", "star", "bndry", $
  format='( "* ", T9,A, T21,A, T33,A, T45,A, T57,A, T69,A, T78, " *" )'  

str = strtrim( string(ghead.npar_file), 2)
print, "file:", str, $
  format='( "* ", A,T9,A, T21,A, T33,A, T45,A, T57,A, T69,A, T78, " *" )'  
 
str = strtrim( string(ghead.npar_snap), 2)
print, "snap:", str, $
  format='( "* ", A,T9,A, T21,A, T33,A, T45,A, T57,A, T69,A, T78, " *" )'  

print, "mass:", ghead.mass, $
  format='( "* ", A,T9,E9.3, T21, E9.3, T33, E9.3, T45, E9.3, T57, E9.3, T69, E9.3, T78, " *" )'  

str = strtrim( string(ghead.npar_hw), 2)
print, "nphw:", str, $
  format='( "* ", A,T9,A, T21,A, T33,A, T45,A, T57,A, T69,A, T78, " *" )'  

print, format='( "* ", T78, " *" )'

print, "time:", ghead.time, "z:", ghead.z, "h:", ghead.h, $
  format='( "* ", A,T10,E12.5, T25,A,T32,E12.5, T50,A,T57,E12.5, T78, " *" )'

print, "boxlen:", ghead.boxlen, $
       "omegaM:", ghead.omegaM, $
       "omegaL:", ghead.omegaL, $
  format='( "* ", A,T10,E12.5, T25,A,T32,E12.5, T50,A,T57,E12.5, T78, " *" )'

print, format='( "* ", T78, " *" )'

print, "|gadget flags|",          $
  "  sfr:", ghead.flag_sfr,       $
  "  feed:", ghead.flag_feedback, $
  "  cool:", ghead.flag_cool,     $
  "  age:", ghead.flag_age,       $
  "   metal:", ghead.flag_metals,  $
  "  entr:", ghead.flag_entr_ics, $
  format='( "* ", A, 6(A,I1), T78, " *" )'




tags = tag_names(ghead)
checksphray = total( strmatch(tags,"flag_helium",/fold_case), /integer )
if (  checksphray eq 1 ) then begin
    print, "|sphray flags|",         $
      "  Hmf:", ghead.flag_Hmf,      $
      "  Hemf:", ghead.flag_Hemf,    $
      "  He:", ghead.flag_helium,    $
      "    G_HI:", ghead.flag_gammaHI, $
      "  cldy:", ghead.flag_cloudy,  $
      "   eos:", ghead.flag_eos,      $
      format='("* ", A, 6(A,I1), T78, " *")'

    print, format='( "* ", T78, " *" )'

    print, "nrays:", ghead.nrays, $
      "omegaB:", ghead.omegaB, $
      format='( "* ", A,T10,I12, T25,A,T32,E12.5, T78, " *" )'


endif






print, "* ", format='( A, 75("-"), " *" )'   

end
