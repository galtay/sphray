function gadget_header_define, flag_sphray=flag_sphray
;+
; NAME: 
;    GADGET_HEADER_DEFINE
;
; PURPOSE: 
;    Creates an empty header structure.
;
; CATEGORY:
;    IDL Hands Library - Gadget Routines
;
; CALLING SEQUENCE: 
;    Result = gadget_header_define( [ /flag_sphray ] )
;
; KEYWORD PARAMETERS:
;    /flag_sphray = If set, then the additional variables omegaB,
;    nrays, flag_Hmf, flag_Hemf, flag_helium, flag_gammaHI,
;    flag_cloudy, and flag_eos are defined as part of the header.
;
; OUTPUTS:
;    Returns an empty header structure.
;
; EXAMPLE:
;    ghead = gadget_header_define()    
;
; MODIFICATION HISTORY:
;    Version 1, Gabriel Altay, CMU, 14 Oct 2009.
;    Version 2, Gabriel Altay, ICC, 3 Feb 2010.
;-



; check keyword arguments
;---------------------------

if (n_elements(flag_sphray) eq 0) then flag_sphray=0

; define header
;---------------

if (flag_sphray) then begin
    
    header = create_struct( name="sphray_header",         $
                            "npar_file",     lonarr(6),   $
                            "mass",          dblarr(6),   $
                            "time",          0.0d0,       $
                            "z",             0.0d0,       $  
                            "flag_sfr",      0L,          $
                            "flag_feedback", 0L,          $
                            "npar_snap",     lonarr(6),   $
                            "flag_cool",     0L,          $
                            "nfiles",        0L,          $
                            "boxlen",        0.d0,        $
                            "omegaM",        0.d0,        $
                            "omegaL",        0.d0,        $
                            "h",             0.d0,        $
                            "flag_age",      0L,          $
                            "flag_metals",   0L,          $
                            "npar_hw",       lonarr(6),   $
                            "flag_entr_ics", 0L,          $
                            "omegaB",        0.d0,        $
                            "nrays",         0LL,         $
                            "flag_Hmf",      0L,          $
                            "flag_Hemf",     0L,          $                    
                            "flag_helium",   0L,          $
                            "flag_gammaHI",  0L,          $
                            "flag_cloudy",   0L,          $
                            "flag_eos",      0L,          $
                            "unused",        fltarr(5)    )

  endif else begin

     header = create_struct( name="gadget_header",        $
                            "npar_file",     lonarr(6),   $
                            "mass",          dblarr(6),   $
                            "time",          0.0d0,       $
                            "z",             0.0d0,       $  
                            "flag_sfr",      0L,          $
                            "flag_feedback", 0L,          $
                            "npar_snap",     lonarr(6),   $
                            "flag_cool",     0L,          $
                            "nfiles",        0L,          $
                            "boxlen",        0.d0,        $
                            "omegaM",        0.d0,        $
                            "omegaL",        0.d0,        $
                            "h",             0.d0,        $
                            "flag_age",      0L,          $
                            "flag_metals",   0L,          $
                            "npar_hw",       lonarr(6),   $
                            "flag_entr_ics", 0L,          $
                            "unused",        fltarr(15)   )

  endelse

  return, header


END
