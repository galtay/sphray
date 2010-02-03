function gadget_header_read, snapbase, silent=silent, flag_sphray=flag_sphray
;+
; NAME:
;    gadget_header_read
;
; PURPOSE:
;    Reads gadget header(s) into a gadget header structure.
;
; CATEGORY:
;    IDL Hands Library - Gadget Routines
;
; CALLING SEQUENCE:
;    RESULT = gadget_header_read( snapbase, [ /silent, /flag_sphray ] ) 
;
; INPUTS:
;    snapbase - Full file name (including path) to the gadget
;               snapshot with no file number extension. 
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;    silent      - If set, the header information is not reported to
;                  screen. 
;    flag_sphray - If set, the extra data in a Sphray header is defined in
;                  the header structure (see gadget_header_define.pro). 
;
; OUTPUTS:
;    Header structure with data from snapbase.  Will be an
;    array of header structures for snapshots spread over multiple
;    files (e.g. ghead[3].npar_file[1] is the number of halo
;    particles in the 4th file ) 
;
; RESTRICTIONS:
;    Relys on readu not padding the header structure. 
;
; EXAMPLE:
;    snapbase = "/home/galtay/Data/DMG_N64_L1/snapshot_002"
;    gheads = gadget_header_read(snapbase, /silent)
;    gadget_header_print, gheads[4]
;
; MODIFICATION HISTORY:
;    Version 1, Gabriel Altay, CMU, 14 Oct 2009.
;-


; check arguments
;-----------------
on_error, 2

if (n_elements(snapbase) eq 0) then begin
    message, "usage: gheads = gadget_header_read(snapbase)"
endif

if (n_elements(silent) eq 0) then silent=0

if (n_elements(flag_sphray) eq 0) then flag_sphray=0

; define a scalar header structure
;------------------------------------------
tmphead = gadget_header_define(flag_sphray=flag_sphray)

; determine how many files in the snapshot
;------------------------------------------
snapfile = snapbase
openr, lun, snapfile, /f77_unformatted, /get_lun, error=err

; if passed snapbase exists as a file (single file snapshot)
;------------------------------------------------------------
if (err eq 0) then begin
    readu, lun, tmphead
    free_lun, lun
    if (not silent) then gadget_header_print, tmphead 
    return, tmphead
    
; if passed snapbase requires an extension (multiple file snapshot)
;--------------------------------------------------------------------
endif else begin
    snapfile = snapbase + ".0"
    openr, lun, snapfile, /f77_unformatted, /get_lun, error=err
    if (err ne 0) then begin
        message, 'error opening, ', /continue
        message, snapfile
    endif
    
    readu, lun, tmphead    
    free_lun, lun
    gheads = replicate(tmphead, tmphead.nfiles)
    
    for i = 0, tmphead.nfiles-1 do begin
        snapfile = snapbase + "." + strtrim(i,2)
        openr, lun, snapfile, /f77_unformatted, /get_lun, error=err
        
        if (err ne 0) then begin
            message, 'error opening, ', /continue
            message, snapfile
        endif

        readu, lun, tmphead
        free_lun, lun
        gheads[i] = tmphead
        
    endfor
    
    if (not silent) then gadget_header_print, gheads[0] 
    return, gheads
    
endelse




end
