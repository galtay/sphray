pro altay_oplot, x, y, ct=ct, _extra=extra

on_error,2


if (!D.name eq "PS") then begin

    thick=6.0

endif else if (!D.name eq "X") then begin

    thick=2.0

endif else begin

    message, "  Plotting device not recognized"

endelse


if n_elements(ct) ne 0 then loadct, ct
oplot, x, y, thick=thick, _strict_extra=extra


end
