PRO SCREEN_to_png, filename
IF N_PARAMS() EQ 0 THEN filename="screen_dump.png"
                                 
; Get the screen dump of the current graphics window.
image = tvrd(true=1)
write_png,filename,image


END
