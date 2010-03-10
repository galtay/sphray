; read cloudy table for redshifts
;----------------------------------
ion_table_file   = "../../ionization_tables/h1.hdf5"
ion_table_data   = h5_parse(ion_table_file,/read)
ion_table_header = ion_table_data.header
bckgnd_spectrum  = ion_table_header.spectrum
zfile            = bckgnd_spectrum.redshift._data
ghi              = bckgnd_spectrum.gammahi._data

numin = 1.0d0
numax = 5.0d0


file = 'hm01_table.txt'

openw, lun, file, /get_lun

printf, lun, 'HM 2001 Quasars + Galaxies Table'
printf, lun
printf, lun, 'z : GHI_table : GHI_int : GHI_table / GHI_int : n_gamma'

for i = 0, n_elements(zfile)-1 do begin
   make_sphray_source_from_iontable, zfile[i], numin, numax, ngamma, gammaHI, /doplot, /dofile
   printf,  lun, string(zfile[i], format='(F4.2)' ), $
           ghi[i], gammaHI, ghi[i]/gammaHI, ngamma
endfor

free_lun, lun


end


