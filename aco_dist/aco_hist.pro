;read in the file
flines=file_lines('tab.dat')
tDat_raw=fltarr(2,flines)
openr, lun, 'tab.dat', /get_lun
readf, lun, tDat_raw
close, lun
free_lun, lun

;make the distribution data
distro=fltarr(total(tDat_raw(1,*)))
plc=0
for i=0, flines-1 do begin

   distro(plc:plc+tDat_raw(1,i)-1)=tDat_raw(0,i)
   plc+=tDat_raw(1,i)

endfor

cghistoplot, distro, binsize=0.1, mininput=0.1, maxinput=100

end
