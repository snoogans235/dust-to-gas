function gamma_dist, X, P

;X=[alpha, beta]

gam_dist=1 / (gamma(1+P(0)) * P(1)^(1+P(1))) * X^(P(0)) * exp(-1*X / P(1))

return, gam_dist

end

;*************************************************************
function exponential_dist, X, P

;x=[th]

exp_dist= exp(-1* X / P(0) ) / P(0)

return, exp_dist

end

;*************************************************************
pro aco_hist, file
;read in the file
flines=file_lines(file)
tDat_raw=fltarr(2,flines)
openr, lun, file, /get_lun
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

;check out some different distributions and their fits
ind_var=findgen(300) * (30. - 0.1) / (1000-1)+0.1
histo=cghistogram(distro, binsize=0.1, min=0.1, max=30)
histo=float(histo) / max(histo,/nan)
cgplot, ind_var, histo, psym=10

;try one: gamma distribution
parinfo=replicate({limited:[0,0], limits:[0.,0.]},2)
parinfo[*].limited(0)=1.
parinfo[0].limits(0)=-1
parinfo[1].limits(0)=0.

gam_fit = mpfitfun('gamma_dist', ind_var, histo, fltarr(300)+1., [1., 1.], bestnorm=bestnorm, perror=perror, /nan, /quiet)

gam_build=gamma_dist(ind_var, gam_fit)
cgplot, ind_var, gam_build, color='red', /overplot

print, 'Gamma Distribution: ', bestnorm, gam_fit

;try two: exponential distribuito
parinfo=replicate({limited:[0,0], limits:[0.,0.]},1)
parinfo.limited(0)=1
parinfo.limits(0)=0.

exp_fit = mpfitfun('exponential_dist', ind_Var, histo, fltarr(300)+1., [1.], bestnorm, /nan, /quiet)

exp_build = exponential_dist(ind_var, exp_fit)
cgplot, ind_var, exp_build, color='blue', /overplot
print, 'Exponential Distribution: ', bestnorm, exp_fit

end
