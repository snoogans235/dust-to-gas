function gamma_dist, X, P

;X=[alpha, beta]

gam_dist= 1. / (gamma(1+P(0)) * P(1)^(1+P(1))) * X^(P(0)) * exp(-1*X / P(1))

;normalize
gam_dist = gam_dist / total(gam_dist)

return, gam_dist

end

;*************************************************************
function exponential_dist, X, P

;x=[th]

exp_dist= exp(-1* X / P(0) ) / P(0)

;normalize
exp_dist = exp_dist / total( exp_dist)

return, exp_dist

end

;*************************************************************
function maxwell_dist, X, P

;x=[a] (energy)
max_dist=sqrt(2/!pi) * X^2 * exp(-1* X^2 / (2*P(0)^2)) / P(0)^3

;normalize
max_dist = max_dist / total(max_dist)

return, max_dist

end

;*************************************************************
function comp_simp_int, X, Y

sum=0.

for i=0, n_elements(Y)-2 do begin
  if i mod 2 eq 0 then sum+=2*Y(i)
  if i mod 2 ne 0 then sum+=4*Y(i)
endfor

sum+=Y(0) + Y(i)
area= 2.*(X(i) - X(0)) / i * sum

return, area

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
histo=float(cghistogram(distro, binsize=0.1, min=0.1, max=30))
histo=float(histo) / total(histo,/nan)
cgplot, ind_var, histo, psym=10

;try one: gamma distribution
parinfo=replicate({limited:[0,0], limits:[0.,0.]},2)
parinfo[*].limited(0)=1.
parinfo[0].limits(0)=-1
parinfo[1].limits(0)=0.

gam_fit = mpfitfun('gamma_dist', ind_var, histo, 1/sqrt(histo), [1., 1.], bestnorm=bestnorm, perror=perror, /nan, /quiet)

gam_build=gamma_dist(ind_var, gam_fit)
cgplot, ind_var, gam_build, color='red', /overplot

print, 'Gamma Distribution: ', bestnorm, gam_fit

;try two: exponential distribution
parinfo=replicate({limited:[0,0], limits:[0.,0.]},1)
parinfo.limited(0)=1
parinfo.limits(0)=0.

exp_fit = mpfitfun('exponential_dist', ind_var, histo, 1/sqrt(histo), [1.], bestnorm=bestnorm, perror=perror, /nan, /quiet)

exp_build = exponential_dist(ind_var, exp_fit)
cgplot, ind_var, exp_build, color='blue', /overplot
print, 'Exponential Distribution: ', bestnorm, exp_fit

;try three: boltzman distribution
parinfo=replicate({limited:[0,0], limits:[0.,0.]},1)
parinfo.limited(0)=1
parinfo.limits(0)=0.

max_fit = mpfitfun('maxwell_dist', ind_var, histo, 1/sqrt(histo), [1.], bestnorm=bestnorm, perror=perror, /nan, /quiet)

max_build=maxwell_dist(ind_var, max_fit)
cgplot, ind_var, max_build, color='dark green', /overplot
print, 'Maxwell Distribution: ', bestnorm, max_fit

stop

end
