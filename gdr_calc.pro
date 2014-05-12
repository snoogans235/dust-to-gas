function hi_mass, ihi, hdr, beam

dist=9.4  ;Mpc

;convert summ from Jy/beam to just good ole jansky, beam size of things is 
ppb=1.133*beam/(float(sxpar(hdr, 'CDELT1'))*3600)^2
summ=ihi/ppb
mass=2.36e5*dist^2*summ
sd = mass / (dist * sxpar(hdr, 'CDELT1')*!pi/180.)^2 / 1e12

return, sd

end

;*****************************************************************
function mcmc, ai, siga, d, hi, co

  sz=size(d)
  msk = where(finite(d) eq 1, nel)
  chnsz=50000
  tol = 0.001
  brat=fltarr(sz(1),sz(2))+10e10
  brat(where(finite(d) ne 1)) = -1*!values.f_nan
  bval=fltarr(sz(1),sz(2))+10e-10
  bval(where(finite(d) ne 1)) = -1*!values.f_nan
  chn_i = 0  
  chain=fltarr(sz(1),sz(2),chnSz)  
  xarr=fltarr(sz(1),sz(2),1000)
  for i=0,sz(1)-1 do begin
    for j=0,sz(2)-1 do xarr(i,j,*)=findgen(1000)
  endfor

  while abs(mean(brat, /nan)) gt tol do begin

    ;create n+1 point between 0.01 and 100
    at = ai + siga*randomn(x,[sz(1),sz(2)])
    stop

    ;calculate the sprad
    vari = biweight_mean(d(msk) / (hi(msk) + ai(msk) * co(msk)),sigmai)
    vart = biweight_mean(d(msk) / (hi(msk) + at(msk) * co(msk)),sigmat)

    ;set up limits for minimum
    min_tst = exp((sigmai - sigmat)/2.)
    alph = min([1., min_tst])
    u=randomu(z)

    ;if values are good then save them if not repeat step
    if u le alph then begin
      ai = at
      chain(*,*,chn_i)=ai
      ++chn_i
      plc=chn_i-1
    endif 

    ;fit a line to the data to determine whether the line chain has converged
    if plc gt 0 and plc mod(1000) eq 0 then begin
      for i=0, sz(1)-1 do begin
        ht=where(finite(d(i,*)) eq 1, htsz)
        for j=0, htsz-1 do begin
          fit=linfit(xarr(i,ht(j),*), chain(i,ht(j),plc-1000:plc-1))
          brat(i,ht(j))=fit(0)/bval(i,ht(j))
          bval(i,ht(j))=fit(0)
        endfor
      endfor
    endif

  endwhile

  ;print, 'Fraction of steps:  ' + string(chnSz / num, format='(F6.4)')
  return, chain

end


;*****************************************************************
pro gdr_main, mdp, mhip, icop

;read in files
md=mrdfits(mdp, 0, hdrmd)
ihi=mrdfits(mhip, 0, hdrmhi)
ico=mrdfits(icop, 0, hdrico)
sz=size(md)

;dust map has two pixels that are out of place
md(45:55 ,120:130)=-1*!values.f_nan

mask=where(finite(md) ne 1, nel)
md(mask)=!values.f_nan
ihi(mask)=!values.f_nan
ico(mask)=!values.f_nan

;determine mhi
mhi=hi_mass(ihi, hdrmhi, 24.9^2)

;run the mcmc chain
chain = mcmc(fltarr(sz(1),sz(2))+5., 1, md, mhi, ico)
;chain = mcmc(chain(*,*,49999), 0.001, md, mhi, ico)
dgr = md / (mhi - chain(*,*,n_elements(chain(1,1,*))-1)*ico)

stop

plot, alog10(ico/ihi), alog10(dgr), psym=5, xrange=[-1,1.5], xtitle='Log(I!ICO!N \ !4R!3!IHI!N)', ytitle='Log(DGR)'

stop

;set up the density color palette
cgLoadCT, 33
TVLCT, cgColor('grey', /Triple), 0
TVLCT, r, g, b, /Get
palette = [ [r], [g], [b] ]

cgimage, alog10(dgr(58:80 ,45:78)), /axes, palette=palette, bottom=0, /keep_aspect_ratio, minvalue=min(alog10(dgr(58:80,45:78)),/nan), maxvalue=max(alog10(dgr(58:80,45:78)),/nan)

cgcolorbar, range=[min(alog10(dgr(58:80,45:78)),/nan),max(alog10(dgr(58:80,45:78)),/nan)], /vertical

stop

end
