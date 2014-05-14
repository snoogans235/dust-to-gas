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
  chnsz=50000.
  seg=10000
  tol = 0.01
  brat=10e10
  bval=10e-10
  chn_i = 0.
  chain=fltarr(sz(1),sz(2),chnSz)  
  xarr=fltarr(sz(1),sz(2),seg)
  for i=0,sz(1)-1 do begin
    for j=0,sz(2)-1 do xarr(i,j,*)=findgen(seg)
  endfor

  !p.multi=[0,2,1]

  while abs(mean(brat, /nan)) gt tol do begin

    ;create n+1 point and check to make sure between 0.01 and 100
    at = ai + siga*randomn(x,[sz(1),sz(2)])
    
    ;check to make sure above lower bound
    lw=where(at(msk) lt 0.01,lwsz)
    scale=1.
    while lwsz gt 0 do begin
      fill=ai+scale*siga*randomn(y,lwsz)
      for i=0, lwsz-1 do at(msk(lw(i)))=fill(i)
      lw=where(at(msk) lt 0.01,lwsz)
      scale+=scale*siga
    endwhile

    ;check to make sure below upper bound
    hg=where(at(msk) gt 100, hgsz)
    scale=1.
    while hgsz gt 0 do begin
      fill=ai+scale*siga*randomn(z,hgsz)
      for i=0, hgsz-1 do at(msk(hg(i)))=fill(i)
      hg=where(at(msk) gt 100, hgsz)
      scale+=scale*siga
    endwhile

                                ;it would be interesting to break up the galaxy
                                ;into 2x2 or 4x4 regions and then use those to
                                ;determine if any of the regions are out of
                                ;whack and save the good ones rather than
                                ;reseting the entire galaxy.

    ;calculate the sprad
;    vari = biweight_mean(d(msk) / (hi(msk) + ai(msk) * co(msk)),sigmai)
;    vart = biweight_mean(d(msk) / (hi(msk) + at(msk) * co(msk)),sigmat)
    dgri = d(msk) / (hi(msk) + ai(msk) * co(msk))
    dgrt = d(msk) / (hi(msk) + at(msk) * co(msk))

    vari = biweight_mean(dgri,sigmai)
    vart = biweight_mean(dgrt, sigmat)

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

      ;fit a line to the data to determine whether the line chain has converged
      if plc gt 0 and plc mod(seg) eq 0 then begin
        for i=0, sz(1)-1 do begin
          ht=where(finite(d(i,*)) eq 1, htsz)
          for j=0, htsz-1 do begin
            fit=linfit(xarr(i,ht(j),*), chain(i,ht(j),plc-seg:plc-1))
            brat=fit(0)/bval
            bval=fit(0)
          endfor
       endfor
        print, mean(brat,/nan), fit(1), mean(chain(71,71,plc-seg:plc-1),/nan), sigmat
        plot, chain(71,71,*), yrange=[0.01, 100], /ylog
        plot, alog10(co(msk)/hi(msk)), alog10(dgrt), psym=5, xrange=[-2,2]
      endif
    endif

    ;reset the chain size if too large
    if chn_i ge chnsz-1 then break;chn_i=0.

  endwhile

  ;print, 'Fraction of steps:  ' + string(chnSz / num, format='(F6.4)')
  return, chain

!p.multi=[0]

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

;determine mhi
mhi=hi_mass(ihi, hdrmhi, 24.9^2)

mask=where(finite(md) ne 1, nel)
md(mask)=!values.f_nan
ihi(mask)=!values.f_nan
ico(mask)=!values.f_nan

;run the mcmc chain
chain = mcmc(fltarr(sz(1),sz(2))+10., .01, md, mhi, ico)
;chain = mcmc(chain(*,*,49999), 0.001, md, mhi, ico)
dgr = md / (mhi - chain(*,*,n_elements(chain(1,1,*))-1)*ico)

stop

plot, alog10(ico/mhi), alog10(dgr), psym=5, xrange=[-1,1.5], xtitle='Log(I!ICO!N \ !4R!3!IHI!N)', ytitle='Log(DGR)'

stop

;set up the density color palette
cgLoadCT, 33
TVLCT, cgColor('grey', /Triple), 0
TVLCT, r, g, b, /Get
palette = [ [r], [g], [b] ]

cgimage, alog10(dgr(58:80 ,45:78)), /axes, palette=palette, bottom=0, /keep_aspect_ratio, minvalue=min(alog10(dgr(58:80,45:78)),/nan), maxvalue=max(alog10(dgr(58:80,45:78)),/nan)

cgcolorbar, range=[min(alog10(dgr(58:80,45:78)),/nan),max(alog10(dgr(58:80,45:78)),/nan)], /vertical

stop

;for output purposes it will be useful to show a_co, dgr, and h_2 mass

end
