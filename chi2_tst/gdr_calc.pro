function hi_mass, ihi, hdr, beam

dist=9.4  ;Mpc

;convert summ from Jy/beam to just good ole jansky
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

;  !p.multi=[0,2,1]

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
      hg=where(at(msk) gt 100, hgsz1)
      scale+=scale*siga
    endwhile

                                ;it would be interesting to break up the galaxy
                                ;into 2x2 or 4x4 regions and then use those to
                                ;determine if any of the regions are out of
                                ;whack and save the good ones rather than
                                ;reseting the entire galaxy.

    ;set the unneeded pixels to blank values
    at(where(finite(d) eq 0))=!values.f_nan
    ai(where(finite(d) eq 0))=!values.f_nan

    ;calculate the spread
    dgri = d(msk) / (hi(msk) + ai(msk) * co(msk))
    dgrt = d(msk) / (hi(msk) + at(msk) * co(msk))

    vari = biweight_mean(dgri(msk),sigmai)
    vart = biweight_mean(dgrt(msk),sigmat)

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
        ;print, mean(brat,/nan), fit(1), mean(chain(71,71,plc-seg:plc-1),/nan), sigmat
        ;plot, chain(71,71,*), yrange=[0.01, 100], /ylog
        ;plot, alog10(co(msk)/hi(msk)), alog10(dgrt), psym=5, xrange=[-2,2]
      endif
    endif

    ;reset the chain size if too large
    if chn_i gt chnsz-1 then break;chn_i=0.

  endwhile

  ;print, 'Fraction of steps:  ' + string(chnSz / num, format='(F6.4)')
  return, chain

;!p.multi=[0]

end

;*****************************************************************
pro dgr_output, dgr, aco, mh2, ico_shi, chain


;set up the density color palette
cgLoadCT, 33
TVLCT, cgColor('grey', /Triple), 0
TVLCT, r, g, b, /Get
palette = [ [r], [g], [b] ]

set_plot, 'ps'

device, filename='mcmc_check.ps', /inches, xsize=9, ysize=9
  !p.multi=[0,2,2]

  !p.position=[0.15, 0.1, 0.5, 0.45]
  mean=biweight_mean(alog10(dgr(where(finite(dgr) eq 1))),sigma)
  cgplot, alog10(ico_shi), alog10(dgr), psym=5, xtitle='Log(I!ICO!N \ !4R!3!IHI!N)', ytitle='Log(DGR)'
  cgplot, alog10(ico_shi), fltarr(5000)+mean, linestyle=0, color='red', /overplot
  cgplot, alog10(ico_shi), fltarr(5000)+mean-sigma, linestyle=2, color='red', /overplot
  cgplot, alog10(ico_shi), fltarr(5000)+mean+sigma, linestyle=2, color='red', /overplot

  !p.position=[0.6, 0.1, 0.95, 0.45]
  cgplot, findgen(n_elements(chain)) * 1e-4, chain, ytitle='!4a!3!ICO!N', xtitle='Chain Number x 10!E4!N';, yrange=[0.01, 100], /ylog

  !p.position=[0.15, 0.6, 0.5, 0.95]
  cghistoplot, aco, nbins=25, xtitle='!4a!3!ICO!N', ytitle='Frequency', /nan
 
  !p.position=[0.6, 0.6, 0.95, 0.95]
  cghistoplot, alog10(dgr), nbins=25, xtitle='Log(DGR)', ytitle='Frequency', /nan
  

device, /close

!p.multi=[0,3,1]
device, filename='gdr_output.ps', /inches, xsize=14, ysize=5

  ;will need to add ra and dec to axis

  ;plot the gas to dust ratio
  !p.position=[0.02, 0.05, 0.27, 0.95] 
  cgimage, alog10(dgr), /axes, palette=palette, bottom=0, scale=1, minValue=min(alog10(dgr),/nan), maxvalue=max(alog10(dgr),/nan), /keep_aspect_ratio, title='Log(DGR)'

  dgrcon=cgconlevels(alog10(dgr),nlevels=10,minvalue=min(alog10(dgr),/nan))
  cgcontour, alog10(dgr), levels=dgrcon, /onimage, label=0

  cgcolorbar, range=[min(alog10(dgr),/nan), max(alog10(dgr),/nan)], TLocation='right', /vertical, position=[0.29, 0.05, 0.31, 0.95]

  ;plot aco values
  !p.position=[0.33, 0.05, 0.58, 0.95]
   cgimage, aco, /axes, palette=palette, bottom=0, scale=1, minValue=min(aco,/nan), maxvalue=max(aco, /nan), /keep_aspect_ratio, title='!4a!3!ICO!N'

   acocon=cgconlevels(aco,nlevels=10,minvalue=min(aco,/nan))
   cgcontour, aco, levels=acocon, /onimage, label=0

   cgcolorbar, range=[min(aco,/nan), max(aco,/nan)], TLocation='right', /vertical, position=[0.60,0.05,0.62,0.95], title='M!I!9!Z(6E)!X!N pc!E-2!N (K km s!E-1!N)!E-1!N'

   ;plot mh2 values
   !p.position=[0.68, 0.05, 0.91,0.95]
   cgimage, mh2, /axes, palette=palette, bottom=0, scale=1, minValue=min(mh2,/nan), maxvalue=max(mh2,/nan), title='M!IH!I2!N'

   mh2con=cgconlevels(mh2, nlevels=10,minvalue=min(mh2,/nan))
   cgcontour, mh2, levels=mh2con, /onimage, label=0

   cgcolorbar, range=[min(mh2,/nan),max(mh2,/nan)], TLocation='right', /vertical, position=[0.95, 0.05, 0.97,0.95], title='M!I!9!Z(6E)!X!N pc!E-2!N'

device, /close
set_plot, 'x'
!p.multi=0

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
chain = mcmc(fltarr(sz(1),sz(2))+50., .01, md, mhi, ico)
aco = mean(chain(*,*,n_elements(chain(1,1,*))-10000:n_elements(chain(1,1,*))-1),dimension=3)
dgr = md / (mhi + aco*ico)

dgr_output, dgr[58:80,45:78], aco[58:80,45:78], aco[58:80,45:78]*ico[58:80,45:78], ico[58:80,45:78] / mhi[58:80,45:78], chain(71,71,*)

end
