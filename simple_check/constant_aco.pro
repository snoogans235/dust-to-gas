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
pro constant_aco, mdp, mhip, icop

  ;read in the files
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
  mhi(mask)=!values.f_nan
  ico(mask)=!values.f_nan
  mask=where(finite(md) eq 1)

  aco_num=1000
  aco_min=0.01
  aco_max=100.
  aco=findgen(aco_num)*(aco_max - aco_min) / (aco_num - 1) + aco_min
  mean_g = fltarr(aco_num)
  sigm_g = fltarr(aco_num)
  mean_T = fltarr(aco_num)
  sigm_T = fltarr(aco_num)

  ;make a dgr map
  for i = 0, aco_num - 1 do begin
    dgr = alog10(md / (mhi + ico * aco(i)))
    mean_g(i) = mean(dgr(mask),/nan)
    sigm_g(i) = sqrt(variance(dgr(mask),/nan))
    mean_t(i) = biweight_mean(dgr(mask),sig)
    sigm_t(i) = sig
  endfor

  aco_bst=aco(where(sigm_g eq min(sigm_g)))
  dgr = md / (mhi + ico*aco_bst(0))

  ;plot map
  cgLoadCT, 33
  TVLCT, cgColor('grey', /Triple), 0
  TVLCT, r, g, b, /Get
  palette = [ [r], [g], [b] ]
 
  set_plot, 'ps'
  device, filename='constant_aco.ps', /inches, xsize=10, ysize=10

    cgimage, alog10(dgr[58:80,45:78]), /axes, palette=palette, bottom=0, scale=1, minValue=min(alog10(dgr),/nan), maxvalue=max(alog10(dgr),/nan), /keep_aspect_ratio, title='Log(DGR)'

    dgrcon=cgconlevels(alog10(dgr[58:80,45:78]),nlevels=10,minvalue=min(alog10(dgr),/nan))
    cgcontour, alog10(dgr[58:80,45:78]), levels=dgrcon, /onimage, label=0

    cgcolorbar, range=[min(alog10(dgr[58:80,45:78]),/nan), max(alog10(dgr),/nan)], TLocation='right', /vertical;, position=[0.29, 0.05, 0.31, 0.95]

  device,/close

  !p.multi=[0,2,1]
  device, filename='dgr_plots.ps', /inches, xsize=10, ysize=6
    
    cgplot, alog10(ico/mhi), alog10(dgr), ytitle='Log(DGR)', xtitle='Log(I!ICO!N / !4R!3!IHI!N)', psym=5

    cgplot, alog10((findgen(100)+1)/100), fltarr(100)+alog10(mean(dgr,/nan)), color='red', /overplot

    cgplot, aco, sigm_g, xtitle='!4a!3!ICO!N', ytitle='!4r!3!IDGR!N'

  device,/close
  !p.multi=0
  set_plot, 'x'

  print, 'Tukey Mean:', aco(where(sigm_t eq min(sigm_t)))
  print, 'Tukey Standard dev:', sigm_t(where(sigm_t eq min(sigm_t)))
  print, 'Gaussian Mean:', aco_bst
  print, 'Gaussian Standard dev:', sigm_g(where(sigm_g eq min(sigm_t)))

end
