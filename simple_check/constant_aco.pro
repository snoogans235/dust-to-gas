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
pro constant_aco, mdp, mhip, icop, aco

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
  ihi(mask)=!values.f_nan
  ico(mask)=!values.f_nan

  ;make a dgr map
  dgr = md / (ihi + ico * aco)

  ;plot map
  cgLoadCT, 33
  TVLCT, cgColor('grey', /Triple), 0
  TVLCT, r, g, b, /Get
  palette = [ [r], [g], [b] ]

  cgimage, alog10(dgr[58:80,45:78]), /axes, palette=palette, bottom=0, scale=1, minValue=min(alog10(dgr),/nan), maxvalue=max(alog10(dgr),/nan), /keep_aspect_ratio, title='Log(DGR)'

  dgrcon=cgconlevels(alog10(dgr[58:80,45:78]),nlevels=10,minvalue=min(alog10(dgr),/nan))
  cgcontour, alog10(dgr[58:80,45:78]), levels=dgrcon, /onimage, label=0

  cgcolorbar, range=[min(alog10(dgr[58:80,45:78]),/nan), max(alog10(dgr),/nan)], TLocation='right', /vertical;, position=[0.29, 0.05, 0.31, 0.95]

end
