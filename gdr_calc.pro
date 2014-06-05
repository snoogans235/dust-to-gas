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
function gal_grid, img, grd_sz

  sz=size(img)
  grd=create_struct('grd_num', 0L, 'plc_vls', fltarr(grd_sz(0)*grd_sz(1)))  
  grd_app=create_struct('grd_num', 0L, 'plc_vls', fltarr(grd_sz(0)*grd_sz(1)))

  ;will go through an image and generate a grid with a pixel overlay
  overlay=0
  if overlay ge min(grd_sz,/nan) then begin
    print, 'Grid is smaller than overlay.  Try increasing grid size.'
    grd.grd_num=!values.f_nan
    return, grd
  endif

;plt_msk_regi=intarr(sz(1),sz(2))
;plt_msk_full=intarr(sz(1),sz(2))

  ;since idl's where function returns a 1D array, I will need to have the grid 
  ;values represent their incremental position in the array, this way I can
  ;group what pixels go into each grid segment and identify them by their
  ;1D position value
  grd_cnt=findgen(sz(1),sz(2))
  grd_cnt(where(finite(img) eq 0)) = !values.f_nan
  ;grd_cnt = grd_cnt - min(grd_cnt, /nan)
  
  ;use a counter for the grid numbers
  grd_num = 0

  ;establish crazy high top left corner value
  bl_corn=[10e10, 10e10]

  ;run through the image rows and start to lay down the grid
  for i=0, sz(2)-1 do begin

    ;check to see if any values are in the selected rows
    ht=where(finite(grd_cnt(*,i)) eq 1,htsz)

    for j=0, htsz-1, grd_sz(1)-overlay do begin

      ;set up a counter to know where the bottom left corner is at.  So every
      ;grd_sz - 1 rows a new box will be made if the x position of the top left
      ;corner is the same (i.e. the min value of hit doesn't change)
      if min(ht, /nan) ne bl_corn(0) or i ge bl_corn(1) or cont_flag eq 1 then begin
       if j eq 0 then cont_flag=1
        ;reset bottom left corner
        bl_corn=[ht(0), i+grd_sz(0)-1]

        ;generate an array that contains the one dimensional position value 
        ;for pixels in the target image
        elem=fltarr(grd_sz(0)*grd_sz(1))
        for l = 0, grd_sz(0)-1 do elem(l*grd_sz(1):l*grd_sz(1)+grd_sz(1)-1)=grd_cnt(ht(j),i+l)+findgen(grd_sz(1))

        ;filter out nan's from the grid and any check to see if the grid
        ;has enough elements.  A grid us useable if the number of elements is 
        ;greater than the number of elements of the grid size.  So for a 4x4
        ;grid the region would need more than 4 elements
        nan_flag = where(finite(img(elem)) eq 1, nan_flag_sz)
        

        ;append the structures to each other.  If the first entry, will need 
        ;to fill grd first
        if nan_flag_sz gt min(grd_sz,/nan) then begin
          if grd_num eq 0 then begin
            grd.grd_num=grd_num++
            grd.plc_vls[0:nan_flag_sz-1]=elem(nan_flag)
            if nan_flag_sz ne grd_sz(0)*grd_sz(1) then grd[grd_num-1].plc_vls(nan_flag_sz:grd_sz(0)*grd_sz(1)-1)=!values.f_nan
          endif else begin
            grd_app.grd_num=grd_num++
            grd_app.plc_vls[0:nan_flag_sz-1]=elem(nan_flag)
            if nan_flag_sz ne grd_sz(0)*grd_sz(1) then grd_app.plc_vls(nan_flag_sz:grd_sz(0)*grd_sz(1)-1)= !values.f_nan
            grd=struct_append(grd, grd_app)
          endelse
;plt_msk_full(grd[grd_num-1].plc_vls)=1 
;plt_msk_regi(grd[grd_num-1].plc_vls)=17
        endif

       

;cgLoadCT, 33
;TVLCT, cgColor('grey', /Triple), 0
;TVLCT, r, g, b, /Get
;palette = [ [r], [g], [b] ]

;cgimage, img[58:80,43:78], /axes, palette=palette, bottom=0, scale=1, minValue=min(img,/nan)-0.1*min(img,/nan), maxvalue=max(img,/nan), /keep_aspect_ratio, oposition=oposi
;cgimage, plt_msk_full[58:80,43:78], transparent=50, alphafgpos=oposi, minValue=-2, maxValue=-2
;cgimage, plt_msk_regi[58:80,43:78], transparent=50, alphafgpos=oposi, minValue=-2, maxValue=20

;plt_msk_regi(grd[grd_num-1].plc_vls)=0
;stop

      endif 

      if j gt 0 and cont_flag eq 0 then j = htsz
  
    endfor
    cont_flag=0
  endfor

  return, grd

end

;*****************************************************************
function mcmc, ai, siga, d, hi, co, grid, mean

  ;it might be less intense memory wise to keep the random numbers cut down to
  ;one element arrays

  ;set up the bounds for desired distribution
  mu=0.
  sig=1e-5

  sz=size(d)
  ;grdsz=max(grid.grd_num)-1
  flag=intarr(sz(1),sz(2))
  msk = where(finite(d) eq 1, nel)
  chnsz=50000.
  chain=fltarr(sz(1),sz(2),chnSz)
  chn_i = 0L
  badvals=fltarr(sz(1),sz(2))+1
  brat=10
  tol=1
  cnt=0L

  while chn_i lt chnsz do begin

    ;create n+1 point and check to make sure between 0.01 and 100
     
;since I am begin caught in the limit check of the code it might be worth while
;to look into making my own random number generator with boundaries.  Why this
;isn't an option in idl is beyond me...
   at = ai + badvals*siga*randomn(x,[sz(1),sz(2)])
;    at = ai + badvals*siga*randomn(x,[sz(1),sz(2)], gamma(mean))-mean

    
    ;check to make sure above lower bound
    ;put low and high check into one loop
    lw=where(at(msk) lt 0.01, lwsz)
    hg=where(at(msk) gt 100, hgsz)
    badval=lwsz+hgsz
    whilecnt=0

    while badval gt 0 do begin
     if hgsz gt 0 then fill_hg=ai+siga*randomn(z1,hgsz)
     if lwsz gt 0 then fill_lw=ai+siga*randomn(z2,lwsz)

     if lwsz gt hgsz then loop_cnt = lwsz else loop_cnt=hgsz
     for i=0, loop_cnt - 1 do begin
       if i lt lwsz then at(msk(lw(i)))=fill_lw(i)
       if i lt hgsz then at(msk(hg(i)))=fill_hg(i)
     endfor
if whilecnt gt 10 then stop
     lw=where(at(msk) lt 0.01, lwsz)
     hg=where(at(msk) gt 100, hgsz)
     badval=lwsz+hgsz
     ++whilecnt
   endwhile

    ;calculate the variance for each grid region and determine which grid
    ;regions are useable.
    dgri = d / (hi + ai * co)
    dgrt = d / (hi + at * co)
    
    ;calculate the variance for the entire galaxy
    avg_i=biweight_mean(dgri(msk),sig_i)
    avg_t=biweight_mean(dgrt(msk),sig_t)
 
    ;determine if the variance is less than what should be excuded
    min_tst = exp((-1*(sig_t)^2 + sig_i^2 + 2 * (sig_t - sig_i)*mu) / (2 * sig^2))
    alph = min([1,min_tst])
    u=randomu(z)

    ;need to look and see if at is worthy of keeping.  If not then go into
    ;grid_tst and determine a way to select what is a good grid region and what
    ;isn't.  It might be wise to look at whether the mean falls in a certain
    ;range?  No this will bias the data towards the initial condition.  It 
    ;would be useful to look at the dgr variance in each grid region and if the
    ;variance is too large in that case then try to reposition?
    
    ;reverse conditions for saving.  This will begin determining which grid
    ;cells should be changed.
    if u le alph then begin
      ai=at
      chain(*,*,chn_i)=ai
      ++chn_i
    endif else begin
      ;instead of doing the mcmc check in grid_test it would be useful to pick
      ;out where the values are less and where they are more then send those
      ;into the grid to be changed?  I would need some sort of map of which
      ;pixels to change then feed that into the random generation.  My only
      ;fear is that this will build towards finding a mean of the starting
      ;value.

      ;look at the variance in each dgri grid location
      badvals = grid_tst(dgrt, grid)

      ;gtavg=where(dgri ge avg_i + sig_i)
      ;ltavg=where(dgri le avg_i - sig_i)
      ;flag(gtavg)=1
      ;flag(ltavg)=1

      ;badvals = grid_tst(flag, grid)
    endelse

    ++cnt
    if chn_i mod 1000 eq 0 then print, chn_i, cnt

  endwhile
stop
  print, 'Fraction of steps:  ' + string(1. * chn_i / cnt, format='(F6.4)')
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
function grid_tst, img, grid

  ;should return a 2d array that will have which pixels need to be generated
  img_sz = size(img)
  grd_el = n_elements(grid)
  vari = fltarr(grd_el)
  flag=intarr(img_sz(1), img_sz(2))

  for i=0, grd_el - 1 do vari(i)=variance(img(grid[i].plc_vls),/nan)

  var_avg=mean(vari,/nan)
  var_sig=sqrt(variance(vari,/nan))

  ht_hi=where(vari gt var_avg+var_sig, hisz)
  ht_lo=where(vari lt var_sig-var_avg, losz)

  for i=0, hisz-1 do flag(grid[ht_hi(i)].plc_vls)=1
  for i=0, losz-1 do flag(grid[ht_lo(i)].plc_vls)=1

  flag(where(finite(flag) eq 0))=!values.f_nan

  return, flag

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

;establish a grid to use
grid=gal_grid(md, [3,3]) ;[y,x]
if grid[0].grd_num eq long(!values.f_nan) then stop

;run the mcmc chain
chain = mcmc(fltarr(sz(1),sz(2))+1.25, 0.1, md, mhi, ico, grid)
;chain = mcmc(chain(*,*,n_elements(chain(1,1,*))-1),0.01, md, mhi, ico)
aco = mean(chain(*,*,n_elements(chain(1,1,*))-10000:n_elements(chain(1,1,*))-1),dimension=3)
dgr = md / (mhi + aco*ico)
aco(mask)=!values.f_nan

dgr_output, dgr[58:80,45:78], aco[58:80,45:78], aco[58:80,45:78]*ico[58:80,45:78], ico[58:80,45:78] / mhi[58:80,45:78], chain(71,71,*)

end
