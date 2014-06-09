function regAn, regFil, param
;function regAn, regFil, img, pn

;the idea is to drop a circle and then run though and see which values are
;outside of the region or -nan values

;generate an array of values to that contain useable coordinates
nlines=file_lines(regFil)
corn=intarr(2, nlines)
openr, lun, regFil, /get_lun
readf, lun, corn
close, lun
free_lun, lun

;build an array that encompases all of the points in the region
ec=nlines-1

;create a region based on the minimum and maximum corners
mxCorn=[min(corn(0,1:ec)), max(corn(0,1:ec)),0]
myCorn=[min(corn(1,1:ec)), max(corn(1,1:ec)),0]
mxCorn(2) = mxCorn(1) - mxCorn(0)
myCorn(2) = myCorn(1) - myCorn(0)
valx=fltarr(mxCorn(2)+1, myCorn(2)+1)
valy=fltarr(mxCorn(2)+1, myCorn(2)+1)

for i=0, myCorn(2) do valx(*,i)=dindgen(mxCorn(2)+1)+mxCorn(0)
for i=0, mxCorn(2) do valy(i,*)=dindgen(myCorn(2)+1)+myCorn(0)

refx = dindgen(mxCorn(2)+1)+mxCorn(0)
refy = dindgen(myCorn(2)+1)+myCorn(0)

hitarr=intarr(mxCorn(2)+1, myCorn(2)+1)

;param=mrdfits(img, 0, hdr)
;sp=strsplit(pn, '_', /extract)
;sp=strsplit(img, '/', /extract)
;sp=strsplit(sp(3), '_', /extract)
;sp=strsplit(sp(3), '.', /extract)
;if sp(0) eq 'mass' then param=param / 1e5
param_mini=param(mxCorn(0):mxCorn(1), myCorn(0):myCorn(1))

;I am going to define the perimiter of each object from corn(*,i) to
;corn(*,i+1), and then randomly generate a pixel value (x_r, y_r) from a
;uniform distribution then find its corresponding y cooridinate and determine
;if the x value lies between the perimeter values.

peri=fltarr(mxCorn(2)+1, myCorn(2)+1)
;construct the perimeter of the region
for i=1, ec - 1 do begin

  case 1 of
    ;if y's are the same
    corn(1,i+1) eq corn(1,i): begin
      htx=[[where(valx(*,0) eq corn(0,i))],[where(valx(*,0) eq corn(0,i+1))]]
      hty=where(valy(0,*) eq corn(1,i))
      htx1=min(htx)
      htx2=max(htx)
      peri(htx1(0):htx2(0),hty)=10.

    end
    ;if x's are the same
    corn(0,i+1) eq corn(0,i): begin
      hty=[[where(valy(0,*) eq corn(1,i))], [where(valy(0,*) eq corn(1,i+1))]]
      htx=where(valx(*,0) eq corn(0,i))
      hty1=min(hty)
      hty2=max(hty)
      peri(htx,hty1(0):hty2(0))=10.

    end
    else: begin
      sml=min([corn(0,i), corn(0,i+1)])
      
      ;perimeter line will need to have the number of elements of the greatest
      ;corner difference

      case 1 of 
        abs(corn(0,i)-corn(0,i+1)) ge abs(corn(1,i)-corn(1,i+1)): begin
          smlel=abs(corn(0,i)-corn(0,i+1))
          mxarr=1.*(corn(0,i+1)-corn(0,i)) / smlel
          xarr = mxarr*findgen(smlel+1)+corn(0,i)
          m=1.*(corn(1,i+1)-corn(1,i))/(1.*(corn(0,i+1)-corn(0,i)))
          b=corn(1,i) - m * corn(0,i)
          line = m*xarr+b

        end

        abs(corn(0,i)-corn(0,i+1)) lt abs(corn(1,i)-corn(1,i+1)): begin
          smlel=abs(corn(1,i)-corn(1,i+1))
          mxarr=1.*(corn(0,i+1)-corn(0,i)) / smlel
          xarr = mxarr*findgen(smlel+1)+corn(0,i)
          m = (corn(1,i+1)-corn(1,i)) / (xarr(smlel) - xarr(0))
          b = corn(1,i) - m*corn(0,i)
          line = m*xarr+b     
        end
      endcase

      ;save the line values
      lnsz=n_elements(line)
      for cnt=0, lnsz-1 do begin
        htx=where(valx(*,0) ge round(xarr(cnt)) and valx(*,0) le round(xarr(cnt)))
        hty=where(valy(0,*) ge round(line(cnt)) and valy(0,*) le round(line(cnt)))
        peri(htx(0),hty(0))=10.
      endfor

    end
   endcase 

endfor

;randomly generate points to fill the area
rnd=round([[mxCorn(2)*randomn(s, 1000000)+Corn(0,0)], [myCorn(2)*randomn(s,1000000)+corn(1,0)]])

htval=where(rnd(*,0) lt mxCorn(1) and rnd(*,0) gt mxCorn(0) and rnd(*,1) lt myCorn(1) and rnd(*,1) gt myCorn(0),htvalsz)

for i=0, htvalsz-1 do begin

  ;determine the position in valx and valy of the random point
  htx=where(valx(*,0) eq rnd(htval(i),0)) ;pixel value
  hty=where(valy(0,*) eq rnd(htval(i),1)) ;pixel value
  
  htp=where(peri(*,hty) ge 10.,htpsz) ;perimeter boundaries in x
    if rnd(htval(i),0) gt valx(htp(0),hty) and rnd(htval(i),0) lt valx(htp(htpsz-1),hty) then peri(htx,hty)=1.

endfor

hits=where(peri lt 1.,htsz)
if htsz gt 0 then begin
  valx(hits)=-1*!values.f_nan
  valy(hits)=-1*!values.f_nan
endif

;account for nans in the data
hits=where(finite(param_mini) eq 0)
valx(hits)=-1*!values.f_nan
valy(hits)=-1*!values.f_nan

;make the radius array.  Will need to have the same dimensions as the image
;being compared to i.e. same dimensions as param

sz=size(param)
rad = fltarr(sz(1), sz(2))
hits=where(finite(param) eq 0)

for i = 0, sz(1)-1 do begin
  for j =0, sz(2)-1 do begin

    if i ge mxcorn(0) and i lt mxcorn(1) and j ge mycorn(0) and j lt mycorn(1) then begin
      if finite(valx(i-mxcorn(0),j-mycorn(0))) eq 0 then rad(i,j)=-1*!values.f_nan else rad(i,j)=param(i,j)
      ;rad(i,j)=sqrt((valx(i-mxcorn(0),j-mycorn(0))-corn(0,0))^2+(valy(i-mxcorn(0),j-mycorn(0))-corn(1,0))^2)

    endif else begin
      rad(i,j)=-1*!values.f_nan
    endelse
  endfor
endfor

ht=where(finite(rad(mxcorn(0):mxcorn(1),mycorn(0):mycorn(1))) eq 0)
param_mini(ht) = -1*!values.f_nan

return, rad

end

;************************************************************
function aco_calc, md, mhi, ico, regi

  aco_num=1000
  aco_min=0.01
  aco_max=100.
  aco=findgen(aco_num)*(aco_max - aco_min) / (aco_num - 1) + aco_min
  mean_g = fltarr(aco_num)
  sigm_g = fltarr(aco_num)
;  mean_T = fltarr(aco_num)
;  sigm_T = fltarr(aco_num)

  ;make a dgr map
  for i = 0, aco_num - 1 do begin
    dgr = md / (mhi + ico * aco(i))
    mean_g(i) = mean(dgr,/nan)
    sigm_g(i) = sqrt(variance(dgr,/nan))
;    mean_t(i) = biweight_mean(dgr(mask),sig)
;    sigm_t(i) = sig
  endfor

  aco_bst=aco(where(sigm_g eq min(sigm_g)))
  dgr = md / (mhi + ico*aco_bst(0))

  dgr_output, dgr, mhi, aco_bst(0)*ico, aco_bst, sigm_g, regi

  return, aco

end

;************************************************************

pro dgr_output, dgr, hi, h2, aco, sig, reg


  sp = strsplit(reg, '/', /extract)
  sp = strsplit(sp(1), '.', /extract)
  filename=sp(0) + 'aco_output.ps'

  !p.multi=[0,2,1]
  set_plot, 'ps'
  device, filename=filename, /inches, xsize=10, ysize=6

    ;show scatter for sig_hi, sig_h2 vs log dgr
    cgplot, alog10(hi / h2), alog10(dgr), psym=5, ytitle='Log(DGR)', xtitle='Log(!4R!3!IH!I2!N!N / !4R!3!IHI!N)'
    cgplot, alog10((findgen(100)+1)/10), fltarr(100)+alog10(mean(dgr,/nan)), color='red', /overplot
    cgplot, alog10(aco), alog10(sig), xtitle='!4a!3!ICO!N', ytitle='!4r!3!IDGR!N'

  device, /close

end

;***********************************************************

pro conReg_main, mdp, mhip, icop, regip

;read in images
sd_d = mrdfits(mdp, 0, hdr_d)
sd_hi = mrdfits(mhip, 0, hdr_hi)
sd_co = mrdfits(icop, 0, hdr_co)

;mask out unuseable values
msk=where(finite(sd_d) ne 1, nel)
sd_hi(msk) = -1*!values.f_nan
sd_co(msk) = -1*!values.f_nan

;might want to make a blank file for the region stuff, or just pass in an image

;create an array for the regions
spawn, 'ls '+regip+'/*.dat', regis

;cycle through regions 
for i = 0, n_elements(regis)-1 do begin

  reg_d = regAn(regis(i), sd_d)
  reg_hi = regAn(regis(i), sd_hi)
  reg_co = regAn(regis(i), sd_co)

  aco = aco_calc(reg_d, reg_hi, reg_co, regis(i))

  ;stop

  ;dgr_output, 

endfor

end
