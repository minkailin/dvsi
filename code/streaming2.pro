FUNCTION logticks_exp, axis, index, value
   ; Determine the base-10 exponent
   exponent   = LONG( ALOG10( value ) )
   ; Construct the tickmark string based on the exponent
;   tickmark = '10!E' + STRTRIM( STRING( exponent ), 2 ) + '!N'
    tickmark = textoidl('10^{' + STRTRIM( STRING( exponent ), 2 )+'}')
   ; Return the formatted tickmark string
   RETURN, tickmark
END

pro streaming2, eta=eta, kx=kx, kz=kz, dgratio=dgratio, tstop=tstop, title=title

!p.font=0

nmodes=5

if not keyword_set(title) then title = ''

ii    = dcomplex(0d0, 1d0)
  dfrac = dgratio/(1d0 + dgratio)

  nkz      =64
  kz_axis = 10d0^(alog10(kz[0]) + alog10(kz[1]/kz[0])*dindgen(nkz)/(nkz-1d0)) ;kz[0]+ (kz[1]-kz[0])*dindgen(nkz)/(nkz-1d0) 

  rate = dblarr(nmodes,nkz) 
  freq =  dblarr(nmodes,nkz) 

  for n = 0, nkz-1 do begin 

   kz = kz_axis(n)

  ksq   = kx^2d0 + kz^2d0
  smallh = 0.05
  kappa2 = 1d0 - 6d0*(1d0-dfrac)*eta*smallh

  kzsq  = kz^2d0/ksq 

  c0 = 2d0*kx*tstop*(1d0 - dfrac)*(1d0 - 2d0*dfrac)*kzsq*kappa2 
  
  c1 = kzsq*kappa2

  c2 = -ii*tstop*( dfrac*kappa2 - 2d0*ii*kx*(1d0-dfrac)^2d0 )

  c3 = eta^2d0*(kappa2 - 2d0*ii*kx*(1d0 - dfrac))/(ksq*(1d0-dfrac)) + 1d0;*( 1d0 + 2d0*ii*eta^2d0/kx )
  c3 =-c3
  
  c4 = ii*dfrac*tstop

  c5 = eta^2d0/(1d0-dfrac)/ksq 

  coeffs = [c0, c1, c2, c3, c4, c5] 

;  coeffs=[c0,c1,c2,c3,c4]

  roots = fz_roots(coeffs, /double, eps=1d-15)
  
  rfreq  =-real_part(roots)
  growth = imaginary(roots)

  for j=0, nmodes-1 do begin
     rate[j,n] = max(growth[j]) 
     freq[j,n] = max(rfreq[j])
  endfor
  
endfor

filter = where(rate le 0d0)
rate(filter) = -1d0

loadct,10,/silent
color_arr = dindgen(5)*256d0/5.


xtitle  =tex2idl('$K_{!Xz}$') +'!X'
ytitle = tex2idl('$s/\Omega_{!XK}$') + '!X'

set_plot,'ps'
file = strcompress('streaming2.ps',/remove_all)
device, filename=file $
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color 
plot, kz_axis, rate[0,*],xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
      ,charsize=1.5, thick=4, xrange=xrange, title=tex2idl(title), xtitle=xtitle,$
      linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2, yrange=[0d0,max(rate)], /xlog, xtickformat='logticks_exp', $
      psym=2,symsize=1
       
for n=1, nmodes-1 do begin
   oplot, kz_axis, rate[n,*], thick=4, color=color_arr[0], psym=2,symsize=1
endfor 

xyouts, 1500, 0.7, 'streaming instability', charsize=1.5,align=1
xyouts, 0.01, 0.1, 'overstable epicycle', charsize=1.5

if keyword_set(legend) then begin
   x0=legend(0)
   x1=legend(1)
   y0=legend(2)
   dy=legend(3)
   for j=0, n_elements(label)-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
     endfor
endif
device,/close

ytitle = tex2idl('$\omega/\Omega_{!XK}$') + '!X'

set_plot,'ps'
file = strcompress('streaming2_re.ps',/remove_all)
device, filename=file $
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color 
plot, kz_axis, freq[0,*],xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
      ,charsize=1.5, thick=4, xrange=xrange, title=title, xtitle=xtitle,$
      linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2, yrange=[min(freq),max(freq)], /xlog, xtickformat='logticks_exp', $
      psym=2,symsize=1
       
for n=1, nmodes-1 do begin
oplot, kz_axis, freq[n,*], thick=4, color=color_arr[0], psym=2,symsize=1
endfor 

if keyword_set(legend) then begin
   x0=legend(0)
   x1=legend(1)
   y0=legend(2)
   dy=legend(3)
   for j=0, n_elements(label)-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
     endfor
endif
device,/close

end
