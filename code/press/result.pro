FUNCTION logticks_exp, axis, index, value
   ; Determine the base-10 exponent
   exponent   = LONG( ALOG10( value ) )
   ; Construct the tickmark string based on the exponent
;   tickmark = '10!E' + STRTRIM( STRING( exponent ), 2 ) + '!N'
    tickmark = textoidl('10^{' + STRTRIM( STRING( exponent ), 2 )+'}')
   ; Return the formatted tickmark string
   RETURN, tickmark
END

pro result, loc=loc, xrange=xrange, yrange=yrange, mode=mode 

  !p.font = 0 

  if not keyword_set(loc) then begin 
     loc = './'
  endif else begin
     location=strcompress(loc,/remove_all)
  endelse

  params = dblarr(5,1)
  openr,1,filepath('params.dat',root_dir='.',subdir=location)
  readf,1,params
  close,1
  nz      = fix(params(0,0))
  smallhg = params(1,0)
  kx      = params(2,0)
  Hd      = params(3,0)
  smallq  = params(4,0) 


  fname = 'basic.dat'
  nz    = file_lines(filepath(fname,root_dir='.',subdir=location))
  basic = dblarr(7,nz)
  openr,1, filepath(fname,root_dir='.',subdir=location)
  readf,1,basic
  close,1

  zaxis  = basic(0,*) 
  lnrho = basic(1,*)
  eps    = basic(2,*)
  tstop  = basic(3,*)	
  omega2 = basic(4,*)
  kappa2 = basic(5,*)
  vshear = basic(6,*)

  lnrhog  = lnrho + alog(1d0 - eps) 
  dgratio = eps/(1d0 - eps) 

  xtitle  = tex2idl('$z/H_g$') + '!X'
  ytitle  = tex2idl('$10^3(\Omega^2/\Omega_K^2 - 1)$') + '!X'
  set_plot, 'ps'
  fname = 'omega2.ps'
  device, filename=filepath(fname,root_dir='.',subdir=location) $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, (omega2-1d0)*1d3,xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=2, thick=4, xrange=xrange, yrange=yrange $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $ 
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close	

   ytitle  = tex2idl('$10^3(\kappa^2/\Omega^2 - 1)$') + '!X'
   fname = 'kappa2.ps'
  set_plot, 'ps'
  device, filename=filepath(fname,root_dir='.',subdir=location) $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, (kappa2/omega2-1d0)*1d3,xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=2, thick=4, xrange=xrange, yrange=yrange $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close

  ytitle  = tex2idl('$t_{!Xstop}\Omega_{!XK}$') + '!X'
   fname = 'tstop.ps'
  set_plot, 'ps'
  device, filename=filepath(fname,root_dir='.',subdir=location) $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, tstop,xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=2, thick=4, xrange=xrange, yrange=yrange $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close













 
  ytitle  = tex2idl('$(r\partial$'+'!X'+'$_z$'+'$\Omega^2$!X)/$h_{!Xg}$$\Omega_K^2$') + '!X'
  fname = 'vshear.ps'
  set_plot, 'ps'
  device, filename=filepath(fname,root_dir='.',subdir=location)  $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, vshear/smallhg,xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=2, thick=4, xrange=xrange, yrange=yrange $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $ 
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close 

  ytitle  = tex2idl('$ln(\rho$'+'!X'+'$_g$'+'$/\rho$'+'!X'+'$_{g0})$') + '!X'
  fname = 'rhog.ps'
  set_plot, 'ps'
  device, filename=filepath(fname,root_dir='.',subdir=location) $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, lnrhog,xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=2, thick=4, xrange=xrange, yrange=yrange $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  oplot, zaxis, -zaxis^2/2d0, thick=8, linestyle=1
  device,/close	

  ytitle  = tex2idl('$\rho$'+'!X'+'$_d$'+'$/\rho$'+'!X'+'$_g$') + '!X'
  set_plot, 'ps'
  fname = 'dgratio.ps'
  device, filename= filepath(fname,root_dir='.',subdir=location)$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, dgratio,xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=2, thick=4, xrange=xrange, yrange=yrange $
        , xtitle=xtitle, ytitle=ytitle,/ylog, ytickformat='logticks_exp' $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close

  ytitle  = tex2idl('$\rho_{!Xd}/\rho_{!Xd0}$') + '!X'
  set_plot, 'ps'
  fname = 'rhod.ps'
  device, filename= filepath(fname,root_dir='.',subdir=location)$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, dgratio*exp(lnrhog),xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=2, thick=4, xrange=xrange, yrange=yrange $
        , xtitle=xtitle, ytitle=ytitle $;,/ylog, ytickformat='logticks_exp' $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close












;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  fname = 'eigenvalues.dat'
  nmodes       = file_lines(filepath(fname,root_dir='.',subdir=location))
  eigenvalues  = dblarr(2,nmodes)
  openr,1,filepath(fname,root_dir='.',subdir=location)
  readf,1,eigenvalues
  close,1
  
  growth = eigenvalues(0,*)/smallhg
  freq   = eigenvalues(1,*)/smallhg

  ytitle  = tex2idl('$s/(h$'+'!X'+'$_g$'+'!X'+'$\Omega$'+'$_K)$') + '!X'
  xtitle  = tex2idl('$\omega$'+'!X'+'$/(h$'+'!X'+'$_g$'+'!X'+'$\Omega$'+'$_K)$') + '!X'


  loadct,5,/silent
  fname = 'eigenvalues.ps'
  set_plot, 'ps'
  device, filename=filepath(fname,root_dir='.',subdir=location) $
         ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
;          ,/color, bits_per_pixel=8,xsize=28.44444, ysize=16
  plot,  freq, growth,xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=1  $
        ,charsize=2, thick=4, psym=2, symsize=2, xrange=xrange, yrange=yrange $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1

;  color_arr = dindgen(2)*256d0/2.
;  oplot, [1,1]*freq[0], [1,1]*growth[0], psym=2,symsize=1.5,color=color_arr(1)
;  device,/close

  nlines = nmodes*nz + 0L  
  eigenvectors  = dblarr(10,nlines)
  fname = 'eigenvectors.dat'
  openr,1,filepath(fname,root_dir='.',subdir=location)
  readf,1,eigenvectors 
  close,1

 

  if not keyword_set(mode) then begin
  temp = max(growth, ngrid) 
;  temp = min(abs(freq),ngrid)
  mode = ngrid + 1 
  endif else begin
  if(n_elements(mode eq 1))then begin
  ngrid = mode - 1 
  endif else begin
  rate_target = mode[0]
  freq_target = mode[1]
  eigen_target = dcomplex(rate_target, freq_target)
  temp = min( abs(dcomplex(growth,freq)-eigen_target), ngrid )
  mode = ngrid + 1 
  endelse
  endelse 
  

;  oplot, [1,1]*freq[ngrid], [1,1]*growth[ngrid], psym=2,symsize=2,color=color_arr(1),thick=4 

  device,/close 

  print, 'mode no., growth, freq', mode, growth(ngrid), freq(ngrid)

  nbeg = ngrid*nz 
  nend = nbeg + nz - 1 

  bigW = dcomplex(eigenvectors(0, nbeg:nend), eigenvectors(1, nbeg:nend))
  dfrac= dcomplex(eigenvectors(2, nbeg:nend), eigenvectors(3, nbeg:nend))
  vx   = dcomplex(eigenvectors(4, nbeg:nend), eigenvectors(5, nbeg:nend))
  vy   = dcomplex(eigenvectors(6, nbeg:nend), eigenvectors(7, nbeg:nend))
  vz   = dcomplex(eigenvectors(8, nbeg:nend), eigenvectors(9, nbeg:nend))

  vrad = vx 
  vmag = sqrt(abs(vx)^2.+abs(vy)^2.+abs(vz)^2.)
  vh   = sqrt(abs(vx)^2.+abs(vy)^2.) 
  vmeri2 = abs(vz)^2. + abs(vx)^2.

  ;; bigW  *= conj(bigW(nz-1))
  ;; dfrac *= conj(dfrac(nz-1))
  ;; vx    *= conj(vx(nz-1)) 
  ;; vy    *= conj(vy(nz-1))
  ;; vz    *= conj(vz(nz-1))

  bigW   /= max(abs(bigW))
  dfrac  /= max(abs(dfrac))
  vx     /= max(abs(vx))
  vy     /= max(abs(vy))
  vz     /= max(abs(vz))
  

  xtitle  = tex2idl('$z/H_g$') + '!X'
  
  fname = 'eigenvec.ps'
  set_plot,'ps'
;  file = strcompress('eigenvec.ps',/remove_all)
  device, filename=filepath(fname,root_dir='.',subdir=location) $
          ,bits_per_pixel=8,xsize=8, ysize=9,xoffset=0,yoffset=0,/inches,/color
  multiplot,/reset
  multiplot,[1,3],margins=[0.15,0.1,0.02,0.01],rowspacing=0.05
  ytitle  = tex2idl('$\delta$'+'!X'+'v'+'$_z$') + '!X'
  plot, zaxis, real_part(vz),xmargin=[8.3,1.7],ymargin=[5,0], ystyle=1, xstyle=1 $
        ,charsize=2, thick=4,yrange=[-1,1], ytitle=ytitle, $
        linestyle = 0, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2
  oplot, zaxis, imaginary(vz), thick=4, linestyle=1 
  multiplot
  ytitle  = tex2idl('$\delta\epsilon$') + '!X'
  plot, zaxis, real_part(dfrac),xmargin=[8.3,1.7], ystyle=2, xstyle=1 $
        ,charsize=2, thick=4,yrange=[-1,1], ytitle=ytitle $
        ,linestyle = 0, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2
  oplot, zaxis, imaginary(dfrac), thick=4, linestyle=1
  multiplot
  ytitle  = tex2idl('$\delta\rho$'+'!X'+'/'+'$\rho$') + '!X'
  plot, zaxis,real_part(bigW),xmargin=[8.3,1.7], xstyle=1, ystyle=1 $
        ,charsize=2, thick=4,yrange=[-1,1], xtitle=xtitle $
        ,linestyle = 0, ytitle = ytitle, xtickinterval=xtickinterval,charthick=2
  oplot, zaxis, imaginary(bigW), thick=4, linestyle=1
  multiplot,/reset
  device,/close
  
  




  







  ytitle  = tex2idl('(|'+'$\delta$'+'!X'+'v'+'$_x$'+'!X'+'|'+'$^2$'+'+|'+'$\delta$'+'!X'+'v'+'$_y$'+'!X'+'|'+'$^2$'+')'+'$^{1/2}$'+'!X/|'+'$\delta$'+'!Xv'+'!X|') + '!X'
  fname = 'eigenvec_vh.ps'
  set_plot, 'ps'
  device, filename=filepath(fname,root_dir='.',subdir=location) $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis,  vh/vmag,xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=1  $
        ,charsize=2, thick=4, xrange=xrange $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close
































  err = dblarr(5, nlines)
  fname = 'error.dat'
  openr,1,filepath(fname,root_dir='.',subdir=location)
  readf,1,err
  close,1

  err1 = err(0,nbeg:nend) ;mass error
  err2 = err(1,nbeg:nend) ;vx error
  err3 = err(2,nbeg:nend) ;vy error
  err4 = err(3,nbeg:nend) ;vz error
  err5 = err(4,nbeg:nend) ;energy error 

  ytitle = 'normalized error'
  fname = 'eigenvec_err.ps'
  set_plot, 'ps'
  device, filename=filepath(fname,root_dir='.',subdir=location) $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis(1:nz-2),  err1(1:nz-2),xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=1  $
        ,charsize=2, thick=4, xrange=xrange $
        , xtitle=xtitle, ytitle=ytitle, yrange=[0.,max(err(*,nbeg+1:nend-1))] $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  oplot, zaxis(1:nz-2), err2(1:nz-2), thick=4, linestyle=1
  oplot, zaxis(1:nz-2), err3(1:nz-2), thick=4, linestyle=2
  oplot, zaxis(1:nz-2), err4(1:nz-2), thick=4, linestyle=3
  oplot, zaxis(1:nz-2), err5(1:nz-2), thick=4, linestyle=4
  device,/close
  
  nonadia = dblarr(4, nlines)
  fname = 'nonadia.dat'
  openr,1,filepath(fname,root_dir='.',subdir=location)
  readf,1,nonadia
  close,1

  divv = dcomplex(nonadia(0,nbeg:nend), nonadia(1,nbeg:nend))
  dC   = dcomplex(nonadia(2,nbeg:nend), nonadia(3,nbeg:nend))

  integrand1 =imaginary(-conj(divv)*dC*exp(lnrho))
  integrand2 =imaginary(-(1d0-eps)*exp(lnrho)*conj(divv)*vrad*smallhg*smallq)

  denom     =int_tabulated(zaxis, exp(lnrho)*vmeri2)

  result = int_tabulated(zaxis, integrand1 + integrand2)

  result/=-2d0*eigenvalues(1,ngrid)*denom

  print, 'growth from int relation', result/smallhg

  
  ytitle = 'non-adia. growth'
  fname = 'eigenvec_nonadia.ps'
  set_plot, 'ps'
  device, filename=filepath(fname,root_dir='.',subdir=location) $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis,  abs(integrand1),xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=1  $
        ,charsize=2, thick=4, xrange=xrange $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close

 
end

