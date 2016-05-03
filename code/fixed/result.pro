FUNCTION logticks_exp, axis, index, value
   ; Determine the base-10 exponent
   exponent   = LONG( ALOG10( value ) )
   ; Construct the tickmark string based on the exponent
;   tickmark = '10!E' + STRTRIM( STRING( exponent ), 2 ) + '!N'
    tickmark = textoidl('10^{' + STRTRIM( STRING( exponent ), 2 )+'}')
   ; Return the formatted tickmark string
   RETURN, tickmark
END

pro result, xrange=xrange, yrange=yrange, mode=mode 

  !p.font = 0 

  nz    = file_lines('basic.dat')
  basic = dblarr(7,nz)
  openr,1,'basic.dat'
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
  device, filename='omega2.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, (omega2-1d0)*1d3,xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=2, thick=4, xrange=xrange, yrange=yrange $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $ 
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close	

   ytitle  = tex2idl('$10^3(\kappa^2/\Omega^2 - 1)$') + '!X'
  set_plot, 'ps'
  device, filename='kappa2.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, (kappa2/omega2-1d0)*1d3,xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=2, thick=4, xrange=xrange, yrange=yrange $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close
 
  ytitle  = tex2idl('$10^3(r\partial$'+'!X'+'$_z$'+'$\Omega^2)/\Omega_K^2$') + '!X'
  set_plot, 'ps'
  device, filename='vshear.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, vshear*1d3,xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=2, thick=4, xrange=xrange, yrange=yrange $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $ 
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close 

  ytitle  = tex2idl('$ln(\rho$'+'!X'+'$_g$'+'$/\rho$'+'!X'+'$_{g0})$') + '!X'
  set_plot, 'ps'
  device, filename='rhog.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, lnrhog,xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=2, thick=4, xrange=xrange, yrange=yrange $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close	

  ytitle  = tex2idl('$\rho$'+'!X'+'$_d$'+'$/\rho$'+'!X'+'$_g$') + '!X'
  set_plot, 'ps'
  device, filename='dgratio.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, dgratio,xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=2, thick=4, xrange=xrange, yrange=yrange $
        , xtitle=xtitle, ytitle=ytitle,/ylog, ytickformat='logticks_exp' $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  params = dblarr(2,1)
  openr,1,'params.dat'
  readf,1,params
  close,1
  nz      = fix(params(0,0))
  smallhg = params(1,0)

  nmodes       = file_lines('eigenvalues.dat')
  eigenvalues  = dblarr(2,nmodes)
  openr,1,'eigenvalues.dat'
  readf,1,eigenvalues
  close,1
  
  growth = eigenvalues(0,*)/smallhg
  freq   = eigenvalues(1,*)/smallhg

  ytitle  = tex2idl('$s/(h$'+'!X'+'$_g$'+'!X'+'$\Omega$'+'$_K)$') + '!X'
  xtitle  = tex2idl('$\omega$'+'!X'+'$/(h$'+'!X'+'$_g$'+'!X'+'$\Omega$'+'$_K)$') + '!X'
  set_plot, 'ps'
  device, filename='eigenvalues.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot,  freq, growth,xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=2  $
        ,charsize=2, thick=4, psym=2, symsize=2, xrange=xrange, yrange=yrange $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=2
  device,/close

  nlines = nmodes*nz + 0L  
  eigenvectors  = dblarr(10,nlines)
  openr,1,'eigenvectors.dat'
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
  

  print, 'mode no., growth, freq', mode, growth(ngrid), freq(ngrid)

  nbeg = ngrid*nz 
  nend = nbeg + nz - 1 

  bigW = dcomplex(eigenvectors(0, nbeg:nend), eigenvectors(1, nbeg:nend))
  dfrac= dcomplex(eigenvectors(2, nbeg:nend), eigenvectors(3, nbeg:nend))
  vx   = dcomplex(eigenvectors(4, nbeg:nend), eigenvectors(5, nbeg:nend))
  vy   = dcomplex(eigenvectors(6, nbeg:nend), eigenvectors(7, nbeg:nend))
  vz   = dcomplex(eigenvectors(8, nbeg:nend), eigenvectors(9, nbeg:nend))

  bigW  *= conj(bigW(nz-1))
  dfrac *= conj(dfrac(nz-1))
  vx    *= conj(vx(nz-1)) 
  vy    *= conj(vy(nz-1))
  vz    *= conj(vz(nz-1))

 bigW   /= max(abs(bigW))
 dfrac  /= max(abs(dfrac))
 vx     /= max(abs(vx))
 vy     /= max(abs(vy))
 vz     /= max(abs(vz))


  xtitle  = tex2idl('$z/H_g$') + '!X'
  ytitle  = tex2idl('$\delta$'+'!X'+'v'+'$_z$') + '!X'
  set_plot, 'ps'
  device, filename='eigenvec_vz.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, real_part(vz),xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=2, thick=4, xrange=xrange, yrange=[-1,1] $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  oplot, zaxis, imaginary(vz), thick=4, linestyle=1 
  device,/close

  ytitle  = tex2idl('$\delta\epsilon$') + '!X'
  set_plot, 'ps'
  device, filename='eigenvec_eps.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, real_part(dfrac),xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=2, thick=4, xrange=xrange, yrange=[-1,1] $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  oplot, zaxis, imaginary(dfrac), thick=4, linestyle=1
  device,/close


end

