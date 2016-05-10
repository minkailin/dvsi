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
        ,charsize=2, thick=4, psym=2, symsize=1.5, xrange=xrange, yrange=yrange $
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

  vmag = sqrt(abs(vx)^2.+abs(vy)^2.+abs(vz)^2.)
  vh   = sqrt(abs(vx)^2.+abs(vy)^2.) 
  vmeri2 = abs(vz)^2. + abs(vx)^2.

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

  set_plot,'ps'
  file = strcompress('eigenvec.ps',/remove_all)
  device, filename=file $
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
  set_plot, 'ps'
  device, filename='eigenvec_vh.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis,  vh/vmag,xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=1  $
        ,charsize=2, thick=4, xrange=xrange $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close
































  err = dblarr(5, nlines)
  openr,1,'error.dat'
  readf,1,err
  close,1

  err1 = err(0,nbeg:nend) ;mass error
  err2 = err(1,nbeg:nend) ;vx error
  err3 = err(2,nbeg:nend) ;vy error
  err4 = err(3,nbeg:nend) ;vz error
  err5 = err(4,nbeg:nend) ;energy error 

  ytitle = 'normalized error'
  set_plot, 'ps'
  device, filename='eigenvec_err.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis,  err1,xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=1  $
        ,charsize=2, thick=4, xrange=xrange $
        , xtitle=xtitle, ytitle=ytitle, yrange=[0.,max(err(*,nbeg:nend))] $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  oplot, zaxis, err2, thick=4, linestyle=1
  oplot, zaxis, err3, thick=4, linestyle=2
  oplot, zaxis, err4, thick=4, linestyle=3
  oplot, zaxis, err5, thick=4, linestyle=4
  device,/close
  
  nonadia = dblarr(4, nlines)
  openr,1,'nonadia.dat'
  readf,1,nonadia
  close,1

  divv = dcomplex(nonadia(0,nbeg:nend), nonadia(1,nbeg:nend))
  dC   = dcomplex(nonadia(2,nbeg:nend), nonadia(3,nbeg:nend))

  integrand =imaginary(-conj(divv)*dC*exp(lnrho))
  denom     =int_tabulated(zaxis, exp(lnrho)*vmeri2)

  result = int_tabulated(zaxis, integrand)
  result/=-2d0*eigenvalues(1,ngrid)*denom

  print, 'growth from int relation', result/smallhg

  
  ytitle = 'non-adia. growth'
  set_plot, 'ps'
  device, filename='eigenvec_nonadia.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis,  abs(integrand),xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=1  $
        ,charsize=2, thick=4, xrange=xrange $
        , xtitle=xtitle, ytitle=ytitle $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close

 
end

