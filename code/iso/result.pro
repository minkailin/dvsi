FUNCTION logticks_exp, axis, index, value
   ; Determine the base-10 exponent
   exponent   = LONG( ALOG10( value ) )
   ; Construct the tickmark string based on the exponent
;   tickmark = '10!E' + STRTRIM( STRING( exponent ), 2 ) + '!N'
    tickmark = textoidl('10^{' + STRTRIM( STRING( exponent ), 2 )+'}')
   ; Return the formatted tickmark string
   RETURN, tickmark
END

pro result

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


end

