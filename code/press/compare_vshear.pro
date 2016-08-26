FUNCTION logticks_exp, axis, index, value
   ; Determine the base-10 exponent
   exponent   = LONG( ALOG10( value ) )
   ; Construct the tickmark string based on the exponent
;   tickmark = '10!E' + STRTRIM( STRING( exponent ), 2 ) + '!N'
    tickmark = textoidl('10^{' + STRTRIM( STRING( exponent ), 2 )+'}')
   ; Return the formatted tickmark string
   RETURN, tickmark
END

function get_data, loc 

  location=strcompress(loc,/remove_all)

  params = dblarr(5,1)
  openr,1, filepath('params.dat',root_dir='.',subdir=location)
  readf,1,params
  close,1
  
   params = dblarr(5,1)
   openr,1,filepath('params.dat',root_dir='.',subdir=location)
   readf,1,params
   close,1
   nz      = fix(params(0,0))
   
   smallhg = params(1,0)
   kx      = params(2,0)
   Hd      = params(3,0)
   smallq  = params(4,0)
   
   nz    = file_lines(filepath('basic.dat',root_dir='.',subdir=location))
   basic = dblarr(7,nz)
   openr,1,filepath('basic.dat',root_dir='.',subdir=location)
   readf,1,basic
   close,1
   
   zaxis  = basic(0,*)
   lnrho = basic(1,*)
   eps    = basic(2,*)
   tstop  = basic(3,*)
   omega2 = basic(4,*)
   kappa2 = basic(5,*)
   vshear = basic(6,*)
 
   dgratio = eps/(1d0 - eps)
   lnrhog  = lnrho + alog(1d0 - eps)  
  
   dlnrhog = deriv(zaxis, lnrhog)
   ddgratio= deriv(zaxis, dgratio)
   Nzsq    = dlnrhog*ddgratio/(1d0+dgratio)^2 

   output = dblarr(4, nz) 
   output(0,*) = zaxis
   output(1,*) = sqrt(omega2)
   output(2,*) = vshear/smallhg 
   output(3,*) = Nzsq/smallhg 

   return, output 
end


pro compare_vshear, loc=loc, xrange=xrange, yrange=yrange, label=label, legend=legend, title=title, ct=ct

  !p.font = 0 

  ncases = n_elements(loc)
  
  location=strcompress(loc,/remove_all)

  output = get_data(loc[0])
  zaxis  = output(0,*)
  vshear = output(2,*) ;this is scaled by h_g

  ytitle  = tex2idl('$(r\partial$'+'!X'+'$_z$'+'$\Omega^2$!X)/$h_{!Xg}$$\Omega_K^2$') + '!X'
  xtitle  = tex2idl('$z/H_g$') + '!X'
  
  col_arr = 256d0*dindgen(ncases)/ncases

  set_plot, 'ps'

  if not keyword_set(title) then title = '' 
  
  loadct,0,/silent
  device, filename='compare_vshear.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot,  zaxis, vshear,xmargin=[8,1.5],ymargin=[3.25,1.75], ystyle=1   $
         ,charsize=2, thick=4, psym=2, symsize=2, xrange=xrange, yrange=yrange $
         , xtitle=xtitle, ytitle=ytitle, title = tex2idl(title+'!X') $
         ,ytickinterval=ytickinterval  $
         ,xtickinterval=xtickinterval, xstyle=1, color = col_arr[0],/nodata

  if not keyword_set(ct) then ct = 11
  loadct,ct,/silent
  for n=0, ncases-1 do begin
     
    output = get_data(loc[n])
  zaxis  = output(0,*)
  vshear = output(2,*) ;this is scaled by h_g 
     
     oplot, zaxis, vshear, thick=4, color=col_arr[n] 
  endfor

  if keyword_set(legend) then begin
     x0=legend(0)
     x1=legend(1)
     y0=legend(2)
     dy=legend(3)
     for j=0, ncases-1 do begin
;        ynew = y0 - dy*j
;        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, color=color_arr(j)
        xyouts, x0, y0-dy*j,tex2idl(label(j)),charsize=2,color=col_arr(j)
     endfor
  endif


  device,/close 

  output = get_data(loc[0])
  zaxis  = output(0,*)
  omega  = output(1,*)
  vshear = output(2,*) ;this is scaled by h_g
  Nzsq   = output(3,*) ;this is scaled by h_g

  vshear /= 2d0*omega 


  ytitle  = tex2idl('$|r\partial$'+'!X'+'$_z$'+'$\Omega$!X|/$h_{!Xg}$$\Omega_K$!X (solid) !C!C $N^2_{!Xz}$/$h_{!Xg}$$\Omega_K^2$!X (dotted)') + '!X'
  xtitle  = tex2idl('$z/H_g$') + '!X' 

  set_plot, 'ps'

  if not keyword_set(title) then title = ''

  loadct,0,/silent
  device, filename='compare_vshear_Nz2.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot,  zaxis, abs(vshear),xmargin=[8,1.5],ymargin=[3.25,1.75], ystyle=1   $
         ,charsize=2, thick=4, psym=2, symsize=2, xrange=xrange, yrange=yrange $
         , xtitle=xtitle, ytitle=ytitle, title = tex2idl(title+'!X') $
         ,ytickinterval=1  $
         ,xtickinterval=xtickinterval, xstyle=1, color = col_arr[0],/nodata
;  oplot, [-1,1]*1d10, [0,0], thick=1
  if not keyword_set(ct) then ct = 11
  loadct,ct,/silent
  for n=0, ncases-1 do begin

  output = get_data(loc[n])
  zaxis  = output(0,*)
  omega  = output(1,*)
  vshear = output(2,*) 
  Nzsq   = output(3,*)

  vshear /= 2d0*omega

  oplot, zaxis, abs(vshear), thick=4, color=col_arr[n]
  oplot, zaxis, Nzsq, thick=4, color=col_arr[n], linestyle=1 
  endfor

  if keyword_set(legend) then begin
     x0=legend(0)
     x1=legend(1)
     y0=legend(2)
     dy=legend(3)
     for j=0, ncases-1 do begin
;        ynew = y0 - dy*j
;        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, color=color_arr(j)
        xyouts, x0, y0-dy*j,tex2idl(label(j)),charsize=2,color=col_arr(j)
     endfor
  endif


  device,/close
  
end

