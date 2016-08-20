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
  
  nz      = fix(params(0,0))
  smallhg = params(1,0)
  kx      = params(2,0)
  Hd      = params(3,0)
  smallq  = params(4,0) 

  nmodes       = file_lines(filepath('eigenvalues.dat',root_dir='.',subdir=location))
  eigenvalues  = dblarr(2,nmodes)
  openr,1, filepath('eigenvalues.dat',root_dir='.',subdir=location)
  readf,1,eigenvalues
  close,1



   array = dblarr(3, nmodes)
  openr,1,filepath('Hd.dat',root_dir='.',subdir=location)
  readf,1, array
  close,1
  Hdust       = array(0,*)
  dgratio     = array(1,*)
  kaxis       = array(2,*)

  growth = eigenvalues(0,*)/smallhg
  freq   = eigenvalues(1,*)/smallhg

  output = dblarr(5, nmodes)
  output(0,*) = dgratio 
  output(1,*) = Hdust 
  output(2,*) = kaxis
  output(3,*) = growth 
  output(4,*) = freq 

  return, output
end


pro compare_eigenvals, loc=loc, xrange=xrange, yrange=yrange, label=label, legend=legend, title=title, ct=ct

  !p.font = 0 

  ncases = n_elements(loc)
  
  location=strcompress(loc,/remove_all)

  output       = get_data(loc[0])
  Hdust        = output(1,*)
  growth       = output(3,*)


  ytitle  = tex2idl('$s_{max}/(h$'+'!X'+'$_g$'+'!X'+'$\Omega$'+'$_K)$') + '!X'
;  xtitle  = tex2idl('$\rho$!X$_{!Xd0}$!X$/\rho$!X$_{!Xg0}$')+'!X'
  xtitle  = tex2idl('$H$!X$_{!Xd}$!X$/H$!X$_{!Xg}$')+'!X'
 
  col_arr = 256d0*dindgen(ncases)/ncases

  set_plot, 'ps'

  if not keyword_set(title) then title = '' 
  
  loadct,0,/silent
  device, filename='compare_eigenvals.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot,  Hdust, growth,xmargin=[8,1.5],ymargin=[3.25,1.75], ystyle=1   $
         ,charsize=2, thick=4, psym=2, symsize=2, xrange=xrange, yrange=yrange $
         , xtitle=xtitle, ytitle=ytitle, title = tex2idl(title+'!X') $
         ,ytickinterval=ytickinterval, xtickformat='logticks_exp'  $
         ,xtickinterval=xtickinterval, xstyle=1, color = col_arr[0],/nodata,/xlog 

  if not keyword_set(ct) then ct = 11
  loadct,ct,/silent
  for n=0, ncases-1 do begin
     
     output       = get_data(loc[n])
     Hdust        = output(1,*)
  growth       = output(3,*)
  

     oplot, Hdust, growth, thick=4, color=col_arr[n];,/xlog,xtickformat='logticks_exp'
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

