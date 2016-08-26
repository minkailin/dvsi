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
  
  growth      = eigenvalues(0,*)/smallhg
  freq        = eigenvalues(1,*)/smallhg
  
  eigen_out = dcomplexarr(nmodes)
  eigen_out = dcomplex(growth, freq)
 
  return, eigen_out 
end


function resolution_test, eigen_lowres, eigen_highres

  common tolerance, err

  tol = err 

  nvals = n_elements(eigen_lowres)
  
  cnt = 0 
  for i=0, nvals-1 do begin
     sigma = eigen_lowres(i)
     temp  = min(abs(sigma - eigen_highres),grid)
     sigmahr = eigen_highres(grid)
     diff = abs( (sigma - sigmahr)/sigmahr )


     if diff gt tol then begin
        eigen_lowres(i) = 0d0 
     endif else begin
        cnt += 1
     endelse
  endfor
  
  eigen_out = eigen_lowres(where(eigen_lowres ne 0d0))
  
  return, eigen_out 
end


pro compare_eigenvals, loc=loc, highres=highres, xrange=xrange, yrange=yrange, label=label, legend=legend, title=title, ct=ct, tol=tol, mode=mode

  common tolerance, err 

  !p.font = 0 

  ncases = n_elements(loc)
  
  location=strcompress(loc,/remove_all)
  
  eigenvalues1 = get_data(loc[0])

  if not keyword_set(tol) then begin
   err = 1d-6 
  endif else begin
   err = tol 
  endelse 


  if keyword_set(highres) then begin
     location2=strcompress(highres,/remove_all)
     eigenvalues2 = get_data(highres[0])

     eigen_out = resolution_test(eigenvalues1, eigenvalues2)
     
  endif else begin
     eigen_out = eigenvalues1
  endelse


  growth      = real_part(eigen_out(*))
  freq        = imaginary(eigen_out(*))

  
  ytitle  = tex2idl('$s/(h$'+'!X'+'$_g$'+'!X'+'$\Omega$'+'$_K)$') + '!X'
  xtitle  = tex2idl('$\omega$'+'!X'+'$/(h$'+'!X'+'$_g$'+'!X'+'$\Omega$'+'$_K)$') + '!X'

 
  col_arr = 256d0*dindgen(ncases)/ncases

  set_plot, 'ps'

  if not keyword_set(title) then title = '' 
  
  loadct,0,/silent
  device, filename='compare_eigenvals.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot,  freq, growth,xmargin=[8,1.5],ymargin=[3.25,1.75], ystyle=1   $
         ,charsize=2, thick=4, psym=2, symsize=2, xrange=xrange, yrange=yrange $
         , xtitle=xtitle, ytitle=ytitle, title = tex2idl(title+'!X') $
         ,ytickinterval=ytickinterval  $
         ,xtickinterval=xtickinterval, xstyle=1, color = col_arr[0],/nodata





  if not keyword_set(mode) then begin
     temp = max(growth, ngrid)
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

   if not keyword_set(ct) then ct = 11
  loadct,ct,/silent
 


 
  for n=0, ncases-1 do begin
     
     eigenvalues1 = get_data(loc[n])

     
     if keyword_set(highres) then begin
        location2=strcompress(highres,/remove_all)
        eigenvalues2 = get_data(highres[n])
        
        eigen_out = resolution_test(eigenvalues1, eigenvalues2)
        
     endif else begin
        eigen_out = eigenvalues1
     endelse
     
    
     growth      = real_part(eigen_out(*))
     freq        = imaginary(eigen_out(*))


     
     oplot, freq, growth, psym=4, symsize=1.5, thick=4, color=col_arr[n]
  endfor

   color_arr = dindgen(2)*256d0/2.
;  oplot, [1,1]*freq[ngrid], [1,1]*growth[ngrid], psym=2,symsize=1.5,color=color_arr(1)







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

