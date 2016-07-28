pro result2d, mode=mode, $
  thick=thick, length=length, xrange=xrange, yrange=yrange, nx=nx, ny=ny, title=title 

  !p.font = 0 
  
  ii = dcomplex(0d0, 1d0)

  if not keyword_set(title) then begin
;     title = '$\rho^{!X1/2}$'+'!X'
;     title+= '(' + '$\delta$' + '!X' + 'v' + '$_x$'+'!X, ' + '$\delta$' +'!Xv'+'$_z$'+'!X), '+'$\delta\epsilon$' + '!X'

     title = '$\rho^{!X1/2}_{tot}$'+'!X'
     title+= '(' + '$\delta$' + '!X' + 'v' + '$_x$'+'!X, ' + '$\delta$' +'!Xv'+'$_z$'+'!X), '+'$\delta$!X($\rho_{!Xdust}$!X/$\rho_{!Xtot}$!X)' +'!X'
     title = tex2idl(title)
  endif
  if not keyword_set(nx) then nx = 24
  if not keyword_set(ny) then ny = 24 
  if not keyword_set(thick) then thick=2
  if not keyword_set(length) then length=2   


  params = dblarr(4,1)
  openr,1,'params.dat'
  readf,1,params
  close,1
  nz      = fix(params(0,0))
  smallhg = params(1,0)
  kx      = params(2,0)
  Hd      = params(3,0)

  nz    = file_lines('basic.dat')
  basic = dblarr(7,nz)
  openr,1,'basic.dat'
  readf,1,basic
  close,1
  
  zaxis_data  = basic(0,*)/Hd
  rho    = exp(basic(1,*))
  eps    = basic(2,*)
  tstop  = basic(3,*)	
  omega2 = basic(4,*)
  kappa2 = basic(5,*)
  vshear = basic(6,*)
  
  if not keyword_set(yrange) then yrange = [-1d0,1d0]*max(zaxis_data)
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if not keyword_set(xrange) then xrange=[-1d0,1d0]*!dpi/kx  

  xaxis = xrange(0) + (xrange(1) - xrange(0))*dindgen(nx)/(nx-1d0)
  xaxis_big = xrange(0) + (xrange(1) - xrange(0))*dindgen(nz)/(nz-1d0)
  
  zaxis = yrange(0) + (yrange(1) - yrange(0))*dindgen(ny)/(ny-1d0)
  zaxis_big = yrange(0) + (yrange(1) - yrange(0))*dindgen(nz)/(nz-1d0)

  nmodes       = file_lines('eigenvalues.dat')
  eigenvalues  = dblarr(2,nmodes)
  openr,1,'eigenvalues.dat'
  readf,1,eigenvalues
  close,1
  
  growth = eigenvalues(0,*)/smallhg
  freq   = eigenvalues(1,*)/smallhg

  nlines = nmodes*nz + 0L  
  eigenvectors  = dblarr(10,nlines)
  openr,1,'eigenvectors.dat'
  readf,1,eigenvectors 
  close,1



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
  
  print, 'mode no., growth, freq', mode, growth(ngrid), freq(ngrid)

  nbeg = ngrid*nz 
  nend = nbeg + nz - 1 

 
  eps    = dcomplex(eigenvectors(2, nbeg:nend), eigenvectors(3, nbeg:nend))
  vx     = dcomplex(eigenvectors(4, nbeg:nend), eigenvectors(5, nbeg:nend))
  vz     = dcomplex(eigenvectors(8, nbeg:nend), eigenvectors(9, nbeg:nend))

;  eps *=dcomplex(0d0,-1d0) ;eps(nz-1)
;  vx  *=dcomplex(0d0,-1d0) ;eps(nz-1)
;  vz  *=dcomplex(0d0,-1d0) ;eps(nz-1)


  eps_re = real_part(eps)
  eps_im = imaginary(eps)
  
  vx_re  = real_part(vx)
  vx_im  = imaginary(vx)

  vz_re  = real_part(vz)
  vz_im  = imaginary(vz)

  re_eps = spline(zaxis_data,eps_re,zaxis_big)
  im_eps = spline(zaxis_data,eps_im,zaxis_big)

  re_vx = spline(zaxis_data,vx_re,zaxis)
  im_vx = spline(zaxis_data,vx_im,zaxis)
 
  re_vz = spline(zaxis_data,vz_re,zaxis)
  im_vz = spline(zaxis_data,vz_im,zaxis)
 
  rho1 =  spline(zaxis_data, rho, zaxis)

  sqrtrho = sqrt(rho1)
  
  vx2d   = dblarr(nx,ny)
  vz2d   = dblarr(nx,ny)
  eps2d  = dblarr(nz,nz)
  
  for i=0, nx-1 do begin
     for k=0, ny-1 do begin
        vx2d(i,k) =sqrtrho(k)*(re_vx(k)*cos(kx*xaxis(i)) - im_vx(k)*sin(kx*xaxis(i)))
        vz2d(i,k) =sqrtrho(k)*(re_vz(k)*cos(kx*xaxis(i)) - im_vz(k)*sin(kx*xaxis(i)))
     endfor
  endfor

  for i=0, nz-1 do begin
     for k=0, nz-1 do begin
        eps2d(i,k) = (re_eps(k)*cos(kx*xaxis_big(i)) - im_eps(k)*sin(kx*xaxis_big(i)))
     endfor
  endfor
  
  eps2d /= max(abs(eps2d))
  

;  xtitle = '!X'+ tex2idl('x/H'+'$_d$') + '!X'
;  ytitle = '!X'+ tex2idl('z/H'+'$_d$') + '!X'

  
  xtitle = '!X'+ tex2idl('x/H'+'$_{dust}$') + '!X'
  ytitle = '!X'+ tex2idl('z/H'+'$_{dust}$') + '!X'

  loadct,3,/silent
  set_plot, 'ps'
  device, filename='result2d.ps'$
          ,/color, bits_per_pixel=8,xsize=16, ysize=16
  
  levels = -1d0 + 2d0*dindgen(24)/23
  contour, eps2d, xaxis_big, zaxis_big,/fill,levels=levels,title=textoidl(title) $
           ,ymargin=[4,2],xmargin=[6,8], charsize=2, xtickinterval=xtickinterval, ytickinterval=xtickinterval, $
           xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, xtitle=xtitle, ytitle=ytitle
  colorbar, position=[0.8, 0.177, 0.85, 0.913],/vertical,/right,range=[-1,1],format='(f5.2)',CHARSIZE=2, text_color='k'
  
  loadct, 34,/silent
  velovect, vx2d, vz2d, xaxis, zaxis,/overplot,length=length,thick=thick, color=128
  
  device,/close
  
 
end

