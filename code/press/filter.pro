pro filter, loc=loc, tol = tol  

  !p.font = 0 

  if not keyword_set(loc) then begin 
     loc = './'
  endif else begin
     location=strcompress(loc,/remove_all)
  endelse

  if not keyword_set(tol) then tol = 1d-6

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

  fname = 'eigenvalues.dat'
  nmodes       = file_lines(filepath(fname,root_dir='.',subdir=location))
  eigenvalues  = dblarr(2,nmodes)
  openr,1,filepath(fname,root_dir='.',subdir=location)
  readf,1,eigenvalues
  close,1
  
  print, 'nmodes=', nmodes 

  nlines = nmodes*nz + 0L  
  eigenvectors  = dblarr(10,nlines)
  fname = 'eigenvectors.dat'
  openr,1,filepath(fname,root_dir='.',subdir=location)
  readf,1,eigenvectors 
  close,1
  
  
  nonadia = dblarr(4, nlines)
  fname = 'nonadia.dat'
  openr,1,filepath(fname,root_dir='.',subdir=location)
  readf,1,nonadia
  close,1

  err = dblarr(5, nlines)
  fname = 'error.dat'
  openr,1,filepath(fname,root_dir='.',subdir=location)
  readf,1,err
  close,1


  openw,10,filepath('eigenvalues.dat',root_dir='.',subdir=location)
  openw,20,filepath('eigenvectors.dat',root_dir='.',subdir=location)
  openw,30,filepath('nonadia.dat',root_dir='.',subdir=location)
  openw,40,filepath('error.dat',root_dir='.',subdir=location)

  cnt = 0 
  for ngrid=0, nmodes-1 do begin
     
     growth = eigenvalues(0,ngrid)

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
  
     divv = dcomplex(nonadia(0,nbeg:nend), nonadia(1,nbeg:nend))
     dC   = dcomplex(nonadia(2,nbeg:nend), nonadia(3,nbeg:nend))

     integrand1 =imaginary(-conj(divv)*dC*exp(lnrho))
     integrand2 =imaginary(-(1d0-eps)*exp(lnrho)*conj(divv)*vrad*smallhg*smallq)

     denom     =int_tabulated(zaxis, exp(lnrho)*vmeri2)

     result = int_tabulated(zaxis, integrand1 + integrand2)

     result/=-2d0*eigenvalues(1,ngrid)*denom

     error = abs((result - growth)/growth)


     err1 = err(0,nbeg:nend)    ;mass error
     err2 = err(1,nbeg:nend)    ;vx error
     err3 = err(2,nbeg:nend)    ;vy error
     err4 = err(3,nbeg:nend)    ;vz error
     err5 = err(4,nbeg:nend)    ;energy error 
     

     if error lt tol then begin
        cnt += 1
        printf,10, eigenvalues(0,ngrid), eigenvalues(1,ngrid), format='(2(e22.15,x))'
        print, eigenvalues(0,ngrid), eigenvalues(1,ngrid), format='(2(e22.15,x))'
        for j=0, nz-1 do begin
           printf,20, real_part(bigW(j)), imaginary(bigW(j)), real_part(dfrac(j)), imaginary(dfrac(j)), $
                  real_part(vx(j)), imaginary(vx(j)), real_part(vy(j)), imaginary(vy(j)), real_part(vz(j)), $
                  imaginary(vz(j)), format='(10(e22.15,x))' 
           printf, 30, real_part(divv(j)), imaginary(divv(j)), real_part(dC(j)), imaginary(dC(j)), format='(4(e22.15,x))'
           printf, 40, abs(err1(j)),  abs(err2(j)), abs(err3(j)),  abs(err4(j)), abs(err5(j)), format='(5(e22.15,x))'

        endfor
     end

  endfor
  close,10
  close,20
  close,30
  close,40
  print, 'nmodes, out=', cnt 
end

  ;; err = dblarr(5, nlines)
  ;; fname = 'error.dat'
  ;; openr,1,filepath(fname,root_dir='.',subdir=location)
  ;; readf,1,err
  ;; close,1

  ;; err1 = err(0,nbeg:nend) ;mass error
  ;; err2 = err(1,nbeg:nend) ;vx error
  ;; err3 = err(2,nbeg:nend) ;vy error
  ;; err4 = err(3,nbeg:nend) ;vz error
  ;; err5 = err(4,nbeg:nend) ;energy error 

  ;; ytitle = 'normalized error'
  ;; fname = 'eigenvec_err.ps'
  ;; set_plot, 'ps'
  ;; device, filename=filepath(fname,root_dir='.',subdir=location) $
  ;;         ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  ;; plot, zaxis(1:nz-2),  err1(1:nz-2),xmargin=[8.5,1.5],ymargin=[3.5,0.5], ystyle=1  $
  ;;       ,charsize=2, thick=4, xrange=xrange $
  ;;       , xtitle=xtitle, ytitle=ytitle, yrange=[0.,max(err(*,nbeg+1:nend-1))] $
  ;;       ,ytickinterval=ytickinterval  $
  ;;       ,xtickinterval=xtickinterval, xstyle=1
  ;; oplot, zaxis(1:nz-2), err2(1:nz-2), thick=4, linestyle=1
  ;; oplot, zaxis(1:nz-2), err3(1:nz-2), thick=4, linestyle=2
  ;; oplot, zaxis(1:nz-2), err4(1:nz-2), thick=4, linestyle=3
  ;; oplot, zaxis(1:nz-2), err5(1:nz-2), thick=4, linestyle=4
  ;; device,/close
