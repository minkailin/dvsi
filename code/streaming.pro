pro streaming, eta=eta, kx=kx, kz=kz, dgratio=dgratio, tstop=tstop

  ii    = dcomplex(0d0, 1d0)
  dfrac = dgratio/(1d0 + dgratio)
  ksq   = kx^2d0 + kz^2d0

  kzsq  = kz^2d0/ksq 

  c0 = 2d0*kx*tstop*(1d0 - dfrac)*(1d0 - 2d0*dfrac)*kzsq
  
  c1 = kzsq

  c2 = -ii*tstop*( dfrac - 2d0*ii*kx*(1d0-dfrac)^2d0 )

  c3 = eta^2d0*(1d0 - 2d0*ii*kx*(1d0 - dfrac))/(ksq*(1d0-dfrac)) + 1d0
  c3 =-c3
  
  c4 = ii*dfrac*tstop

  c5 = eta^2d0/(1d0-dfrac)/ksq 

  coeffs = [c0, c1, c2, c3, c4, c5]
 
;  c0 = 2.0*kx*tstop*(1d0 - dfrac)*(1d0 - 2d0*dfrac)
;  c1 = 1d0
;  c2 = -ii*dfrac*tstop*(ksq/kz^2)
;  coeffs = [c0,c1,c2]

  roots = fz_roots(coeffs, /double, eps=1d-15)
  
  freq   = real_part(roots)
  growth = imaginary(roots)
  
 ; filter = where(abs(freq) gt 1d0)
 ; growth(filter) = 0d0 

  temp = max(growth, grid)

  sigma = dcomplex(freq(grid), growth(grid))
  
;;  stop

  lhs1 = eta^2d0*sigma^2d0*(1d0 - sigma^2d0 - 2d0*ii*kx*(1d0-dfrac))
  rhs1 = kz^2d0 - sigma^2d0*ksq

  work1  = lhs1/rhs1 

  rhs2  = ii*(sigma + 2d0*kx*tstop*(1d0-dfrac)*(1d0-2d0*dfrac))
  lhs2  = ii*sigma/(1d0-dfrac) - tstop*dfrac*ksq/eta^2d0

  work2 = rhs2/lhs2 

  work = work1 

  work = imaginary(work)/freq(grid)

  print, 'freq, growth, work', freq(grid), growth(grid), work

   stheory = 4d0*dfrac*tstop^3d0*((1d0-2d0*dfrac)*(1d0-dfrac)*kx)^2d0*ksq/kz^2d0

;; ;  stheory = 4d0*dfrac*(1d0-dfrac)^2*(2d0*dfrac-1d0)^2*kx^4/kz^2
  
   freqtheory = 2d0*(1d0-dfrac)*(2d0*dfrac-1d0)*tstop*kx
  
   print, 's,om theory=', stheory, freqtheory



end
