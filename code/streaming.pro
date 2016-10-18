pro streaming, eta=eta, kx=kx, kz=kz, dgratio=dgratio, tstop=tstop

  ii    = dcomplex(0d0, 1d0)
  dfrac = dgratio/(1d0 + dgratio)
  ksq   = kx^2 + kz^2

  c0 = 2d0*kx*tstop*(1d0 - dfrac)*(1d0 - 2d0*dfrac)*kz^2/ksq
  
  c1 = kz^2/ksq

  c2 = -ii*tstop*( dfrac - 2d0*ii*kx*(1d0-dfrac)^2 )

  c3 = eta^2*(1d0 - 2d0*ii*kx*(1d0 - dfrac))/(ksq*(1d0-dfrac)) + 1d0
  c3 =-c3
  
  c4 = ii*dfrac*tstop

  c5 = eta^2/(1d0-dfrac)/ksq 

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

  work1 = sigma^2*(1d0 - sigma^2 - 2d0*ii*kx*(1d0-dfrac))
  work1/= kz^2 - sigma^2*ksq
  work1*= eta^2

  work2  = ii*(sigma + 2d0*kx*tstop*(1d0-dfrac)*(1d0-2d0*dfrac))
  work2 /= ii*sigma/(1d0-dfrac) - tstop*dfrac*ksq/eta^2

print, imaginary(work1), imaginary(work2)

stop


  work = imaginary(work)/freq(grid)

  print, 'freq, growth, work', freq(grid), growth(grid), work

   stheory = 4d0*dfrac*tstop^3*((1d0-2d0*dfrac)*(1d0-dfrac)*kx)^2*ksq/kz^2

;; ;  stheory = 4d0*dfrac*(1d0-dfrac)^2*(2d0*dfrac-1d0)^2*kx^4/kz^2
  
   freqtheory = 2d0*(1d0-dfrac)*(2d0*dfrac-1d0)*tstop*kx
  
   print, 's,om theory=', stheory, freqtheory



end
