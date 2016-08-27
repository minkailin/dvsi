pro streaming, eta=eta, kx=kx, kz=kz, dgratio=dgratio, tstop=tstop

  ii    = dcomplex(0d0, 1d0)
  dfrac = dgratio/(1d0 + dgratio)
  ksq   = kx^2 + kz^2

  c0 = 2d0*kx*tstop*(1d0 - dfrac)*(1d0 - 2d0*dfrac)*kz^2/ksq
  
  c1 = kz^2/ksq

  c2 = -ii*tstop*( dfrac - 2d0*ii*kx*(1d0-dfrac)^2 )

;  c2 = -ii*tstop*dfrac

  c3 = eta^2*(1d0 - 2d0*ii*kx*(1d0 - dfrac))/(ksq*(1d0-dfrac)) + 1d0
  c3 =-c3
  
  c4 = ii*dfrac*tstop

  c5 = eta^2/(1d0-dfrac)/ksq 

  coeffs = [c0, c1, c2, c3, c4, c5]

  roots = fz_roots(coeffs, /double, eps=1d-15)
  
  freq   = real_part(roots)
  growth = imaginary(roots)

  filter = where(abs(freq) gt 1d0)
  growth(filter) = 0d0 

  temp = max(growth, grid)

  sigma = dcomplex(freq(grid), growth(grid))

;  work = sigma^2*(1d0 - sigma^2 - 2d0*ii*kx*(1d0-dfrac))
;  work/= kz^2 - sigma^2*ksq
;  work*= eta^2

  work  = ii*(sigma + 2d0*kx*tstop*(1d0-dfrac)*(1d0-2d0*dfrac))
  work /= ii*sigma/(1d0-dfrac) - tstop*dfrac*ksq/eta^2


  work = imaginary(work)/freq(grid)

  print, 'freq, growth, work', freq(grid), growth(grid), work
;stop
end
