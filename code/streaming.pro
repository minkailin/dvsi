pro streaming, eta=eta, kx=kx, kz=kz, dgratio=dgratio, tstop=tstop

  ii    = dcomplex(0d0, 1d0)
  dfrac = dgratio/(1d0 + dgratio)
  ksq   = kx^2d0 + kz^2d0
  smallh = 0.05
  kappa2 = 1d0 - 6d0*(1d0-dfrac)*eta*smallh

  bta = 0d0 ;2d0*ii*(1d0-dfrac)*kx*eta^2d0/(dfrac*ksq) 


  kzsq  = kz^2d0/ksq 

  c0 = 2d0*kx*tstop*(1d0 - dfrac)*(1d0 - 2d0*dfrac)*kzsq*kappa2 
  
  c1 = kzsq*kappa2

  c2 = -ii*tstop*( dfrac*kappa2*(1d0-bta)  - 2d0*ii*kx*(1d0-dfrac)*(1d0 - dfrac*(1d0+bta)) )

  c3 = eta^2d0*(kappa2 - 2d0*ii*kx*(1d0 - dfrac))/(ksq*(1d0-dfrac)) + 1d0;*( 1d0 + 2d0*ii*eta^2d0/kx )
  c3 =-c3
  
  c4 = ii*dfrac*tstop*(1d0-bta)

  c5 = eta^2d0/(1d0-dfrac)/ksq 

  coeffs = [c0, c1, c2, c3, c4, c5]
 
; coeffs=[c2,c3,c4]


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

  lhs1 = eta^2d0*sigma^2d0*(kappa2 - sigma^2d0 - 2d0*ii*kx*(1d0-dfrac))
  rhs1 = kz^2d0*kappa2 - sigma^2d0*ksq
  
  Q1  = lhs1/rhs1 ;this is Q if W=1

  rhs2  = ii*(sigma + 2d0*kx*tstop*(1d0-dfrac)*(1d0-2d0*dfrac))
  lhs2  = ii*sigma/(1d0-dfrac) - tstop*dfrac*ksq/eta^2d0

  Q2 = rhs2/lhs2  ; this is Q if W = 1 

  vx = -ii*sigma*2d0*(1d0-dfrac)*eta + kx*sigma*Q1/eta 
  vx/= sigma^2d0 - kappa2 

  vz = Q1*kz/eta/sigma 

  ke = abs(vx)^2d0 + abs(vz)^2d0

  xi_x = ii*vx/sigma 
  lagrangian_Q = Q1 + xi_x*(-2d0*(1d0-dfrac)*eta)
  lagrangian_Q2 = Q2 + xi_x*(-2d0*(1d0-dfrac)*eta)
  
  work1 = imaginary(lagrangian_Q)*freq(grid)*1d4 
  work2 = imaginary(lagrangian_Q2)*freq(grid)*1d4 

  print, 'freq, growth, work', freq(grid), growth(grid), work1, work2

  sign = signum(freq(grid))
  
  re_Q = real_part(lagrangian_Q)
  im_Q = imaginary(lagrangian_Q)

  arg = atan2(im_Q, re_Q)  
  
  print, 'phase lag is', sign*arg*360.0/(2d0*!dpi)

  s_check = abs(sigma)^2d0*imaginary(lagrangian_Q) 
  s_check/= 2d0*freq(grid)*ke

  sig_check = kx*vx + kz*vz
  sig_check/=eta 

  print, 's_check', s_check

end
