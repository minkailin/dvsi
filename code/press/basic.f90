real*8 function lnrhog(zhat)
  use global
  implicit none 
  real*8, intent(in) :: zhat 

  lnrhog = -0.5d0*zhat**2d0 - dgratio*delta**2d0*(1d0 - exp(-zhat**2d0/2d0/delta**2d0) ) 
end function lnrhog

real*8 function dlnrhog_dr(zhat)
  use global
  implicit none 
  real*8, intent(in) :: zhat 
 
  dlnrhog_dr = smallh_g*rhog0_power + dlnHg_dlnr*smallh_g*( zhat**2d0 + 2d0*dgratio*delta**2d0*( 1d0 - exp(-zhat**2d0/2d0/delta**2d0) ) ) &
              + smalld*dgratio*smallh_g*delta**2d0*(1d0 - exp(-zhat**2d0/2d0/delta**2d0) )
end function dlnrhog_dr

real*8 function d2lnrhog_dr2(zhat)
  !this is (1/r)*(d/dr)( r*dlnrhog_dr )
  use global
  implicit none 
  real*8, intent(in) :: zhat 
  
  d2lnrhog_dr2 = -2d0*dlnHg_dlnr**2d0*smallh_g**2d0*( zhat**2d0 + 2d0*dgratio*delta**2d0*( 1d0 - exp(-zhat**2d0/2d0/delta**2d0) ) ) &
                 -4d0*smalld*dgratio*smallh_g**2d0*delta**2d0*dlnHg_dlnr*( 1d0 - exp(-zhat**2d0/2d0/delta**2d0) ) &
                 -(smalld*smallh_g*delta)**2d0*dgratio*( 1d0 - exp(-zhat**2d0/2d0/delta**2d0) )
end function d2lnrhog_dr2

real*8 function dlnrhog_dz(zhat)
  use global
  implicit none 
  real*8, intent(in) :: zhat 
  
  dlnrhog_dz = -zhat*(1d0 + dgratio*exp( - zhat**2d0/2d0/delta**2d0 ) )
end function dlnrhog_dz

real*8 function d2lnrhog_dz2(zhat)
  use global
  implicit none 
  real*8, intent(in) :: zhat 
  
  d2lnrhog_dz2 = -1d0 - dgratio*exp(-zhat**2d0/2d0/delta**2d0)*(1d0 - zhat**2d0/delta**2d0)
end function d2lnrhog_dz2

real*8 function del2_rhog(zhat)
  use global
  implicit none 
  real*8, intent(in) :: zhat 
  real*8, external  :: d2lnrhog_dr2, d2lnrhog_dz2

  del2_rhog = d2lnrhog_dr2(zhat) + d2lnrhog_dz2(zhat)
end function del2_rhog

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

real*8 function eps_tilde(zhat)
  !dust to gas ratio as function of z
  use global
  implicit none 
  real*8, intent(in) :: zhat 

  eps_tilde = dgratio*exp( -zhat**2d0/2d0/delta**2d0 )
end function eps_tilde

real*8 function deps_tilde_dr(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat 

  deps_tilde_dr = -smalld*dgratio*smallh_g*exp(-zhat**2d0/2d0/delta**2d0)
end function deps_tilde_dr

real*8 function d2eps_tilde_dr2(zhat)
  !this is (1/r)*(d/dr)*(r*deps_tilde_dr)
  use global
   implicit none
  real*8, intent(in) :: zhat 
  
  d2eps_tilde_dr2 = (smalld*smallh_g)**2d0*dgratio*exp(-zhat**2d0/2d0/delta**2d0) 
end function d2eps_tilde_dr2

real*8 function deps_tilde_dz(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: eps_tilde
  
  deps_tilde_dz = -zhat*eps_tilde(zhat)/delta**2d0
end function deps_tilde_dz

real*8 function d2eps_tilde_dz2(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: eps_tilde
  
  d2eps_tilde_dz2 = -(eps_tilde(zhat)/delta**2d0)*(1d0 - zhat**2d0/delta**2d0)
end function d2eps_tilde_dz2

real*8 function del2_eps_tilde(zhat)
  use global
   implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: d2eps_tilde_dr2, d2eps_tilde_dz2

  del2_eps_tilde = d2eps_tilde_dr2(zhat) + d2eps_tilde_dz2(zhat)
end function del2_eps_tilde

real*8 function del_eps_tilde_sq(zhat)
  use global
   implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: deps_tilde_dr, deps_tilde_dz

  del_eps_tilde_sq = deps_tilde_dr(zhat)**2d0 + deps_tilde_dz(zhat)**2d0
end function del_eps_tilde_sq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real*8 function lnrho(zhat)
  use global
   implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: lnrhog, eps_tilde

  lnrho = lnrhog(zhat) + log(1d0 + eps_tilde(zhat))
end function lnrho

real*8 function dlnrho_dr(zhat)
  use global
   implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: dlnrhog_dr, eps_tilde, deps_tilde_dr
  
  dlnrho_dr = dlnrhog_dr(zhat) + deps_tilde_dr(zhat)/(1d0 + eps_tilde(zhat))
end function dlnrho_dr

real*8 function dlnrho_dz(zhat)
  use global
   implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: dlnrhog_dz, eps_tilde, deps_tilde_dz
  
  dlnrho_dz =dlnrhog_dz(zhat) + deps_tilde_dz(zhat)/(1d0 + eps_tilde(zhat))
end function dlnrho_dz


real*8 function d2lnrho_dz2(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: eps_tilde, deps_tilde_dz, d2eps_tilde_dz2, d2lnrhog_dz2
  
  d2lnrho_dz2 = d2lnrhog_dz2(zhat) + d2eps_tilde_dz2(zhat)/(1d0+eps_tilde(zhat)) -deps_tilde_dz(zhat)**2d0/(1d0 + eps_tilde(zhat))**2d0
end function d2lnrho_dz2

real*8 function del2_lnrho(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: del2_rhog, del2_eps_tilde, del_eps_tilde_sq, eps_tilde
  
  del2_lnrho = del2_rhog(zhat) + del2_eps_tilde(zhat)/(1d0 + eps_tilde(zhat)) - del_eps_tilde_sq(zhat)/(1d0 + eps_tilde(zhat))**2d0
end function del2_lnrho

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real*8 function eps(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: eps_tilde

  eps = eps_tilde(zhat)/(1d0 + eps_tilde(zhat))
end function eps

real*8 function deps_dr(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: eps_tilde, deps_tilde_dr

  deps_dr = deps_tilde_dr(zhat)/(1d0 + eps_tilde(zhat))**2d0
end function deps_dr

real*8 function deps_dz(zhat)
  use global
   implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: eps_tilde, deps_tilde_dz

  deps_dz = deps_tilde_dz(zhat)/(1d0 + eps_tilde(zhat))**2d0
end function deps_dz

real*8 function del2_eps(zhat)
  use global
   implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: eps_tilde, del_eps_tilde_sq, del2_eps_tilde 

  del2_eps = del2_eps_tilde(zhat)/(1d0 + eps_tilde(zhat))**2d0 - 2d0*del_eps_tilde_sq(zhat)/(1d0 + eps_tilde(zhat))**3d0
end function del2_eps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real*8 function del_lnrho_dot_del_eps(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: dlnrho_dr, deps_dr, dlnrho_dz, deps_dz
  
  del_lnrho_dot_del_eps = dlnrho_dr(zhat)*deps_dr(zhat) +  dlnrho_dz(zhat)*deps_dz(zhat)
end function del_lnrho_dot_del_eps

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real*8 function Fr(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: eps, deps_dr, dlnrho_dr

  Fr = -(1d0-eps(zhat))*smallh_g*smallq + deps_dr(zhat) - (1d0-eps(zhat))*dlnrho_dr(zhat)
end function Fr

real*8 function Fz(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: eps, deps_dz, dlnrho_dz
  
  Fz = deps_dz(zhat) - (1d0-eps(zhat))*dlnrho_dz(zhat)
end function Fz

real*8 function div_F(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: eps, deps_dr, del2_eps, dlnrho_dr, del_lnrho_dot_del_eps, del2_lnrho 
  
  div_F = 2d0*deps_dr(zhat)*smallh_g*smallq - (1d0 - eps(zhat))*(smallh_g*smallq)**2d0 + del2_eps(zhat)
  div_F = div_F - (1d0 - eps(zhat))*smallh_g*smallq*dlnrho_dr(zhat) + del_lnrho_dot_del_eps(zhat)
  div_F = div_F - (1d0 - eps(zhat))*del2_lnrho(zhat)
end function div_F

real*8 function F_dot_deleps(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: deps_dz, deps_dr, Fr, Fz
  
  F_dot_deleps = Fr(zhat)*deps_dr(zhat) + Fz(zhat)*deps_dz(zhat)
end function F_dot_deleps


real*8 function F_dot_dellncs2(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: Fr
  
  F_dot_dellncs2 = Fr(zhat)*smallh_g*smallq 
end function F_dot_dellncs2

real*8 function gr(zhat)
  !grad P/rho in r 
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: eps, dlnrhog_dr, Fr

  gr = (1d0 - eps(zhat))*smallh_g*smallq + (1d0 - eps(zhat))*dlnrhog_dr(zhat)
end function Gr

real*8 function gz(zhat)
  !grad P/rho in z 
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: eps, dlnrhog_dz, Fz
  
   gz = (1d0 - eps(zhat))*dlnrhog_dz(zhat)
end function Gz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real*8 function rotation(zhat)
  use global
   implicit none
  real*8, intent(in) :: zhat
  real*8, external  :: dlnrhog_dr, eps_tilde
  
  rotation = 1d0 - (3d0/2d0)*smallh_g**2d0*zhat**2d0 + smallh_g**2d0*(smallq +  dlnrhog_dr(zhat)/smallh_g)/(1d0 + eps_tilde(zhat))
  rotation = sqrt(rotation)
end function rotation

real*8 function kappa_sq(zhat)
  use global
   implicit none
  real*8, intent(in) :: zhat
  real*8 :: om2, term1, term2, term3 
  real*8, external  :: rotation, dlnrhog_dr, d2lnrhog_dr2, eps_tilde, deps_tilde_dr
  
  om2 = rotation(zhat)**2d0

  term1 = -2d0 + 6d0*smallh_g**2d0*zhat**2d0
  term2 = d2lnrhog_dr2(zhat)/(1d0 + eps_tilde(zhat))
  term3 = smallh_g**2d0*smallq/(1d0 + eps_tilde(zhat))**2d0 - smallh_g*deps_tilde_dr(zhat)/(1d0 + eps_tilde(zhat))**2d0 - smallh_g**2d0/(1d0+eps_tilde(zhat))  
  term3 = (smallq + dlnrhog_dr(zhat)/smallh_g)*term3 

  kappa_sq = 3d0*om2 + term1 + term2 + term3 
end function kappa_sq

real*8 function vertical_shear(zhat)
  use global
   implicit none
  real*8, intent(in) :: zhat
  real*8, external  :: dlnrho_dr, deps_dz, dlnrho_dz, deps_dr, eps 

  vertical_shear = -dlnrho_dr(zhat)*deps_dz(zhat) - dlnrho_dz(zhat)*( smallh_g*smallq*(1d0-eps(zhat)) - deps_dr(zhat) )

end function vertical_shear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real*8 function stopping_time(zhat)
  use global
   implicit none
  real*8, intent(in) :: zhat
  real*8, external  :: lnrho 

  stopping_time = tstop*exp(lnrho(0d0))/exp(lnrho(zhat))
end function stopping_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chebyshev_poly(l, zbar, T_l, dT_l, d2T_l)
  implicit none
  real*8, intent(in) :: l, zbar
  real*8, intent(out):: T_l, dT_l, d2T_l
  real*8 :: t, lt, lsq
  lsq = l*l
  t = acos(zbar)
  lt = l*t
  T_l = cos(lt)
  if(abs(zbar).lt.1d0) then
     dT_l = l*sin(lt)/sin(t)
     d2T_l= -lsq*cos(lt)/sin(t)**2 + l*cos(t)*sin(lt)/sin(t)**3
  else
     dT_l =lsq
     d2t_l =lsq*(lsq-1d0)/3d0
     if(zbar.eq.-1d0) then
        dT_l = (-1d0)**(l+1d0)*dT_l
        d2T_l= (-1d0)**(l+2d0)*d2T_l
     endif
  endif
  return
end subroutine chebyshev_poly
