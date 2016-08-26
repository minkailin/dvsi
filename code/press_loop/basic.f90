real*8 function lnrhog(zhat)
  use global
  implicit none 
  real*8, intent(in) :: zhat 

if(vstruct.eq.'mlin') lnrhog = -0.5d0*zhat**2d0 - dgratio*delta**2d0*(1d0 - exp(-zhat**2d0/2d0/delta**2d0) ) 
if(vstruct.eq.'tl02') lnrhog  = -0.5d0*zhat**2d0
end function lnrhog

real*8 function dlnrhog_dr(zhat)
  use global
  implicit none 
  real*8, intent(in) :: zhat 
 
if(vstruct.eq.'mlin') then 
  dlnrhog_dr = smallh_g*rhog0_power + dlnHg_dlnr*smallh_g*( zhat**2d0 + 2d0*dgratio*delta**2d0*( 1d0 - exp(-zhat**2d0/2d0/delta**2d0) ) ) &
              + smalld*dgratio*smallh_g*delta**2d0*(1d0 - exp(-zhat**2d0/2d0/delta**2d0) )
endif

if(vstruct.eq.'tl02')  dlnrhog_dr = smallh_g*rhog0_power + dlnHg_dlnr*smallh_g*zhat**2d0
end function dlnrhog_dr

real*8 function d2lnrhog_dr2(zhat)
  !this is (1/r)*(d/dr)( r*dlnrhog_dr )
  use global
  implicit none 
  real*8, intent(in) :: zhat 
  
if(vstruct.eq.'mlin') then
  d2lnrhog_dr2 = -2d0*dlnHg_dlnr**2d0*smallh_g**2d0*( zhat**2d0 + 2d0*dgratio*delta**2d0*( 1d0 - exp(-zhat**2d0/2d0/delta**2d0) ) ) &
                 -4d0*smalld*dgratio*smallh_g**2d0*delta**2d0*dlnHg_dlnr*( 1d0 - exp(-zhat**2d0/2d0/delta**2d0) ) &
                 -(smalld*smallh_g*delta)**2d0*dgratio*( 1d0 - exp(-zhat**2d0/2d0/delta**2d0) )
endif  

if(vstruct.eq.'tl02') d2lnrhog_dr2 = -2d0*dlnHg_dlnr**2d0*smallh_g**2d0*zhat**2d0
end function d2lnrhog_dr2

real*8 function dlnrhog_dz(zhat)
  use global
  implicit none 
  real*8, intent(in) :: zhat 
  
if(vstruct.eq.'mlin')  dlnrhog_dz = -zhat*(1d0 + dgratio*exp( - zhat**2d0/2d0/delta**2d0 ) )
if(vstruct.eq.'tl02')  dlnrhog_dz = -zhat
end function dlnrhog_dz

real*8 function d2lnrhog_dz2(zhat)
  use global
  implicit none 
  real*8, intent(in) :: zhat 
  
if(vstruct.eq.'mlin')  d2lnrhog_dz2 = -1d0 - dgratio*exp(-zhat**2d0/2d0/delta**2d0)*(1d0 - zhat**2d0/delta**2d0)
if(vstruct.eq.'tl02')  d2lnrhog_dz2 = -1d0
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

if(vstruct.eq.'mlin')  eps_tilde = dgratio*exp( -zhat**2d0/2d0/delta**2d0 )
if(vstruct.eq.'tl02')  eps_tilde = dgratio*exp( - beta**2d0*(exp(zhat**2d0/2d0) - 1d0) )
end function eps_tilde

real*8 function deps_tilde_dr(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external :: eps_tilde

if(vstruct.eq.'mlin')  deps_tilde_dr =(-smalld*dgratio + varHd*(zhat**2d0/delta**2d0)*dlnHg_dlnr)*smallh_g*exp(-zhat**2d0/2d0/delta**2d0)
if(vstruct.eq.'tl02')  deps_tilde_dr = eps_tilde(zhat)*( -smalld*smallh_g + beta**2d0*zhat**2d0*smallh_g*dlnHg_dlnr*exp(zhat**2d0/2d0) )
end function deps_tilde_dr

real*8 function d2eps_tilde_dr2(zhat)
  !this is (1/r)*(d/dr)*(r*deps_tilde_dr)
  use global
   implicit none
  real*8, intent(in) :: zhat 
  real*8, external :: eps_tilde, deps_tilde_dr
  real*8 :: temp, delta2 
  
if(vstruct.eq.'mlin')  then 
!   d2eps_tilde_dr2 = (smalld*smallh_g)**2d0*dgratio*exp(-zhat**2d0/2d0/delta**2d0) 
   d2eps_tilde_dr2  = smalld**2d0*dgratio - varHd*(zhat**2d0/delta**2d0)*dlnHg_dlnr**2d0
   d2eps_tilde_dr2  = d2eps_tilde_dr2 + varHd*( -smalld*dgratio + (zhat**2d0/delta**2d0)*dlnHg_dlnr )*zhat**2d0*dlnHg_dlnr/delta**2d0 
   d2eps_tilde_dr2  = d2eps_tilde_dr2*smallh_g**2d0*exp(-zhat**2d0/2d0/delta**2d0) 
endif

if(vstruct.eq.'tl02') then 
  temp = beta**2d0*zhat**2d0*smallh_g**2d0*(-2d0 - zhat**2d0)*dlnHg_dlnr**2d0*exp(zhat**2d0/2d0)
  d2eps_tilde_dr2 = eps_tilde(zhat)*temp + deps_tilde_dr(zhat)**2d0/eps_tilde(zhat)
endif 
end function d2eps_tilde_dr2

real*8 function deps_tilde_dz(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: eps_tilde
  
if(vstruct.eq.'mlin')  deps_tilde_dz = -zhat*eps_tilde(zhat)/delta**2d0

if(vstruct.eq.'tl02') then
  deps_tilde_dz  = -beta**2d0*exp(zhat**2d0/2d0)*zhat 
  deps_tilde_dz  = deps_tilde_dz*eps_tilde(zhat) 
endif

end function deps_tilde_dz

real*8 function d2eps_tilde_dz2(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat 
  real*8, external  :: eps_tilde, deps_tilde_dz
  real*8 :: temp 

if(vstruct.eq.'mlin')  d2eps_tilde_dz2 = -(eps_tilde(zhat)/delta**2d0)*(1d0 - zhat**2d0/delta**2d0)

if(vstruct.eq.'tl02') then
  temp = -beta**2d0*exp(zhat**2d0/2d0)*(1d0 + zhat**2d0)
  d2eps_tilde_dz2 = deps_tilde_dz(zhat)**2d0/eps_tilde(zhat) + eps_tilde(zhat)*temp 
endif 
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cubint ( ntab, xtab, ftab, ia, ib, result, error )

!*****************************************************************************80
!
!! CUBINT approximates an integral using cubic interpolation of data.
!
!  Discussion:
!
!    The integral to be approximated is
! 
!      Integral ( XTAB(IB) <= X <= XTAB(IA) ) F(X) DX
!
!    The routine estimates the error in integration.
!
!  Modified:
!
!    10 February 2006
!
!  Reference:
!
!    Philip Davis, Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Second Edition,
!    Dover, 2007,
!    ISBN: 0486453391,
!    LC: QA299.3.D28.
!
!    Philip Gill, GF Miller,
!    An algorithm for the integration of unequally spaced data,
!    The Computer Journal, 
!    Number 15, Number 1, 1972, pages 80-83.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NTAB, the number of tabulated points.
!    NTAB must be at least 4.
!
!    Input, real ( kind = 8 ) XTAB(NTAB), contains the points at which the
!    function was tabulated.  XTAB should contain distinct
!    values, given in ascending order.
!
!    Input, real ( kind = 8 ) FTAB(NTAB), contains the tabulated function
!    values, FTAB(I) = F(XTAB(I)).
!
!    Input, integer ( kind = 4 ) IA, the entry of XTAB at which integration
!    is to begin.  IA must be no less than 1 and no greater
!    than NTAB.
!
!    Input, integer ( kind = 4 ) IB, the entry of XTAB at which integration
!    is to end.  IB must be no less than 1 and no greater than
!    NTAB.
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the
!    integral from XTAB(IA) to XTAB(IB) of the function.
!
!    Output, real ( kind = 8 ) ERROR, an estimate of the error in
!    integration.
!
  implicit none

  integer ( kind = 4 ) ntab

  real ( kind = 8 ) c
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) d3
  real ( kind = 8 ) error
  real ( kind = 8 ) ftab(ntab)
  real ( kind = 8 ) h1
  real ( kind = 8 ) h2
  real ( kind = 8 ) h3
  real ( kind = 8 ) h4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) it
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r3
  real ( kind = 8 ) r4
  real ( kind = 8 ) result
  real ( kind = 8 ) s
  real ( kind = 8 ) term
  real ( kind = 8 ) xtab(ntab)

  result = 0.0D+00
  error = 0.0D+00
 
  if ( ia == ib ) then
    return
  end if
 
  if ( ntab < 4 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i8)' ) '  NTAB must be at least 4, but input NTAB = ', ntab
    stop 1
  end if
 
  if ( ia < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i8)' ) '  IA must be at least 1, but input IA = ', ia
    stop 1
  end if
 
  if ( ntab < ia ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i8)' ) '  IA must be <= NTAB, but input IA = ', ia
    stop 1
  end if
 
  if ( ib < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i8)' ) '  IB must be at least 1, but input IB = ', ib
    stop 1
  end if
 
  if ( ntab < ib ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBINT - Fatal error!'
    write ( *, '(a,i8)' ) '  IB must be <= NTAB, but input IB = ', ib
    stop 1
  end if
!
!  Temporarily switch IA and IB, and store minus sign in IND
!  so that, while integration is carried out from low X's
!  to high ones, the sense of the integral is preserved.
!
  if ( ib < ia ) then
    ind = -1
    it = ib
    ib = ia
    ia = it
  else
    ind = 1
  end if
 
  s = 0.0D+00
  c = 0.0D+00
  r4 = 0.0D+00
  j = ntab-2
  if ( ia < ntab-1 .or. ntab == 4 ) then
    j = max ( 3, ia )
  end if

  k = 4
  if ( 2 < ib .or. ntab == 4 ) then
    k = min ( ntab, ib + 2 ) - 1
  end if
 
  do i = j, k
 
    if ( i <= j ) then
 
      h2 = xtab(j-1) - xtab(j-2)
      d3 = ( ftab(j-1) - ftab(j-2) ) / h2
      h3 = xtab(j) - xtab(j-1)
      d1 = ( ftab(j) - ftab(j-1) ) / h3
      h1 = h2 + h3
      d2 = ( d1 - d3 ) / h1
      h4 = xtab(j+1) - xtab(j)
      r1 = ( ftab(j+1) - ftab(j) ) / h4
      r2 = ( r1 - d1 ) / ( h4 + h3 )
      h1 = h1 + h4
      r3 = (r2-d2) / h1
 
      if ( ia <= 1 ) then
        result = h2 * ( ftab(1) + h2 * ( 0.5D+00 * d3 - h2 &
          * ( d2 / 6.0D+00 -(h2+h3+h3)*r3/12.0D+00)))
        s = -h2**3 * (h2*(3.0D+00*h2+5.0D+00*h4)+10.0D+00*h3*h1) / 60.0D+00
      end if
 
    else
 
      h4 = xtab(i+1) - xtab(i)
      r1 = ( ftab(i+1) - ftab(i) ) / h4
      r4 = h4 + h3
      r2 = ( r1 - d1 ) / r4
      r4 = r4 + h2
      r3 = ( r2 - d2 ) / r4
      r4 = ( r3 - d3 ) / ( r4 + h1 )
 
    end if
 
    if ( ia < i .and. i <= ib ) then
 
      term = h3 * ( ( ftab(i) + ftab(i-1) ) * 0.5D+00 &
        -h3 * h3 * ( d2 + r2 + ( h2 - h4 ) * r3 ) / 12.0D+00 )
      result = result + term
      c = h3**3 * ( 2.0D+00 * h3 * h3 &
        + 5.0D+00 * ( h3 * ( h4 + h2 ) + 2.0D+00 * h2 * h4 ) ) / 120.0D+00
      error = error + (c+s)*r4
 
      if ( i /= j ) then
        s = c
      else
        s = s + c + c
      end if
 
    else
 
      error = error + r4 * s
 
    end if
 
    if ( k <= i ) then
 
      if ( ntab <= ib ) then
        term = h4 * ( ftab(ntab) - h4 * ( 0.5 * r1 &
          + h4 * ( r2 / 6.0D+00 + ( h3 + h3 + h4 ) * r3 / 12.0D+00 )))
        result = result + term
        error = error - h4**3 * r4 * &
          ( h4 * ( 3.0D+00 * h4 + 5.0D+00 * h2 ) &
          + 10.0D+00 * h3 * ( h2 + h3 + h4 ) ) / 60.0D+00
      end if
 
      if ( ntab-1 <= ib ) then
        error = error + s * r4
      end if

    else

      h1 = h2
      h2 = h3
      h3 = h4
      d1 = r1
      d2 = r2
      d3 = r3
    end if
 
  end do
!
!  Restore original values of IA and IB, reverse signs
!  of RESULT and ERROR, to account for integration
!  that proceeded from high X to low X.
!
  if ( ind /= 1 ) then
    it = ib
    ib = ia
    ia = it
    result = -result
    error = -error
  end if
 
  return
end
