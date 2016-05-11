module global
  implicit none 
  integer, parameter :: nvar = 5 
  real*8, parameter :: pi = 2d0*acos(0d0)
  complex*16, parameter :: ii = (0d0, 1d0), zero = (0d0, 0d0)

  integer :: nz, bignz, Nzmid
  real*8 :: smallh_g, rhog0_power, smallq, dgratio, delta, smalld, tstop, Hd
  real*8 :: zmax 
  real*8 :: dlnHg_dlnr
  real*8 :: krad, bc_tol 

  real*8, allocatable :: zaxis(:), lnrho_arr(:), eps_arr(:), dlnrhodr_arr(:), depsdr_arr(:),dlnrhodz_arr(:), depsdz_arr(:) 
  real*8, allocatable :: del2lnrho_arr(:), del2eps_arr(:), dellnrho_dot_deleps_arr(:), tstop_arr(:), d2lnrhodz2_arr(:)
  real*8, allocatable :: Fr_arr(:), Fz_arr(:), F_dot_deleps_arr(:), divF_arr(:), F_dot_dellncs2_arr(:)
  real*8, allocatable :: omega_arr(:), kappa2_arr(:), vshear_arr(:) 

  complex*16, allocatable     :: T(:,:), Tp(:,:), Tpp(:,:)

  complex*16, allocatable :: L11(:,:), L12(:,:), L13(:,:), L14(:,:), L15(:,:)
  complex*16, allocatable :: L21(:,:), L22(:,:), L23(:,:), L24(:,:), L25(:,:)
  complex*16, allocatable :: L31(:,:), L32(:,:), L33(:,:), L34(:,:), L35(:,:)
  complex*16, allocatable :: L41(:,:), L42(:,:), L43(:,:), L44(:,:), L45(:,:)
  complex*16, allocatable :: L51(:,:), L52(:,:), L53(:,:), L54(:,:), L55(:,:)
  complex*16, allocatable :: R11(:,:), R12(:,:), R13(:,:), R14(:,:), R15(:,:) 
  complex*16, allocatable :: R21(:,:), R22(:,:), R23(:,:), R24(:,:), R25(:,:)
  complex*16, allocatable :: R31(:,:), R32(:,:), R33(:,:), R34(:,:), R35(:,:)
  complex*16, allocatable :: R41(:,:), R42(:,:), R43(:,:), R44(:,:), R45(:,:) 
  complex*16, allocatable :: R51(:,:), R52(:,:), R53(:,:), R54(:,:), R55(:,:) 
  complex*16, allocatable :: bigW(:), dfrac(:), vx(:), vy(:), vz(:)

  complex*16, allocatable :: bigmatrix_lhs(:,:)
  complex*16, allocatable :: bigmatrix_rhs(:,:)
  complex*16, allocatable :: eigenvec(:), eigenvals(:)

end module global

program dvsi
  use global
  implicit none 
  integer :: i, j, k, lmax 
  real*8 :: z, zbar, m, T_l, Tp_l, Tpp_l, dT_l, d2T_l
  real*8, external  :: lnrho, eps, dlnrho_dr, deps_dr, dlnrho_dz, deps_dz
  real*8, external  :: del2_lnrho, del2_eps, del_lnrho_dot_del_eps, d2lnrho_dz2 
  real*8, external  :: rotation, kappa_sq, vertical_shear
  real*8, external  :: Fr, Fz, F_dot_deleps, div_F, F_dot_dellncs2, gr, gz
  real*8, external  :: stopping_time 
  namelist /params/ smallh_g, rhog0_power, smallq, Hd, dgratio, smalld, tstop 
  namelist /grid/ zmax, nz 
  namelist /mode/ krad

  !read input parameters.
  open(7, file="input")
  read(7, nml=params)
  read(7, nml=grid)
  read(7, nml=mode)
  close(7)

  open(10,file='params.dat') 
  write(10,fmt='(4(e22.15,x))'), dble(nz), smallh_g, krad, Hd
  close(10)

  !check we have odd grid points so one grid point coincides with midplane 
  if(mod(nz,2).eq.0) then
     print*, 'Nz needs to be odd but Nz=', nz
     stop
  else
     Nzmid = (Nz+1)/2
  endif
  
  !calculate secondary parameters 
  bignz      = nvar*nz 
  dlnHg_dlnr = (3d0 + smallq)/2d0 
  delta      = sqrt( 1d0/( 1d0/Hd**2d0 - 1d0 ) )
  

  !allocate physical array
  
  allocate(zaxis(nz))
  allocate(lnrho_arr(nz))
  allocate(dlnrhodr_arr(nz))
  allocate(dlnrhodz_arr(nz))
  allocate(d2lnrhodz2_arr(nz))
  allocate(del2lnrho_arr(nz))
  allocate(eps_arr(nz))
  allocate(depsdr_arr(nz))
  allocate(depsdz_arr(nz))
  allocate(del2eps_arr(nz))
  allocate(dellnrho_dot_deleps_arr(nz))
  allocate(tstop_arr(nz))
  allocate(omega_arr(nz))
  allocate(kappa2_arr(nz))
  allocate(vshear_arr(nz))
  allocate(Fr_arr(nz))
  allocate(Fz_arr(nz))
  allocate(F_dot_deleps_arr(nz))
  allocate(divF_arr(nz))
  allocate(F_dot_dellncs2_arr(nz))

  !setup physical grid 
  
  lmax = nz-1
  zaxis(Nzmid) = 0d0
  do j=Nzmid+1, nz
     zaxis(j) = -zmax*cos(pi*(j-1d0)/lmax)
  enddo
  do j=1, Nzmid-1
     zaxis(j) = -zaxis(Nz - j + 1)
  enddo
  
  do k=1, nz 
     z = zaxis(k)
     
     lnrho_arr(k)               = lnrho(z)
     dlnrhodr_arr(k)            = dlnrho_dr(z)
     dlnrhodz_arr(k)            = dlnrho_dz(z)
     d2lnrhodz2_arr(k)          = d2lnrho_dz2(z)
     del2lnrho_arr(k)           = del2_lnrho(z)
     
     eps_arr(k)                 = eps(z)
     depsdr_arr(k)              = deps_dr(z)
     depsdz_arr(k)              = deps_dz(z)
     del2eps_arr(k)             = del2_eps(z)

     dellnrho_dot_deleps_arr(k) = del_lnrho_dot_del_eps(z)

     tstop_arr(k)  = stopping_time(z)

     omega_arr(k)  = rotation(z)  
     kappa2_arr(k) = kappa_sq(z)
     vshear_arr(k) = vertical_shear(z)

     Fr_arr(k)              = Fr(z)
     Fz_arr(k)              = Fz(z)
     F_dot_deleps_arr(k)    = F_dot_deleps(z)
     divF_arr(k)            = div_F(z)
     F_dot_dellncs2_arr(k)  = F_dot_dellncs2(z) 
  end do
  
  !output basic state 
  open(10,file='basic.dat')
  do k=1, nz
     write(10,fmt='(7(e22.15,x))'), zaxis(k), lnrho_arr(k), eps_arr(k), tstop_arr(k), omega_arr(k)**2d0, kappa2_arr(k), vshear_arr(k) 
  enddo
  close(10)
  
  !allocate matrices for eigenvalue problem 
  
  allocate(T(nz,nz))
  allocate(Tp(nz,nz))
  allocate(Tpp(nz,nz))
  allocate(bigmatrix_rhs(bignz, bignz))
  allocate(bigmatrix_lhs(bignz, bignz))

  allocate(L11(nz,nz)) ;  allocate(L12(nz,nz)) ;  allocate(L13(nz,nz)) ;  allocate(L14(nz,nz)) ;  allocate(L15(nz,nz)) 
  allocate(L21(nz,nz)) ;  allocate(L22(nz,nz)) ;  allocate(L23(nz,nz)) ;  allocate(L24(nz,nz)) ;  allocate(L25(nz,nz))
  allocate(L31(nz,nz)) ;  allocate(L32(nz,nz)) ;  allocate(L33(nz,nz)) ;  allocate(L34(nz,nz)) ;  allocate(L35(nz,nz))
  allocate(L41(nz,nz)) ;  allocate(L42(nz,nz)) ;  allocate(L43(nz,nz)) ;  allocate(L44(nz,nz)) ;  allocate(L45(nz,nz))
  allocate(L51(nz,nz)) ;  allocate(L52(nz,nz)) ;  allocate(L53(nz,nz)) ;  allocate(L54(nz,nz)) ;  allocate(L55(nz,nz))
  
  allocate(R11(nz,nz)) ;  allocate(R12(nz,nz)) ;  allocate(R13(nz,nz)) ;  allocate(R14(nz,nz)) ;  allocate(R15(nz,nz)) 
  allocate(R21(nz,nz)) ;  allocate(R22(nz,nz)) ;  allocate(R23(nz,nz)) ;  allocate(R24(nz,nz)) ;  allocate(R25(nz,nz))
  allocate(R31(nz,nz)) ;  allocate(R32(nz,nz)) ;  allocate(R33(nz,nz)) ;  allocate(R34(nz,nz)) ;  allocate(R35(nz,nz))
  allocate(R41(nz,nz)) ;  allocate(R42(nz,nz)) ;  allocate(R43(nz,nz)) ;  allocate(R44(nz,nz)) ;  allocate(R45(nz,nz))
  allocate(R51(nz,nz)) ;  allocate(R52(nz,nz)) ;  allocate(R53(nz,nz)) ;  allocate(R54(nz,nz)) ;  allocate(R55(nz,nz))

  allocate(bigW(nz))  
  allocate(dfrac(nz))
  allocate(vx(nz))
  allocate(vy(nz))
  allocate(vz(nz))

  !chebyshev polynomials 
  do j=1, nz !jth physical grid
     zbar = zaxis(j)/zmax
     do k=1, nz !kth basis
        m = dble(k-1d0)
        call chebyshev_poly(m,     zbar, T_l, dT_l, d2T_l)
        T(j,k)  =   T_l
        Tp(j,k) =  dT_l/zmax        !convert to deriv wrt physical grid
        Tpp(j,k)= d2T_l/zmax**2d0   !convert to deriv wrt physical grid
     enddo
  enddo
 
  call fill_bigmatrix(krad)
  call eigenvalue_problem

end program dvsi

subroutine fill_bigmatrix(knum)
  use global
  implicit none
  integer :: i, ibeg, iend 
  real*8, intent(in) :: knum 
  real*8 :: kx, z, gr_at_z, gz_at_z
  real*8, external :: gr, gz

  !convert k to code units (in Hgas, input is in H_dust) 

  kx = knum/Hd 

  !fill sub-matrices 
  do i=1, nz 

     z = zaxis(i)
     gr_at_z = gr(z)
     gz_at_z = gz(z)

     !continuity equation 
     L11(i,:) = (0d0, 0d0)
     L12(i,:) = (0d0, 0d0)
     L13(i,:) = -( ii*kx + dlnrhodr_arr(i))*T(i,:)
     L14(i,:) = (0d0, 0d0)
     L15(i,:) =( -Tp(i,:) - dlnrhodz_arr(i)*T(i,:) )
     
     R11(i,:) = T(i,:)
     R12(i,:) = zero
     R13(i,:) = zero 
     R14(i,:) = zero
     R15(i,:) = zero

     !dust continuity equation (equiv. energy equation) 
     L21(i,:) = - tstop_arr(i)*( &
          Fr_arr(i)*( &
          -depsdr_arr(i) - (1d0-eps_arr(i))*dlnrhodr_arr(i) + ii*kx*(1d0-eps_arr(i)) &
          )*T(i,:) &
          + Fz_arr(i)*(&
          -depsdz_arr(i)*T(i,:) + (1d0-eps_arr(i))*Tp(i,:) &
          ) &
          + depsdr_arr(i)*gr_at_z*T(i,:) &
          + depsdz_arr(i)*gz_at_z*T(i,:) &
          + divF_arr(i)*(1d0 - eps_arr(i))*T(i,:) &
          + eps_arr(i)*(&
          (ii*kx - dlnrhodr_arr(i))*gr_at_z*T(i,:) + Tp(i,:)*gz_at_z - T(i,:)*divF_arr(i) &
          ) &
          -0.5d0*Fr_arr(i)*smallh_g*smallq*(1d0-eps_arr(i))*T(i,:) &
          -0.5d0*eps_arr(i)*smallh_g*smallq*gr_at_z*T(i,:) &
          )
          
     L22(i,:) = - tstop_arr(i)*( &
          Fr_arr(i)*( &
          -ii*kx + smallh_g*smallq + dlnrhodr_arr(i) &
          )*T(i,:) &
          + Fz_arr(i)*(&
          -Tp(i,:) &
          ) &
          - depsdr_arr(i)*ii*kx*T(i,:) &
          - depsdz_arr(i)*( Tp(i,:) + dlnrhodz_arr(i)*T(i,:) ) &
          - divF_arr(i)*T(i,:) & 
          + eps_arr(i)*(&
          ii*kx*dlnrhodr_arr(i)*T(i,:) + dlnrhodz_arr(i)*( Tp(i,:) + dlnrhodz_arr(i)*T(i,:) ) &
          + kx*kx*T(i,:) &
          - ( Tpp(i,:) + 2d0*Tp(i,:)*dlnrhodz_arr(i) + T(i,:)*(d2lnrhodz2_arr(i) + dlnrhodz_arr(i)**2d0) ) &
          )&
          +0.5d0*Fr_arr(i)*smallh_g*smallq*T(i,:) &
          +0.5d0*eps_arr(i)*smallh_g*smallq*ii*kx*T(i,:) &
          )
    
     L23(i,:) =-gr_at_z*T(i,:) - ii*kx*(1d0 - eps_arr(i))*T(i,:) + (1d0 - eps_arr(i))*smallh_g*smallq*T(i,:)
     L24(i,:) = (0d0, 0d0)
     L25(i,:) =-gz_at_z*T(i,:) -       (1d0 - eps_arr(i))*Tp(i,:)
     
     R21(i,:) = zero
     R22(i,:) = T(i,:)
     R23(i,:) = zero 
     R24(i,:) = zero
     R25(i,:) = zero

     !vx mom eqn 
     L31(i,:) = gr_at_z*T(i,:) !acting on drho/rho
     L32(i,:) = -ii*kx*T(i,:) !- dlnrhodr_arr(i)*T(i,:)!acting on dP/rho 
     L33(i,:) = (0d0, 0d0)
     L34(i,:) = 2d0*omega_arr(i)*T(i,:)
     L35(i,:) = (0d0, 0d0)

     R31(i,:) = zero
     R32(i,:) = zero
     R33(i,:) = T(i,:)
     R34(i,:) = zero
     R35(i,:) = zero
     
     !vy mom eqn 
     L41(i,:) = (0d0, 0d0)
     L42(i,:) = (0d0, 0d0)
     L43(i,:) = -kappa2_arr(i)*T(i,:)/(2d0*omega_arr(i))
     L44(i,:) = 0d0 
     L45(i,:) = -vshear_arr(i)*T(i,:)/(2d0*omega_arr(i))

     R41(i,:) = zero
     R42(i,:) = zero
     R43(i,:) = zero
     R44(i,:) = T(i,:)
     R45(i,:) = zero

     !vz mom eqn 
     L51(i,:) = gz_at_z*T(i,:)
     L52(i,:) = -Tp(i,:) - dlnrhodz_arr(i)*T(i,:)
     L53(i,:) = (0d0, 0d0)
     L54(i,:) = (0d0, 0d0)
     L55(i,:) = (0d0, 0d0)
        
     R51(i,:) = zero
     R52(i,:) = zero
     R53(i,:) = zero
     R54(i,:) = zero
     R55(i,:) = T(i,:)

  enddo 

  !boundary conditions
  !vz  = 0 
  do i = 1, nz, nz-1 
     L51(i,:) = zero 
     L52(i,:) = zero 
     L53(i,:) = zero 
     L54(i,:) = zero 
     L55(i,:) = T(i,:)
     
     R51(i,:) = zero 
     R52(i,:) = zero 
     R53(i,:) = zero 
     R54(i,:) = zero
     R55(i,:) = zero 
  enddo

!!$  !deps = 0 
!!$  do i = 1, nz, nz-1 
!!$     L21(i,:) = (1d0 - eps_arr(i))*T(i,:)
!!$     L22(i,:) =-T(i,:)
!!$     L23(i,:) = zero 
!!$     L24(i,:) = zero 
!!$     L25(i,:) = zero
!!$     
!!$     R21(i,:) = zero 
!!$     R22(i,:) = zero 
!!$     R23(i,:) = zero 
!!$     R24(i,:) = zero
!!$     R25(i,:) = zero 
!!$  enddo

  !DeltaP = 0 
!!$  do i = 1, nz, nz-1 
!!$     L21(i,:) = zero 
!!$     L22(i,:) = zero 
!!$     L23(i,:) = T(i,:)*Fr_arr(i)
!!$     L24(i,:) = zero 
!!$     L25(i,:) = T(i,:)*Fz_arr(i)
!!$     
!!$     R21(i,:) = zero 
!!$     R22(i,:) = T(i,:)
!!$     R23(i,:) = zero 
!!$     R24(i,:) = zero
!!$     R25(i,:) = zero 
!!$  enddo

  !no vertical flux in energy eq
!!$    do i = 1, nz, nz-1 
!!$     L21(i,:) = (1d0 - 2d0*eps_arr(i))*Fz_arr(i)*T(i,:)
!!$     L22(i,:) = -eps_arr(i)*Tp(i,:) - (Fz_arr(i)+eps_arr(i)*dlnrhodz_arr(i))*T(i,:)
!!$     L23(i,:) = zero
!!$     L24(i,:) = zero 
!!$     L25(i,:) = zero
!!$     
!!$     R21(i,:) = zero 
!!$     R22(i,:) = zero
!!$     R23(i,:) = zero 
!!$     R24(i,:) = zero
!!$     R25(i,:) = zero 
!!$  enddo

!note: no lag pert in rhod give crazy results 
  
  !fill big matrices 
  ibeg = 1
  iend = nz
  bigmatrix_lhs(ibeg:iend, 1:nz)          = L11
  bigmatrix_lhs(ibeg:iend, nz+1:2*nz)     = L12
  bigmatrix_lhs(ibeg:iend, 2*nz+1:3*nz )  = L13
  bigmatrix_lhs(ibeg:iend, 3*nz+1:4*nz )  = L14
  bigmatrix_lhs(ibeg:iend, 4*nz+1:bignz)  = L15
  
  bigmatrix_rhs(ibeg:iend, 1:nz)          = R11
  bigmatrix_rhs(ibeg:iend, nz+1:2*nz)     = R12
  bigmatrix_rhs(ibeg:iend, 2*nz+1:3*nz )  = R13
  bigmatrix_rhs(ibeg:iend, 3*nz+1:4*nz )  = R14
  bigmatrix_rhs(ibeg:iend, 4*nz+1:bignz)  = R15

  ibeg = nz+1
  iend = 2*nz
  bigmatrix_lhs(ibeg:iend, 1:nz)          = L21
  bigmatrix_lhs(ibeg:iend, nz+1:2*nz)     = L22
  bigmatrix_lhs(ibeg:iend, 2*nz+1:3*nz )  = L23
  bigmatrix_lhs(ibeg:iend, 3*nz+1:4*nz )  = L24
  bigmatrix_lhs(ibeg:iend, 4*nz+1:bignz)  = L25
  
  bigmatrix_rhs(ibeg:iend, 1:nz)          = R21
  bigmatrix_rhs(ibeg:iend, nz+1:2*nz)     = R22
  bigmatrix_rhs(ibeg:iend, 2*nz+1:3*nz )  = R23
  bigmatrix_rhs(ibeg:iend, 3*nz+1:4*nz )  = R24
  bigmatrix_rhs(ibeg:iend, 4*nz+1:bignz)  = R25

  ibeg = 2*nz+1
  iend = 3*nz
  bigmatrix_lhs(ibeg:iend, 1:nz)          = L31
  bigmatrix_lhs(ibeg:iend, nz+1:2*nz)     = L32
  bigmatrix_lhs(ibeg:iend, 2*nz+1:3*nz )  = L33
  bigmatrix_lhs(ibeg:iend, 3*nz+1:4*nz )  = L34
  bigmatrix_lhs(ibeg:iend, 4*nz+1:bignz)  = L35
  
  bigmatrix_rhs(ibeg:iend, 1:nz)          = R31
  bigmatrix_rhs(ibeg:iend, nz+1:2*nz)     = R32
  bigmatrix_rhs(ibeg:iend, 2*nz+1:3*nz )  = R33
  bigmatrix_rhs(ibeg:iend, 3*nz+1:4*nz )  = R34
  bigmatrix_rhs(ibeg:iend, 4*nz+1:bignz)  = R35

  ibeg = 3*nz+1
  iend = 4*nz
  bigmatrix_lhs(ibeg:iend, 1:nz)          = L41
  bigmatrix_lhs(ibeg:iend, nz+1:2*nz)     = L42
  bigmatrix_lhs(ibeg:iend, 2*nz+1:3*nz )  = L43
  bigmatrix_lhs(ibeg:iend, 3*nz+1:4*nz )  = L44
  bigmatrix_lhs(ibeg:iend, 4*nz+1:bignz)  = L45
  
  bigmatrix_rhs(ibeg:iend, 1:nz)          = R41
  bigmatrix_rhs(ibeg:iend, nz+1:2*nz)     = R42
  bigmatrix_rhs(ibeg:iend, 2*nz+1:3*nz )  = R43
  bigmatrix_rhs(ibeg:iend, 3*nz+1:4*nz )  = R44
  bigmatrix_rhs(ibeg:iend, 4*nz+1:bignz)  = R45

  ibeg = 4*nz+1
  iend = bignz
  bigmatrix_lhs(ibeg:iend, 1:nz)          = L51
  bigmatrix_lhs(ibeg:iend, nz+1:2*nz)     = L52
  bigmatrix_lhs(ibeg:iend, 2*nz+1:3*nz )  = L53
  bigmatrix_lhs(ibeg:iend, 3*nz+1:4*nz )  = L54
  bigmatrix_lhs(ibeg:iend, 4*nz+1:bignz)  = L55
  
  bigmatrix_rhs(ibeg:iend, 1:nz)          = R51
  bigmatrix_rhs(ibeg:iend, nz+1:2*nz)     = R52
  bigmatrix_rhs(ibeg:iend, 2*nz+1:3*nz )  = R53
  bigmatrix_rhs(ibeg:iend, 3*nz+1:4*nz )  = R54
  bigmatrix_rhs(ibeg:iend, 4*nz+1:bignz)  = R55
  
end subroutine fill_bigmatrix

subroutine eigenvalue_problem
  use global
  implicit none
  character*1, parameter :: jobvl = 'N', jobvr='V'
  integer :: i, j, info, loc(1), cnt 
  integer :: lwork 
  complex*16, allocatable :: work(:)
  complex*16 :: apha(bignz), bta(bignz), eigen(bignz), bigQ(nz), dbigQ(nz), d2bigQ(nz), dbigW(nz), dvz(nz)
  complex*16 :: err1(nz), err2(nz), err3(nz), err4(nz), err5(nz)
  complex*16 :: dFx(nz), dFz(nz), ddivF(nz), ep(nz), depdx(nz), depdz(nz), dC(nz), divv(nz)
  complex*16 :: vl(bignz,bignz), vr(bignz,bignz)
  real*8, parameter :: min_rate = 1d-6, max_rate = 1d0 
  real*8 :: rwork(8*bignz), eigen_re, eigen_im, eigen_mag, vz_max, error_in, error_out 
  real*8 :: kx
  
  kx = krad/Hd 

  lwork = 4*bignz
  allocate(work(lwork))

  !eigenvalue problem 
  call zggev (jobvl, jobvr, bignz, bigmatrix_lhs, bignz, bigmatrix_rhs, bignz, apha, bta, VL, & 
       bignz, VR, bignz, WORK, LWORK, RWORK, INFO) 
  if(info .ne. 0 ) print*, 'eigen failed?', info

  !compute eigenvalues, output data 

  open(20,file='eigenvalues.dat')
  open(30,file='eigenvectors.dat')
  open(40,file='error.dat')
  open(50,file='nonadia.dat')
  
  cnt = 0 
  do i=1, bignz
     eigen(i) = apha(i)/bta(i)
     eigen_re = dble(eigen(i))
     eigen_im = dimag(eigen(i))
     eigen_mag= sqrt(eigen_re**2d0 + eigen_im**2d0)     

     !filter out modes that grow too slow or too fast 
     if((eigen_re.ge.min_rate).and.(eigen_mag.lt.max_rate)) then 
        
        bigW  = matmul(T,vr(1:nz,i))
        dbigW = matmul(Tp,vr(1:nz,i))
        bigQ  = matmul(T,vr(nz+1:2*nz,i)) !this is actually dP/rho 
        dbigQ = matmul(Tp,vr(nz+1:2*nz,i))
        d2bigQ= matmul(Tpp,vr(nz+1:2*nz,i))
        vx    = matmul(T,vr(2*nz+1:3*nz,i))
        vy    = matmul(T,vr(3*nz+1:4*nz,i))
        vz    = matmul(T,vr(4*nz+1:bignz,i))
        dvz   = matmul(Tp,vr(4*nz+1:bignz,i))
    
        dfrac  = (1d0 - eps_arr)*bigW - bigQ
 
        cnt = cnt + 1
        write(20,fmt='(2(e22.15,x))'), eigen_re, eigen_im   
        do j=1,nz
           write(30,fmt='(10(e22.15,x))') dble(bigW(j)), dimag(bigW(j)), dble(dfrac(j)), dimag(dfrac(j)), & 
                                          dble(vx(j)), dimag(vx(j)), dble(vy(j)), dimag(vy(j)), dble(vz(j)), dimag(vz(j))
        enddo
 

        !error tests 
        !mass 
        err1 = eigen(i)*bigW + (ii*kx + dlnrhodr_arr)*vx + dlnrhodz_arr*vz + dvz
        err1 = err1/(eigen(i)*bigW(nz))
        !vx 
        dFx  = -bigW*Fr_arr - ii*kx*bigQ 
        err2 = eigen(i)*vx - 2d0*omega_arr*vy - dFx
        err2 = err2/(eigen(i)*vx(nz))
        !vy 
        err3 = eigen(i)*vy + (kappa2_arr*vx + vshear_arr*vz)/(2d0*omega_arr)
        err3 = err3/(eigen(i)*vy(nz))
        !vz
        dFz  = -bigW*Fz_arr - (dbigQ + bigQ*dlnrhodz_arr)
        err4 = eigen(i)*vz - dFz 
        err4 = err4/(eigen(i)*vz(nzmid))
        !energy
        ep   = (1d0 - eps_arr)*bigW - bigQ
        depdx= ((1d0-eps_arr)*(ii*kx-dlnrhodr_arr)-depsdr_arr)*bigW + (smallh_g*smallq + dlnrhodr_arr - ii*kx)*bigQ
        depdz= (1d0 - eps_arr)*dbigW - depsdz_arr*bigW - dbigQ
        ddivF= ((dlnrhodr_arr - ii*kx)*Fr_arr - divF_arr + 1d0)*bigW + (ii*kx*dlnrhodr_arr+kx*kx)*bigQ + eigen(i)*dvz

!        ddivF = Fr_arr*(dlnrhodr_arr - ii*kx)*bigW - Fz_arr*dbigW - bigW*divF_arr + ii*kx*bigQ*dlnrhodr_arr &
!             +dlnrhodz_arr*(dbigQ + bigQ*dlnrhodz_arr) + kx*kx*bigQ &
!             -(d2bigQ + 2d0*dbigQ*dlnrhodz_arr + bigQ*(d2lnrhodz2_arr + dlnrhodz_arr**2d0))

        dC   = Fr_arr*depdx + Fz_arr*depdz + dFx*depsdr_arr + dFz*depsdz_arr + ep*divF_arr + eps_arr*ddivF
        dC   = dC - 0.5d0*smallh_g*smallq*(Fr_arr*ep + eps_arr*dFx)
        dC   = -tstop_arr*dC 
        
        err5 = eigen(i)*bigQ + (1d0-eps_arr)*(ii*kx*vx + dvz) - vx*Fr_arr - vz*Fz_arr - (1d0-eps_arr)*vx*smallh_g*smallq
        err5 = err5 - dC 
        err5 = err5/(eigen(i)*bigQ(nz))

        do j=1,nz
           write(40,fmt='(5(e22.15,x))') abs(err1(j)),  abs(err2(j)),  abs(err3(j)),  abs(err4(j)), abs(err5(j))
        enddo

        !nonadiabatic term
        divv = ii*kx*vx + dvz
        do j=1,nz
           write(50,fmt='(4(e22.15,x))') dble(divv(j)), dimag(divv(j)), dble(dC(j)), dimag(dC(j))
        enddo
        
     endif
  enddo

  print*, 'found', cnt, 'modes' 

  close(20)
  close(30)
  close(40)
  close(50)

end subroutine eigenvalue_problem
