module global_variables
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0,1d0)

  real(8),parameter :: fs = 1d0/0.024189d0
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: angstrom = 1d0/0.52917721067d0
  real(8),parameter :: clight = 137.035999139d0
  

  real(8),parameter :: Ip = 29.24d0*ev
! M.J. Brunger, et al, Phys. Rev. A 59, 1245 (1999)
! Ip of Ar-3s = 29.24 eV
  integer :: nt
  real(8) :: Tprop
  real(8) :: dt_ini
  real(8) :: dt
  real(8) :: tt_i, tt_f

  real(8),allocatable :: A_IR(:)
  real(8),parameter :: T_IRpulse = 60d0*fs
  real(8),parameter :: omega_IRpulse = 1.55d0*ev
  real(8),parameter :: E0_IRpulse = 5.338d-9*sqrt(1d12) ! W/cm2

  real(8),allocatable :: E_EUV(:)
  real(8),parameter :: T_EUVpulse = 30d0*fs
  real(8),parameter :: omega_EUVpulse = 27d0*1.55d0*ev

  integer,parameter :: ndelay = 40
  real(8),parameter :: Tdelay_i  = -30d0*fs
  real(8),parameter :: Tdelay_f = 30d0*fs
  real(8) :: Tdelay

  integer,parameter :: NEkin = 200
  real(8),parameter :: Ekin_i = 5d0*ev
  real(8),parameter :: Ekin_f = 20d0*ev
  real(8),parameter :: dEkin = (Ekin_f - Ekin_i)/NEkin
  real(8) :: PES_yield(0:NEkin)

  integer,parameter :: Ncos = 10

end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none


  call set_parameters

  call calc_pes

end program main
!-------------------------------------------------------------------------------
subroutine set_parameters
  use global_variables
  implicit none

  dt_ini = 0.32d0


end subroutine set_parameters
!-------------------------------------------------------------------------------
subroutine calc_pes
  use global_variables
  implicit none
  integer :: it_delay, iekin, icos
  real(8) :: Ekin, cos_theta
  real(8) :: pop

  open(20,file='td_pes.out')
  write(20,"(A)")"# Tdelay (au), Ekin (au), PES_yield (arb. unit)"
  do it_delay = 0, ndelay
    write(*,"(A,I7,2x,A,2x,I7)")"it_delay=",it_delay,'/',ndelay
    Tdelay = Tdelay_i + it_delay*(Tdelay_f - Tdelay_i)/dble(ndelay)
    call laser_init

    PES_yield = 0d0

    do iekin = 0, NEkin
      Ekin = Ekin_i + iekin*(Ekin_f - Ekin_i)/NEkin

      do icos = 0, Ncos
        cos_theta = -1d0 + 2d0*icos/dble(Ncos)
        call calc_pes_population(Ekin, cos_theta, pop)
        PES_yield(iekin) = PES_yield(iekin) + pop
      end do

    end do

    call laser_fin

    do iekin = 0, NEkin
      Ekin = Ekin_i + iekin*(Ekin_f - Ekin_i)/NEkin
      write(20,"(999e26.16e3)")Tdelay,Ekin,PES_yield(iekin)
    end do
    write(20,*)
  end do
    


  close(20)

end subroutine calc_pes
!-------------------------------------------------------------------------------
subroutine laser_init
  use global_variables
  implicit none
  integer :: it
  real(8) :: xx

  tt_i = min(-0.5d0*T_IRpulse, -0.5d0*T_EUVpulse + 0.5d0*T_IRpulse + Tdelay)
  tt_f = max(0.5d0*T_IRpulse, 0.5d0*T_EUVpulse + 0.5d0*T_IRpulse + Tdelay)
  write(*,*)"tt_i =",tt_i/fs
  write(*,*)"tt_f =",tt_f/fs
  Tprop = tt_f - tt_i
  nt = aint(Tprop/dt_ini) + 1
  dt = Tprop/nt

  allocate(A_IR(0:nt), E_EUV(0:nt))
  A_IR = 0d0; E_EUV = 0d0
  
! IR field  
  do it = 0, nt
    xx = tt_i + dt*it - 0.5d0*T_IRpulse
    if(abs(xx) < 0.5d0*T_IRpulse)then
      A_IR(it) = -E0_IRpulse/omega_IRpulse*sin(omega_IRpulse*xx) &
        *cos(pi*xx/T_IRpulse)**4
    end if
  end do

! EUV field  
  do it = 0, nt
    xx = tt_i + dt*it - 0.5d0*T_IRpulse - Tdelay
    if(abs(xx) < 0.5d0*T_EUVpulse)then
      E_EUV(it) = sin(omega_EUVpulse*xx)*cos(pi*xx/T_EUVpulse)**4
    end if
  end do

!  open(101,file='EUV.out')
!  do it = 0,nt
!    xx = tt_i + dt*it
!    write(101,"(999e26.16e3)")xx,E_EUV(it),A_IR(it)
!  end do
!  close(101)

end subroutine laser_init
!-------------------------------------------------------------------------------
subroutine laser_fin
  use global_variables
  implicit none

  deallocate(A_IR, E_EUV)

end subroutine laser_fin
!-------------------------------------------------------------------------------
subroutine calc_pes_population(Ekin, cos_theta, pop)
  use global_variables
  implicit none
  real(8),intent(in)  :: Ekin, cos_theta
  real(8),intent(out) :: pop
  real(8) :: Ekin_t, phi
  real(8) :: pz0, pz, dz, pxy2
  complex(8) :: zb
  integer :: it

  pz0 = sqrt(2d0*Ekin)*cos_theta
  pxy2 = 2d0*Ekin*(1d0-cos_theta**2)

  zb = 0d0
  phi = 0d0
  do it = 0, nt
    pz = pz0 + A_IR(it)
    Ekin_t = 0.5d0*(pz**2+pxy2)
    dz = pz/(Ekin_t + Ip)**3

    phi = phi + 0.5d0*dt*(Ekin_t + Ip)

    zb = zb + zI*E_EUV(it)*dz*exp(-zI*phi)

    phi = phi + 0.5d0*dt*(Ekin_t + Ip)

  end do
  zb = zb*dt

  pop = abs(zb)**2



end subroutine calc_pes_population
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
