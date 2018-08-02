module global_variables
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0,1d0)

  real(8),parameter :: fs = 1d0/0.024189d0
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: angstrom = 1d0/0.52917721067d0
  real(8),parameter :: clight = 137.035999139d0
  

  real(8),parameter :: Ip = 21d0*ev

  real(8) :: Tprop
  real(8) :: dt_ini
  real(8) :: dt
  real(8) :: tt_i, tt_f

  real(8),parameter :: T_IRpulse = 30d0*fs
  real(8),parameter :: omega_IRpulse = 1.55d0*ev
  real(8),parameter :: E0_IRpulse = 5.338d-9*sqrt(1d12) ! W/cm2

  real(8),parameter :: T_EUVpulse = 30d0*fs
  real(8),parameter :: omega_EUVpulse = 27d0*1.55d0*ev

  integer,parameter :: ndelay = 120
  real(8),parameter :: Tdelay_i  = -30d0*fs
  real(8),parameter :: T_delay_f = 30d0*fs

  integer,parameter :: Ne = 200
  real(8),parameter :: Ekin_i = 10d0*ev
  real(8),parameter :: Ekin_f = 40d0*ev
  real(8),parameter :: dEkin = (Ekin_f - Ekin_i)/Ne
  real(8) :: PES_yield(0:ne)

end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none


  call set_parameters


end program main
!-------------------------------------------------------------------------------
subroutine set_parameters
  use global_variables
  implicit none

  dt_ini = 0.08d0


end subroutine set_parameters
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
