module load_Bethell_Xray_cross

! Bethell 2011, Table 2

implicit none

integer, parameter, private :: nrow = 16

double precision, dimension(2, nrow), private :: &
  E_r = reshape((/ &
    0.030D0,  0.055D0, &
    0.055D0,  0.100D0, &
    0.100D0,  0.165D0, &
    0.165D0,  0.284D0, &
    0.284D0,  0.400D0, &
    0.400D0,  0.532D0, &
    0.532D0,  0.708D0, &
    0.708D0,  0.867D0, &
    0.867D0,  1.303D0, &
    1.303D0,  1.840D0, &
    1.840D0,  2.471D0, &
    2.471D0,  3.210D0, &
    3.210D0,  4.038D0, &
    4.038D0,  7.111D0, &
    7.111D0,  8.331D0, &
    8.331D0,  10.00D0   /), (/2, nrow/))

double precision, dimension(3, nrow), private :: &
  c_g = reshape((/ &
    14.2D0,   727D0,    -4130D0, &
    22D0,     445D0,    -1550D0, &
    31D0,     263D0,    -614D0,  &
    43.7D0,   112D0,    -165D0,  &
    49D0,     86D0,     -103D0,  &
    58.6D0,   36.9D0,   -39.9D0, &
    48D0,     130D0,    -82.2D0, &
    77.4D0,   46.3D0,   -22D0,   &
    80.1D0,   69.8D0,   -28.3D0, &
    117D0,    7.43D0,   -1.87D0, &
    107D0,    16D0,     -3.75D0, &
    106D0,    13.6D0,   -2.63D0, &
    138D0,    -1.99D0,  -0.179D0,&
    142D0,    -4.7D0,   0.239D0, &
    138D0,    -3.36D0,  0.133D0, &
    88.9D0,   8.15D0,   -0.547D0   /), (/3, nrow/))

double precision, dimension(3, nrow), private :: &
  c_d = reshape((/ &
    0.0344D0, -1.62D0,  88.2D0,  &
    -0.147D0, 4.19D0,   48.1D0,  &
    -0.677D0, 14.9D0,   9.6D0,   &
    -1.12D0,  23.6D0,   -16.2D0, &
    0.188D0,  24.6D0,   -1.09D0, &
    -3.57D0,  55.5D0,   -37.9D0, &
    -8.24D0,  89.6D0,   -48.1D0, &
    57.1D0,   -49.9D0,  52.1D0,  &
    9.11D0,   72.7D0,   -20.8D0, &
    -8.71D0,  106D0,    -25.7D0, &
    34.9D0,   72.4D0,   -11.4D0, &
    23.6D0,   85.1D0,   -11.3D0, &
    116D0,    28.2D0,   -2.55D0, &
    191D0,    -2.92D0,  1.09D0,  &
    812D0,    -74.7D0,  6.49D0,  &
    -33D0,    137D0,    -6.39D0   /), (/3, nrow/))


contains


pure function sigma_Xray_Bethell(E, eps, G, a)
  use phy_const
  double precision sigma_Xray_Bethell
  double precision, intent(in) :: E, eps, G, a
  double precision sigma_dust_per_H, sigma_gas_per_H, tau, f
  integer i, i0
  if ((E .ge. E_r(1, 1)) .and. (E .le. E_r(2, nrow))) then
    do i=1, nrow
      if ((E .ge. E_r(1, i)) .and. (E .le. E_r(2, i))) then
        i0 = i
        exit
      end if
    end do
  else if (E .lt. E_r(1, 1)) then
    i0 = 1
  else
    i0 = nrow
  end if
  sigma_dust_per_H = 1D-24 / (E*E*E) * (c_d(1,i0) + (c_d(2,i0) + c_d(3,i0) * E) * E) * eps
  sigma_gas_per_H  = 1D-24 / (E*E*E) * (c_g(1,i0) + (c_g(2,i0) + c_g(3,i0) * E) * E)
  if ((eps .le. 1D-30) .or. (G .le. 1D-30)) then
    f = 1D0
  else
    tau = sigma_dust_per_H / G * (3D0/phy_2Pi) / (a*a)
    f = 1.5D0 / tau * (1D0 - 2D0/tau/tau * (1D0 - (tau+1D0)*exp(-tau)))
  end if
  sigma_Xray_Bethell = sigma_gas_per_H + f * sigma_dust_per_H
end function sigma_Xray_Bethell


pure function sigma_Xray_Bethell_gas(E)
  use phy_const
  double precision sigma_Xray_Bethell_gas
  double precision, intent(in) :: E
  integer i, i0
  if ((E .ge. E_r(1, 1)) .and. (E .le. E_r(2, nrow))) then
    do i=1, nrow
      if ((E .ge. E_r(1, i)) .and. (E .le. E_r(2, i))) then
        i0 = i
        exit
      end if
    end do
  else if (E .lt. E_r(1, 1)) then
    i0 = 1
  else
    i0 = nrow
  end if
  sigma_Xray_Bethell_gas =  &
    1D-24 / (E*E*E) * (c_g(1,i0) + (c_g(2,i0) + c_g(3,i0) * E) * E)
end function sigma_Xray_Bethell_gas


pure function sigma_Xray_Bethell_dust(E, eps, G, a)
  use phy_const
  double precision sigma_Xray_Bethell_dust
  double precision, intent(in) :: E, eps, G, a
  double precision tau, f
  integer i, i0
  if ((G .le. 0D0) .or. (a .le. 0D0) .or. (eps .le. 0D0)) then
    sigma_Xray_Bethell_dust = 0D0
    return
  end if
  if ((E .ge. E_r(1, 1)) .and. (E .le. E_r(2, nrow))) then
    do i=1, nrow
      if ((E .ge. E_r(1, i)) .and. (E .le. E_r(2, i))) then
        i0 = i
        exit
      end if
    end do
  else if (E .lt. E_r(1, 1)) then
    i0 = 1
  else
    i0 = nrow
  end if
  sigma_Xray_Bethell_dust = 1D-24 / (E*E*E) * (c_d(1,i0) + (c_d(2,i0) + c_d(3,i0) * E) * E) * eps
  tau = sigma_Xray_Bethell_dust / G * (3D0/phy_2Pi) / (a*a)
  f = 1.5D0 / tau * (1D0 - 2D0/tau/tau * (1D0 - (tau+1D0)*exp(-tau)))
  sigma_Xray_Bethell_dust = f * sigma_Xray_Bethell_dust
end function sigma_Xray_Bethell_dust


end module load_Bethell_Xray_cross
