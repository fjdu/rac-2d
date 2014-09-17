module load_Neufeld_cooling_CO

use trivials

implicit none

integer, parameter, private :: &
  n_T_high = 6, &
  n_T_high_vib = 6, &
  n_T_low = 6, &
  n_log10N_high = 10, &
  n_log10N_high_vib = 8, &
  n_log10N_low = 10

type, private :: Neufeld_cooling_CO_params
  double precision, dimension(n_T_high), private :: &
    T_high     = (/100D0,   300D0,   600D0,   1000D0,  1500D0,  2000D0/)
  double precision, dimension(n_T_high_vib), private :: &
    T_high_vib = (/100D0,   200D0,   400D0,   1000D0,  2000D0,  4000D0/)
  double precision, dimension(n_T_low), private :: &
    T_low  = (/10D0,    20D0,    30D0,    50D0,    80D0,    100D0/)

  double precision, dimension(n_log10N_high), private :: &
    log10N_high = (/14.5D0, 15.0D0, 15.5D0, 16.0D0, 16.5D0, 17.0D0, 17.5D0, 18.0D0, 18.5D0, 19.0D0/)
  double precision, dimension(n_log10N_high_vib), private :: &
    log10N_high_vib = (/13.D0,  14.D0,  15.D0,  16.D0,  17.D0,  18.D0,  19.D0,  20.0D0/)
  double precision, dimension(n_log10N_low), private :: &
    log10N_low  = (/14.5, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0/)

  double precision, dimension(n_T_high), private :: &
    log10_L0_high      = (/23.88, 23.40, 23.07, 22.81, 22.61, 22.47/)
  double precision, dimension(n_T_low), private :: &
    log10_L0_low       = (/24.96, 24.56, 24.36, 24.12, 23.88, 23.80/)

  double precision, dimension(n_T_high, n_log10N_high), private :: &
    log10_L_LTE_high = reshape((/ &
        18.81, 17.82, 17.23, 16.82, 16.52, 16.33, &
        18.81, 17.82, 17.23, 16.82, 16.52, 16.33, &
        18.83, 17.83, 17.24, 16.82, 16.52, 16.33, &
        18.87, 17.84, 17.24, 16.83, 16.52, 16.33, &
        18.96, 17.88, 17.26, 16.84, 16.53, 16.33, &
        19.14, 17.98, 17.32, 16.88, 16.56, 16.35, &
        19.40, 18.15, 17.44, 16.96, 16.62, 16.40, &
        19.71, 18.41, 17.65, 17.13, 16.76, 16.52, &
        20.07, 18.72, 17.93, 17.38, 16.98, 16.72, &
        20.45, 19.07, 18.26, 17.69, 17.27, 17.00, &
        20.85, 19.45, 18.62, 18.03, 17.61, 17.33  &
      /), (/n_T_high, n_log10N_high/))

  double precision, dimension(n_T_high_vib, n_log10N_high_vib), private :: &
    log10_X_L_LTE_high_vib = reshape((/ &
        10.83D0, 10.82D0, 10.82D0, 10.80D0, 10.77D0, 10.79D0, &
        10.84D0, 10.83D0, 10.82D0, 10.80D0, 10.77D0, 10.79D0, &
        10.88D0, 10.86D0, 10.84D0, 10.81D0, 10.77D0, 10.80D0, &
        11.20D0, 11.11D0, 11.03D0, 10.92D0, 10.82D0, 10.82D0, &
        11.99D0, 11.86D0, 11.73D0, 11.42D0, 11.17D0, 11.04D0, &
        12.89D0, 12.77D0, 12.61D0, 12.18D0, 11.95D0, 11.75D0, &
        13.82D0, 13.70D0, 13.44D0, 13.03D0, 12.63D0, 12.31D0, &
        14.76D0, 14.65D0, 14.30D0, 13.93D0, 13.46D0, 13.09D0  &
      /), (/n_T_high_vib, n_log10N_high_vib/))

  double precision, dimension(n_T_low, n_log10N_low), private :: &
    log10_L_LTE_low = reshape((/ &
        21.15, 20.39, 19.97, 19.47, 19.02, 18.81, &
        21.25, 20.45, 20.02, 19.50, 19.04, 18.83, &
        21.45, 20.58, 20.12, 19.57, 19.09, 18.87, &
        21.73, 20.81, 20.31, 19.71, 19.20, 18.96, &
        22.08, 21.11, 20.58, 19.95, 19.39, 19.14, &
        22.47, 21.46, 20.91, 20.25, 19.67, 19.40, &
        22.89, 21.85, 21.28, 20.60, 19.99, 19.71, &
        23.32, 22.25, 21.67, 20.97, 20.35, 20.07, &
        23.76, 22.68, 22.08, 21.37, 20.74, 20.45, &
        24.21, 23.11, 22.51, 21.78, 21.14, 20.85  &
      /), (/n_T_low, n_log10N_low/))

  double precision, dimension(n_T_high, n_log10N_high), private :: &
    log10_n_12_high = reshape((/ &
         4.43,  4.95,  5.25,  5.43,  5.57,  5.68, &
         4.42,  4.95,  5.24,  5.43,  5.57,  5.68, &
         4.38,  4.93,  5.23,  5.42,  5.57,  5.68, &
         4.24,  4.88,  5.20,  5.42,  5.56,  5.68, &
         4.01,  4.75,  5.13,  5.36,  5.54,  5.65, &
         3.64,  4.50,  4.95,  5.24,  5.44,  5.58, &
         3.21,  4.12,  4.65,  5.00,  5.25,  5.41, &
         2.72,  3.67,  4.25,  4.63,  4.92,  5.13, &
         2.22,  3.20,  3.77,  4.18,  4.51,  4.73, &
         1.73,  2.69,  3.29,  3.70,  4.03,  4.26, &
         1.23,  2.20,  2.79,  3.20,  3.54,  3.78  &
      /), (/n_T_high, n_log10N_high/))

  double precision, dimension(n_T_low, n_log10N_low), private :: &
    log10_n_12_low = reshape((/ &
         3.31,  3.72,  3.96,  4.24,  4.39,  4.51, &
         3.11,  3.59,  3.87,  4.17,  4.35,  4.47, &
         2.80,  3.34,  3.64,  4.01,  4.23,  4.36, &
         2.37,  2.95,  3.31,  3.73,  3.99,  4.13, &
         1.90,  2.50,  2.87,  3.33,  3.61,  3.78, &
         1.40,  2.02,  2.40,  2.86,  3.17,  3.34, &
         0.90,  1.52,  1.90,  2.38,  2.69,  2.87, &
         0.40,  1.02,  1.41,  1.89,  2.19,  2.38, &
        -0.09,  0.52,  0.91,  1.39,  1.69,  1.88, &
        -0.59,  0.03,  0.41,  0.88,  1.19,  1.38  &
      /), (/n_T_low, n_log10N_low/))

  double precision, dimension(n_T_high, n_log10N_high), private :: &
    alpha_high = reshape((/ &
         0.39,  0.38,  0.38,  0.37,  0.36,  0.35, &
         0.38,  0.38,  0.38,  0.37,  0.36,  0.35, &
         0.37,  0.37,  0.37,  0.37,  0.36,  0.35, &
         0.36,  0.36,  0.36,  0.36,  0.35,  0.34, &
         0.36,  0.35,  0.35,  0.35,  0.34,  0.33, &
         0.38,  0.35,  0.34,  0.34,  0.33,  0.32, &
         0.41,  0.37,  0.35,  0.33,  0.32,  0.31, &
         0.44,  0.40,  0.37,  0.35,  0.32,  0.31, &
         0.47,  0.43,  0.40,  0.38,  0.35,  0.33, &
         0.49,  0.46,  0.43,  0.41,  0.38,  0.36, &
         0.51,  0.48,  0.46,  0.44,  0.41,  0.39  &
      /), (/n_T_high, n_log10N_high/))

  double precision, dimension(n_T_low, n_log10N_low), private :: &
    alpha_low = reshape((/ &
         0.37,  0.41,  0.42,  0.41,  0.38,  0.37, &
         0.38,  0.42,  0.44,  0.42,  0.39,  0.37, &
         0.38,  0.43,  0.44,  0.43,  0.40,  0.39, &
         0.36,  0.41,  0.43,  0.43,  0.40,  0.39, &
         0.34,  0.37,  0.39,  0.39,  0.39,  0.39, &
         0.31,  0.33,  0.35,  0.35,  0.35,  0.34, &
         0.28,  0.31,  0.32,  0.32,  0.31,  0.31, &
         0.26,  0.28,  0.29,  0.28,  0.28,  0.28, &
         0.25,  0.26,  0.26,  0.26,  0.26,  0.26, &
         0.23,  0.25,  0.25,  0.24,  0.25,  0.25  &
      /), (/n_T_low, n_log10N_low/))

  ! The user should provide these five quantities.
  ! G is of order 1.
  ! dv/dz is expressed in km s-1 cm-1.
  ! n is in cm-3.
  ! So the N parameter has dimension cm-2 km s-1.
  double precision T, log10N

  double precision L, L_LTE, L0, n_12, alpha
  double precision L_LTE_vib, L0_vib

end type Neufeld_cooling_CO_params

type(Neufeld_cooling_CO_params) a_Neufeld_cooling_CO_params

double precision, parameter, private :: ln10 = log(10D0)


contains


!function cooling_Neufeld_CO()
!  double precision cooling_Neufeld_CO
!  associate( &
!    L     => a_Neufeld_cooling_CO_params%L, &
!    L0    => a_Neufeld_cooling_CO_params%L0, &
!    L_LTE => a_Neufeld_cooling_CO_params%L_LTE, &
!    n_12  => a_Neufeld_cooling_CO_params%n_12, &
!    alpha => a_Neufeld_cooling_CO_params%alpha, &
!    n_H2  => a_Neufeld_cooling_CO_params%n_H2, &
!    X_M   => a_Neufeld_cooling_CO_params%X_M, &
!    G     => a_Neufeld_cooling_CO_params%G, &
!    dv_dz => a_Neufeld_cooling_CO_params%dv_dz, &
!    log10N=> a_Neufeld_cooling_CO_params%log10N)
!    log10N = log10(G * X_M * n_H2 / dv_dz)
!    L0    = get_L0()
!    L_LTE = get_L_LTE()
!    n_12  = get_n_12()
!    alpha = get_alpha()
!    L = 1D0 / &
!      (1D0/L0 + n_H2/L_LTE + &
!       1D0/L0 * (n_H2/n_12)**alpha * (1D0 - n_12*L0/L_LTE))
!    cooling_Neufeld_CO = L * n_H2 * n_H2 * X_M
!  end associate
!end function cooling_Neufeld_CO


function get_L0()
  double precision get_L0
  integer i
  associate( &
    a => a_Neufeld_cooling_CO_params, &
    T => a_Neufeld_cooling_CO_params%T)
    if (T .GE. a%T_high(1)) then
      associate( &
        x => T, &
        y => a%T_high, &
        n => n_T_high, &
        idx => i)
        if      (y(1) .GE. x) then
          idx = 2
        else if (y(n) .LE. x) then
          idx = n
        else
          do idx=2, n
            if (y(idx) .GT. x) then
              exit
            end if
          end do
          if (idx .GT. n) then
            idx = n
          end if
        end if
      end associate
      associate( &
        k  => (a%log10_L0_high(i) - a%log10_L0_high(i-1)) &
            / (log(a%T_high(i)) - log(a%T_high(i-1))), &
        dx => log(T) - log(a%T_high(i-1)), &
        y0 => a%log10_L0_high(i-1))
        get_L0 = k * dx + y0
      end associate
    else
      associate( &
        x => T, &
        y => a%T_low, &
        n => n_T_low, &
        idx => i)
        if      (y(1) .GE. x) then
          idx = 2
        else if (y(n) .LE. x) then
          idx = n
        else
          do idx=2, n
            if (y(idx) .GT. x) then
              exit
            end if
          end do
          if (idx .GT. n) then
            idx = n
          end if
        end if
      end associate
      associate( &
        k  => (a%log10_L0_low(i) - a%log10_L0_low(i-1)) &
             / (a%T_low(i) - a%T_low(i-1)), &
        dx => T - a%T_low(i-1), &
        y0 => a%log10_L0_low(i-1))
        get_L0 = dx * k + y0
      end associate
    end if
  end associate
  get_L0 = exp(-get_L0 * ln10)
end function get_L0


function get_L_LTE()
  double precision get_L_LTE
  integer i, j
  associate( &
    a => a_Neufeld_cooling_CO_params, &
    T => a_Neufeld_cooling_CO_params%T, &
    log10N => a_Neufeld_cooling_CO_params%log10N)
    if (T .GE. 100D0) then
      associate( &
        x => T, &
        y => a%T_high, &
        n => n_T_high, &
        idx => i)
        if      (y(1) .GE. x) then
          idx = 2
        else if (y(n) .LE. x) then
          idx = n
        else
          do idx=2, n
            if (y(idx) .GT. x) then
              exit
            end if
          end do
          if (idx .GT. n) then
            idx = n
          end if
        end if
      end associate
      associate( &
        x => log10N, &
        y => a%log10N_high, &
        n => n_log10N_high, &
        idx => j)
        if      (y(1) .GE. x) then
          idx = 2
        else if (y(n) .LE. x) then
          idx = n
        else
          do idx=2, n
            if (y(idx) .GT. x) then
              exit
            end if
          end do
          if (idx .GT. n) then
            idx = n
          end if
        end if
      end associate
      associate( &
        x   => log(T), &
        y   => log10N, &
        x1  => log(a%T_high(i-1)), &
        x2  => log(a%T_high(i)), &
        y1  => a%log10N_high(j-1), &
        y2  => a%log10N_high(j), &
        z11 => a%log10_L_LTE_high(i-1,j-1), &
        z12 => a%log10_L_LTE_high(i-1,j), &
        z21 => a%log10_L_LTE_high(i,j-1), &
        z22 => a%log10_L_LTE_high(i,j), &
        z => get_L_LTE)
        z = calc_four_point_linear_interpol(x, y, x1, x2, y1, y2, z11, z12, z21, z22)
        !associate (&
        !  k1 => (z12-z11) / (y2-y1), &
        !  k2 => (z22-z21) / (y2-y1), &
        !  z0_1 => z11, &
        !  z0_2 => z21)
        !  associate (&
        !    k_k  => (k2-k1)/(x2-x1), &
        !    k_0  => k1, &
        !    k_z0 => (z0_2-z0_1)/(x2-x1), &
        !    dx   => x-x1, &
        !    dy   => y-y1)
        !    z    = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
        !  end associate
        !end associate
      end associate
    else
      associate( &
        x => T, &
        y => a%T_low, &
        n => n_T_low, &
        idx => i)
        if      (y(1) .GE. x) then
          idx = 2
        else if (y(n) .LE. x) then
          idx = n
        else
          do idx=2, n
            if (y(idx) .GT. x) then
              exit
            end if
          end do
          if (idx .GT. n) then
            idx = n
          end if
        end if
      end associate
      associate( &
        x => log10N, &
        y => a%log10N_low, &
        n => n_log10N_low, &
        idx => j)
        if      (y(1) .GE. x) then
          idx = 2
        else if (y(n) .LE. x) then
          idx = n
        else
          do idx=2, n
            if (y(idx) .GT. x) then
              exit
            end if
          end do
          if (idx .GT. n) then
            idx = n
          end if
        end if
      end associate
      associate( &
        x   => log(T), &
        y   => log10N, &
        x1  => log(a%T_low(i-1)), &
        x2  => log(a%T_low(i)), &
        y1  => a%log10N_low(j-1), &
        y2  => a%log10N_low(j), &
        z11 => a%log10_L_LTE_low(i-1,j-1), &
        z12 => a%log10_L_LTE_low(i-1,j), &
        z21 => a%log10_L_LTE_low(i,j-1), &
        z22 => a%log10_L_LTE_low(i,j), &
        z   => get_L_LTE)
        z = calc_four_point_linear_interpol(x, y, x1, x2, y1, y2, z11, z12, z21, z22)
        !associate (&
        !  k1 => (z12-z11) / (y2-y1), &
        !  k2 => (z22-z21) / (y2-y1), &
        !  z0_1 => z11, &
        !  z0_2 => z21)
        !  associate (&
        !    k_k  => (k2-k1)/(x2-x1), &
        !    k_0  => k1, &
        !    k_z0 => (z0_2-z0_1)/(x2-x1), &
        !    dx   => x-x1, &
        !    dy   => y-y1)
        !    z    = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
        !  end associate
        !end associate
      end associate
    end if
  end associate
  get_L_LTE = exp(-get_L_LTE * ln10)
end function get_L_LTE


function get_n_12()
  double precision get_n_12
  integer i, j
  associate( &
    a => a_Neufeld_cooling_CO_params, &
    T => a_Neufeld_cooling_CO_params%T, &
    log10N => a_Neufeld_cooling_CO_params%log10N)
    if (T .GE. 100D0) then
      associate( &
        x => T, &
        y => a%T_high, &
        n => n_T_high, &
        idx => i)
        if      (y(1) .GE. x) then
          idx = 2
        else if (y(n) .LE. x) then
          idx = n
        else
          do idx=2, n
            if (y(idx) .GT. x) then
              exit
            end if
          end do
          if (idx .GT. n) then
            idx = n
          end if
        end if
      end associate
      associate( &
        x => log10N, &
        y => a%log10N_high, &
        n => n_log10N_high, &
        idx => j)
        if      (y(1) .GE. x) then
          idx = 2
        else if (y(n) .LE. x) then
          idx = n
        else
          do idx=2, n
            if (y(idx) .GT. x) then
              exit
            end if
          end do
          if (idx .GT. n) then
            idx = n
          end if
        end if
      end associate
      associate( &
        x   => log(T), &
        y   => log10N, &
        x1  => log(a%T_high(i-1)), &
        x2  => log(a%T_high(i)), &
        y1  => a%log10N_high(j-1), &
        y2  => a%log10N_high(j), &
        z11 => a%log10_n_12_high(i-1,j-1), &
        z12 => a%log10_n_12_high(i-1,j), &
        z21 => a%log10_n_12_high(i,j-1), &
        z22 => a%log10_n_12_high(i,j), &
        z => get_n_12)
        z = calc_four_point_linear_interpol(x, y, x1, x2, y1, y2, z11, z12, z21, z22)
        !associate (&
        !  k1 => (z12-z11) / (y2-y1), &
        !  k2 => (z22-z21) / (y2-y1), &
        !  z0_1 => z11, &
        !  z0_2 => z21)
        !  associate (&
        !    k_k  => (k2-k1)/(x2-x1), &
        !    k_0  => k1, &
        !    k_z0 => (z0_2-z0_1)/(x2-x1), &
        !    dx   => x-x1, &
        !    dy   => y-y1)
        !    z    = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
        !  end associate
        !end associate
      end associate
    else
      associate( &
        x => T, &
        y => a%T_low, &
        n => n_T_low, &
        idx => i)
        if      (y(1) .GE. x) then
          idx = 2
        else if (y(n) .LE. x) then
          idx = n
        else
          do idx=2, n
            if (y(idx) .GT. x) then
              exit
            end if
          end do
          if (idx .GT. n) then
            idx = n
          end if
        end if
      end associate
      associate( &
        x => log10N, &
        y => a%log10N_low, &
        n => n_log10N_low, &
        idx => j)
        if      (y(1) .GE. x) then
          idx = 2
        else if (y(n) .LE. x) then
          idx = n
        else
          do idx=2, n
            if (y(idx) .GT. x) then
              exit
            end if
          end do
          if (idx .GT. n) then
            idx = n
          end if
        end if
      end associate
      associate( &
        x   => log(T), &
        y   => log10N, &
        x1  => log(a%T_low(i-1)), &
        x2  => log(a%T_low(i)), &
        y1  => a%log10N_low(j-1), &
        y2  => a%log10N_low(j), &
        z11 => a%log10_n_12_low(i-1,j-1), &
        z12 => a%log10_n_12_low(i-1,j), &
        z21 => a%log10_n_12_low(i,j-1), &
        z22 => a%log10_n_12_low(i,j), &
        z   => get_n_12)
        z = calc_four_point_linear_interpol(x, y, x1, x2, y1, y2, z11, z12, z21, z22)
        !associate (&
        !  k1 => (z12-z11) / (y2-y1), &
        !  k2 => (z22-z21) / (y2-y1), &
        !  z0_1 => z11, &
        !  z0_2 => z21)
        !  associate (&
        !    k_k  => (k2-k1)/(x2-x1), &
        !    k_0  => k1, &
        !    k_z0 => (z0_2-z0_1)/(x2-x1), &
        !    dx   => x-x1, &
        !    dy   => y-y1)
        !    z    = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
        !  end associate
        !end associate
      end associate
    end if
  end associate
  get_n_12 = exp(-get_n_12 * ln10)
end function get_n_12


function get_alpha()
  double precision get_alpha
  integer i, j
  associate( &
    a => a_Neufeld_cooling_CO_params, &
    T => a_Neufeld_cooling_CO_params%T, &
    log10N => a_Neufeld_cooling_CO_params%log10N)
    if (T .GE. 100D0) then
      associate( &
        x => T, &
        y => a%T_high, &
        n => n_T_high, &
        idx => i)
        if      (y(1) .GE. x) then
          idx = 2
        else if (y(n) .LE. x) then
          idx = n
        else
          do idx=2, n
            if (y(idx) .GT. x) then
              exit
            end if
          end do
          if (idx .GT. n) then
            idx = n
          end if
        end if
      end associate
      associate( &
        x => log10N, &
        y => a%log10N_high, &
        n => n_log10N_high, &
        idx => j)
        if      (y(1) .GE. x) then
          idx = 2
        else if (y(n) .LE. x) then
          idx = n
        else
          do idx=2, n
            if (y(idx) .GT. x) then
              exit
            end if
          end do
          if (idx .GT. n) then
            idx = n
          end if
        end if
      end associate
      associate( &
        x   => log(T), &
        y   => log10N, &
        x1  => log(a%T_high(i-1)), &
        x2  => log(a%T_high(i)), &
        y1  => a%log10N_high(j-1), &
        y2  => a%log10N_high(j), &
        z11 => a%alpha_high(i-1,j-1), &
        z12 => a%alpha_high(i-1,j), &
        z21 => a%alpha_high(i,j-1), &
        z22 => a%alpha_high(i,j))
        get_alpha = calc_four_point_linear_interpol(x, y, x1, x2, y1, y2, z11, z12, z21, z22)
        !associate (&
        !  k1 => (z12-z11) / (y2-y1), &
        !  k2 => (z22-z21) / (y2-y1), &
        !  z0_1 => z11, &
        !  z0_2 => z21)
        !  associate (&
        !    k_k  => (k2-k1)/(x2-x1), &
        !    k_0  => k1, &
        !    k_z0 => (z0_2-z0_1)/(x2-x1), &
        !    dx   => x-x1, &
        !    dy   => y-y1)
        !    z    = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
        !  end associate
        !end associate
      end associate
    else
      associate( &
        x => T, &
        y => a%T_low, &
        n => n_T_low, &
        idx => i)
        if      (y(1) .GE. x) then
          idx = 2
        else if (y(n) .LE. x) then
          idx = n
        else
          do idx=2, n
            if (y(idx) .GT. x) then
              exit
            end if
          end do
          if (idx .GT. n) then
            idx = n
          end if
        end if
      end associate
      associate( &
        x => log10N, &
        y => a%log10N_low, &
        n => n_log10N_low, &
        idx => j)
        if      (y(1) .GE. x) then
          idx = 2
        else if (y(n) .LE. x) then
          idx = n
        else
          do idx=2, n
            if (y(idx) .GT. x) then
              exit
            end if
          end do
          if (idx .GT. n) then
            idx = n
          end if
        end if
      end associate
      associate( &
        x   => log(T), &
        y   => log10N, &
        x1  => log(a%T_low(i-1)), &
        x2  => log(a%T_low(i)), &
        y1  => a%log10N_low(j-1), &
        y2  => a%log10N_low(j), &
        z11 => a%alpha_low(i-1,j-1), &
        z12 => a%alpha_low(i-1,j), &
        z21 => a%alpha_low(i,j-1), &
        z22 => a%alpha_low(i,j))
        get_alpha = calc_four_point_linear_interpol(x, y, x1, x2, y1, y2, z11, z12, z21, z22)
        !associate (&
        !  k1 => (z12-z11) / (y2-y1), &
        !  k2 => (z22-z21) / (y2-y1), &
        !  z0_1 => z11, &
        !  z0_2 => z21)
        !  associate (&
        !    k_k  => (k2-k1)/(x2-x1), &
        !    k_0  => k1, &
        !    k_z0 => (z0_2-z0_1)/(x2-x1), &
        !    dx   => x-x1, &
        !    dy   => y-y1)
        !    z    = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
        !  end associate
        !end associate
      end associate
    end if
  end associate
end function get_alpha


function get_L0_vib()
  ! Neufeld 1993, Table 5
  double precision get_L0_vib
  double precision t1
  associate( &
    T => a_Neufeld_cooling_CO_params%T)
    t1 = exp(-log(T)/3D0)
    get_L0_vib = &
      1.83D-26 * T * exp(-68D0 * t1 - 3080D0/T)
  end associate
end function get_L0_vib


function get_L_LTE_vib()
  double precision get_L_LTE_vib
  integer i, j
  associate( &
    a      => a_Neufeld_cooling_CO_params, &
    T      => a_Neufeld_cooling_CO_params%T, &
    log10N => a_Neufeld_cooling_CO_params%log10N)
    ! Find the index.
    associate( &
      x => T, &
      y => a%T_high_vib, &
      n => n_T_high_vib, &
      idx => i)
      if      (y(1) .GE. x) then
        idx = 2
      else if (y(n) .LE. x) then
        idx = n
      else
        do idx=2, n
          if (y(idx) .GT. x) then
            exit
          end if
        end do
        if (idx .GT. n) then
          idx = n
        end if
      end if
    end associate
    associate( &
      x => log10N, &
      y => a%log10N_high_vib, &
      n => n_log10N_high_vib, &
      idx => j)
      if      (y(1) .GE. x) then
        idx = 2
      else if (y(n) .LE. x) then
        idx = n
      else
        do idx=2, n
          if (y(idx) .GT. x) then
            exit
          end if
        end do
        if (idx .GT. n) then
          idx = n
        end if
      end if
    end associate
    associate( &
      x   => log(T), &
      y   => log10N, &
      x1  => log(a%T_high_vib(i-1)), &
      x2  => log(a%T_high_vib(i)), &
      y1  => a%log10N_high_vib(j-1), &
      y2  => a%log10N_high_vib(j), &
      z11 => a%log10_X_L_LTE_high_vib(i-1,j-1), &
      z12 => a%log10_X_L_LTE_high_vib(i-1,j), &
      z21 => a%log10_X_L_LTE_high_vib(i,j-1), &
      z22 => a%log10_X_L_LTE_high_vib(i,j), &
      z => get_L_LTE_vib)
      z = calc_four_point_linear_interpol(x, y, x1, x2, y1, y2, z11, z12, z21, z22)
      !associate (&
      !  k1 => (z12-z11) / (y2-y1), &
      !  k2 => (z22-z21) / (y2-y1), &
      !  z0_1 => z11, &
      !  z0_2 => z21)
      !  associate (&
      !    k_k  => (k2-k1)/(x2-x1), &
      !    k_0  => k1, &
      !    k_z0 => (z0_2-z0_1)/(x2-x1), &
      !    dx => x-x1, &
      !    dy => y-y1)
      !    z = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
      !  end associate
      !end associate
    end associate
    get_L_LTE_vib = exp(-get_L_LTE_vib*ln10 - 3080D0/T)
  end associate
end function get_L_LTE_vib


end module load_Neufeld_cooling_CO
