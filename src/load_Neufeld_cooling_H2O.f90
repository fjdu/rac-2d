module load_Neufeld_cooling_H2O

use trivials

implicit none

integer, parameter, private :: &
  n_T_high = 6, &
  n_T_high_vib = 6, &
  n_T_low_ortho  = 6, &
  n_T_low_para  = 6, &
  n_log10N_high = 10, &
  n_log10N_high_vib = 8, &
  n_log10N_low_ortho = 10, &
  n_log10N_low_para = 10

type, private :: Neufeld_cooling_H2O_params
  double precision, dimension(n_T_high), private :: &
    T_high      = (/100D0,   200D0,   400D0,   1000D0,  2000D0,  4000D0/)
  double precision, dimension(n_T_low_ortho), private :: &
    T_low_ortho = (/10D0,    20D0,    30D0,    50D0,    80D0,    100D0/)
  double precision, dimension(n_T_low_para), private :: &
    T_low_para  = (/10D0,    20D0,    30D0,    50D0,    80D0,    100D0/)
  double precision, dimension(n_T_high), private :: &
    T_high_vib  = (/100D0,   200D0,   400D0,   1000D0,  2000D0,  4000D0/)

  double precision, dimension(n_log10N_high), private :: &
    log10N_high      = (/10.0D0, 11.0D0, 12.0D0, 13.0D0, 14.0D0, 15.0D0, 16.0D0, 17.0D0, 18.0D0, 19.0D0/)
  double precision, dimension(n_log10N_high_vib), private :: &
    log10N_high_vib  = (/13.0D0, 14.0D0, 15.0D0, 16.0D0, 17.0D0, 18.0D0, 19.0D0, 20D0/)
  double precision, dimension(n_log10N_low_ortho), private :: &
    log10N_low_ortho = (/10.0D0, 11.0D0, 12.0D0, 13.0D0, 14.0D0, 15.0D0, 16.0D0, 17.0D0, 18.0D0, 19.0D0/)
  double precision, dimension(n_log10N_low_para), private :: &
    log10N_low_para  = (/10.0D0, 11.0D0, 12.0D0, 13.0D0, 14.0D0, 15.0D0, 16.0D0, 17.0D0, 18.0D0, 19.0D0/)

  double precision, dimension(n_T_high), private :: &
    log10_L0_high      = (/24.35D0, 23.87D0, 23.42D0, 22.88D0, 22.50D0, 22.14D0/)
  double precision, dimension(n_T_low_ortho), private :: &
    log10_L0_low_ortho = (/26.81D0, 25.88D0, 25.43D0, 24.96D0, 24.58D0, 24.41D0/)
  double precision, dimension(n_T_low_para), private :: &
    log10_L0_low_para  = (/27.01D0, 25.73D0, 25.24D0, 24.75D0, 24.38D0, 24.22D0/)

  double precision, dimension(n_T_high, n_log10N_high), private :: &
    log10_L_LTE_high = reshape((/ &
      14.59D0, 13.85D0, 13.16D0, 12.32D0, 11.86D0, 11.64D0, &
      14.59D0, 13.86D0, 13.16D0, 12.32D0, 11.86D0, 11.64D0, &
      14.60D0, 13.86D0, 13.16D0, 12.32D0, 11.86D0, 11.64D0, &
      14.68D0, 13.88D0, 13.17D0, 12.32D0, 11.86D0, 11.64D0, &
      14.98D0, 14.05D0, 13.25D0, 12.34D0, 11.87D0, 11.65D0, &
      15.53D0, 14.46D0, 13.53D0, 12.49D0, 11.97D0, 11.72D0, &
      16.22D0, 15.05D0, 14.02D0, 12.87D0, 12.35D0, 12.06D0, &
      17.00D0, 15.74D0, 14.63D0, 13.46D0, 12.97D0, 12.66D0, &
      17.83D0, 16.50D0, 15.32D0, 14.16D0, 13.69D0, 13.36D0, &
      18.70D0, 17.31D0, 16.07D0, 14.94D0, 14.46D0, 14.13D0  &
      /), (/n_T_high, n_log10N_high/))

  double precision, dimension(n_T_high_vib, n_log10N_high_vib), private :: &
    log10_X_L_LTE_high_vib = reshape((/ &
      10.98D0, 11.05D0, 11.07D0, 10.88D0, 10.67D0, 10.67D0, &
      10.99D0, 11.05D0, 11.07D0, 10.88D0, 10.67D0, 10.67D0, &
      11.06D0, 11.09D0, 11.09D0, 10.88D0, 10.67D0, 10.67D0, &
      11.36D0, 11.31D0, 11.23D0, 10.93D0, 10.69D0, 10.68D0, &
      11.91D0, 11.82D0, 11.65D0, 11.15D0, 10.79D0, 10.72D0, &
      12.59D0, 12.50D0, 12.28D0, 11.61D0, 11.11D0, 10.90D0, &
      13.29D0, 13.20D0, 12.99D0, 12.25D0, 11.65D0, 11.30D0, &
      13.90D0, 13.90D0, 13.76D0, 13.01D0, 12.34D0, 11.88D0  &
      /), (/n_T_high_vib, n_log10N_high_vib/))

  double precision, dimension(n_T_low_ortho, n_log10N_low_ortho), private :: &
    log10_L_LTE_low_ortho = reshape((/ &
      17.94D0, 16.71D0, 16.08D0, 15.41D0, 14.85D0, 14.60D0, &
      17.96D0, 16.72D0, 16.09D0, 15.42D0, 14.86D0, 14.60D0, &
      18.14D0, 16.86D0, 16.19D0, 15.47D0, 14.88D0, 14.62D0, &
      18.77D0, 17.36D0, 16.58D0, 15.72D0, 15.02D0, 14.73D0, &
      19.70D0, 18.11D0, 17.25D0, 16.27D0, 15.47D0, 15.11D0, &
      20.67D0, 18.93D0, 18.05D0, 17.01D0, 16.12D0, 15.72D0, &
      21.58D0, 19.81D0, 18.91D0, 17.82D0, 16.88D0, 16.45D0, &
      22.53D0, 20.71D0, 19.80D0, 18.69D0, 17.70D0, 17.25D0, &
      23.50D0, 21.64D0, 20.72D0, 19.58D0, 18.56D0, 18.10D0, &
      24.41D0, 22.58D0, 21.65D0, 20.49D0, 19.46D0, 18.99D0  &
      /), (/n_T_low_ortho, n_log10N_low_ortho/))

  double precision, dimension(n_T_low_para, n_log10N_low_para), private :: &
    log10_L_LTE_low_para = reshape((/ &
      17.72D0, 16.60D0, 16.12D0, 15.43D0, 14.86D0, 14.60D0, &
      17.76D0, 16.63D0, 16.13D0, 15.43D0, 14.86D0, 14.60D0, &
      18.07D0, 16.85D0, 16.24D0, 15.48D0, 14.88D0, 14.62D0, &
      18.83D0, 17.41D0, 16.61D0, 15.72D0, 15.03D0, 14.73D0, &
      19.68D0, 18.11D0, 17.26D0, 16.28D0, 15.47D0, 15.11D0, &
      20.50D0, 18.94D0, 18.06D0, 17.01D0, 16.12D0, 15.72D0, &
      21.37D0, 19.83D0, 18.93D0, 17.82D0, 16.87D0, 16.45D0, &
      22.28D0, 20.75D0, 19.81D0, 18.69D0, 17.70D0, 17.25D0, &
      23.22D0, 21.69D0, 20.73D0, 19.58D0, 18.56D0, 18.10D0, &
      24.19D0, 22.60D0, 21.65D0, 20.49D0, 19.45D0, 18.99D0  &
      /), (/n_T_low_para, n_log10N_low_para/))

  double precision, dimension(n_T_high, n_log10N_high), private :: &
    log10_n_12_high = reshape((/ &
       9.00D0,  9.04D0,  9.19D0,  9.50D0,  9.67D0,  9.60D0, &
       8.99D0,  9.04D0,  9.19D0,  9.50D0,  9.67D0,  9.60D0, &
       8.96D0,  9.03D0,  9.19D0,  9.50D0,  9.66D0,  9.59D0, &
       8.74D0,  8.89D0,  9.11D0,  9.47D0,  9.65D0,  9.59D0, &
       8.11D0,  8.37D0,  8.73D0,  9.31D0,  9.56D0,  9.53D0, &
       7.20D0,  7.51D0,  7.95D0,  8.74D0,  9.15D0,  9.20D0, &
       6.22D0,  6.53D0,  6.99D0,  7.87D0,  8.38D0,  8.50D0, &
       5.22D0,  5.57D0,  6.03D0,  6.94D0,  7.48D0,  7.64D0, &
       4.24D0,  4.59D0,  5.09D0,  6.02D0,  6.59D0,  6.78D0, &
       3.21D0,  3.58D0,  4.10D0,  5.08D0,  5.69D0,  5.89D0  &
      /), (/n_T_high, n_log10N_high/))

  double precision, dimension(n_T_low_ortho, n_log10N_low_ortho), private :: &
    log10_n_12_low_ortho = reshape((/ &
       8.81D0,  8.90D0,  9.03D0,  9.14D0,  9.19D0,  9.20D0, &
       8.79D0,  8.88D0,  9.01D0,  9.13D0,  9.20D0,  9.20D0, &
       8.60D0,  8.73D0,  8.88D0,  9.03D0,  9.13D0,  9.15D0, &
       7.96D0,  8.14D0,  8.31D0,  8.55D0,  8.78D0,  8.85D0, &
       7.01D0,  7.21D0,  7.40D0,  7.69D0,  8.00D0,  8.11D0, &
       6.02D0,  6.20D0,  6.41D0,  6.71D0,  7.04D0,  7.17D0, &
       5.03D0,  5.21D0,  5.42D0,  5.70D0,  6.06D0,  6.19D0, &
       4.02D0,  4.21D0,  4.41D0,  4.71D0,  5.05D0,  5.19D0, &
       3.02D0,  3.22D0,  3.41D0,  3.71D0,  4.06D0,  4.20D0, &
       2.03D0,  2.21D0,  2.42D0,  2.70D0,  3.06D0,  3.20D0  &
      /), (/n_T_low_ortho, n_log10N_low_ortho/))

  double precision, dimension(n_T_low_para, n_log10N_low_para), private :: &
    log10_n_12_low_para = reshape((/ &
       9.30D0,  9.11D0,  8.94D0,  8.70D0,  8.55D0,  8.50D0, &
       9.25D0,  9.06D0,  8.91D0,  8.69D0,  8.54D0,  8.49D0, &
       8.94D0,  8.82D0,  8.71D0,  8.56D0,  8.46D0,  8.43D0, &
       8.17D0,  8.16D0,  8.16D0,  8.12D0,  8.10D0,  8.11D0, &
       7.21D0,  7.24D0,  7.28D0,  7.30D0,  7.32D0,  7.35D0, &
       6.20D0,  6.24D0,  6.31D0,  6.33D0,  6.35D0,  6.38D0, &
       5.21D0,  5.25D0,  5.30D0,  5.32D0,  5.35D0,  5.39D0, &
       4.22D0,  4.26D0,  4.31D0,  4.33D0,  4.36D0,  4.40D0, &
       3.23D0,  3.24D0,  3.31D0,  3.33D0,  3.37D0,  3.41D0, &
       2.21D0,  2.25D0,  2.30D0,  2.32D0,  2.37D0,  2.41D0  &
      /), (/n_T_low_para, n_log10N_low_para/))

  double precision, dimension(n_T_high, n_log10N_high), private :: &
    alpha_high = reshape((/ &
       0.43D0,  0.42D0,  0.39D0,  0.36D0,  0.34D0,  0.34D0, &
       0.43D0,  0.42D0,  0.39D0,  0.36D0,  0.34D0,  0.34D0, &
       0.42D0,  0.41D0,  0.39D0,  0.36D0,  0.34D0,  0.34D0, &
       0.41D0,  0.39D0,  0.37D0,  0.35D0,  0.33D0,  0.33D0, &
       0.42D0,  0.38D0,  0.34D0,  0.33D0,  0.32D0,  0.32D0, &
       0.45D0,  0.38D0,  0.34D0,  0.32D0,  0.30D0,  0.30D0, &
       0.47D0,  0.40D0,  0.35D0,  0.32D0,  0.29D0,  0.30D0, &
       0.50D0,  0.42D0,  0.36D0,  0.32D0,  0.28D0,  0.29D0, &
       0.52D0,  0.44D0,  0.37D0,  0.31D0,  0.27D0,  0.28D0, &
       0.53D0,  0.45D0,  0.39D0,  0.31D0,  0.27D0,  0.27D0  &
      /), (/n_T_high, n_log10N_high/))

  double precision, dimension(n_T_low_ortho, n_log10N_low_ortho), private :: &
    alpha_low_ortho = reshape((/ &
       0.71D0,  0.49D0,  0.48D0,  0.46D0,  0.45D0,  0.45D0, &
       0.64D0,  0.49D0,  0.48D0,  0.45D0,  0.45D0,  0.45D0, &
       0.57D0,  0.50D0,  0.47D0,  0.45D0,  0.44D0,  0.44D0, &
       0.56D0,  0.53D0,  0.48D0,  0.44D0,  0.42D0,  0.41D0, &
       0.59D0,  0.60D0,  0.53D0,  0.47D0,  0.45D0,  0.43D0, &
       0.72D0,  0.64D0,  0.58D0,  0.52D0,  0.49D0,  0.47D0, &
       0.85D0,  0.68D0,  0.61D0,  0.56D0,  0.52D0,  0.50D0, &
       0.86D0,  0.70D0,  0.63D0,  0.58D0,  0.55D0,  0.53D0, &
       0.93D0,  0.72D0,  0.66D0,  0.61D0,  0.57D0,  0.55D0, &
       0.87D0,  0.73D0,  0.67D0,  0.62D0,  0.59D0,  0.56D0  &
      /), (/n_T_low_ortho, n_log10N_low_ortho/))

  double precision, dimension(n_T_low_para, n_log10N_low_para), private :: &
    alpha_low_para = reshape((/ &
       0.49D0,  0.72D0,  0.69D0,  0.53D0,  0.46D0,  0.44D0, &
       0.52D0,  0.65D0,  0.66D0,  0.53D0,  0.46D0,  0.44D0, &
       0.49D0,  0.65D0,  0.63D0,  0.51D0,  0.45D0,  0.43D0, &
       0.57D0,  0.73D0,  0.65D0,  0.49D0,  0.43D0,  0.41D0, &
       0.75D0,  0.68D0,  0.64D0,  0.52D0,  0.45D0,  0.42D0, &
       0.76D0,  0.70D0,  0.68D0,  0.55D0,  0.47D0,  0.43D0, &
       0.79D0,  0.73D0,  0.69D0,  0.58D0,  0.48D0,  0.44D0, &
       0.80D0,  0.75D0,  0.71D0,  0.58D0,  0.50D0,  0.46D0, &
       0.82D0,  0.77D0,  0.73D0,  0.60D0,  0.52D0,  0.48D0, &
       0.91D0,  0.80D0,  0.74D0,  0.62D0,  0.53D0,  0.50D0  &
      /), (/n_T_low_para, n_log10N_low_para/))

  ! The user should provide these five quantities.
  ! G is of order 1.
  ! dv/dz is expressed in km s-1 cm-1.
  ! n is in cm-3.
  ! So the N parameter has dimension cm-2 km s-1.
  double precision T, log10N

  double precision L, L_LTE, L0, n_12, alpha
  double precision L_LTE_vib, L0_vib

  double precision :: ortho_ratio=0.75D0, para_ratio=0.25D0

end type Neufeld_cooling_H2O_params

type(Neufeld_cooling_H2O_params) a_Neufeld_cooling_H2O_params

double precision, parameter, private :: ln10 = log(10D0)


contains


function get_L0()
  double precision get_L0
  integer i
  associate( &
    a => a_Neufeld_cooling_H2O_params, &
    T => a_Neufeld_cooling_H2O_params%T)
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
            !/ (a%T_high(i) - a%T_high(i-1)), &
        dx => log(T) - log(a%T_high(i-1)), &
        !dx => T - a%T_high(i-1), &
        y0 => a%log10_L0_high(i-1))
        get_L0 = k * dx + y0
      end associate
    else
      associate( &
        x => T, &
        y => a%T_low_ortho, &
        n => n_T_low_ortho, &
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
        k1  => (a%log10_L0_low_ortho(i) - a%log10_L0_low_ortho(i-1)) &
             / (a%T_low_ortho(i) - a%T_low_ortho(i-1)), &
        k2  => (a%log10_L0_low_para(i)  - a%log10_L0_low_para(i-1)) &
             / (a%T_low_para(i)  - a%T_low_para(i-1)), &
        dx1 => T - a%T_low_ortho(i-1), &
        dx2 => T - a%T_low_para(i-1), &
        y0_1 => a%log10_L0_low_ortho(i-1), &
        y0_2 => a%log10_L0_low_para(i-1))
        get_L0 = a%ortho_ratio * (dx1 * k1 + y0_1) + &
                 a%para_ratio  * (dx2 * k2 + y0_2)
      end associate
    end if
  end associate
  get_L0 = exp(-get_L0*ln10)
end function get_L0


function get_L_LTE()
  double precision get_L_LTE
  double precision tmp1, tmp2
  integer i, j
  associate( &
    a => a_Neufeld_cooling_H2O_params, &
    T => a_Neufeld_cooling_H2O_params%T, &
    log10N => a_Neufeld_cooling_H2O_params%log10N)
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
        !associate( &
        !  k1 => (z12-z11) / (y2-y1), &
        !  k2 => (z22-z21) / (y2-y1), &
        !  z0_1 => z11, &
        !  z0_2 => z21)
        !  associate( &
        !    k_k  => (k2-k1)/(x2-x1), &
        !    k_0  => k1, &
        !    k_z0 => (z0_2-z0_1)/(x2-x1), &
        !    dx => x-x1, &
        !    dy => y-y1)
        !    z = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
        !  end associate
        !end associate
      end associate
    else
      associate( &
        x => T, &
        y => a%T_low_ortho, &
        n => n_T_low_ortho, &
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
        y => a%log10N_low_ortho, &
        n => n_log10N_low_ortho, &
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
        x1  => log(a%T_low_ortho(i-1)), &
        x2  => log(a%T_low_ortho(i)), &
        y1  => a%log10N_low_ortho(j-1), &
        y2  => a%log10N_low_ortho(j), &
        z11 => a%log10_L_LTE_low_ortho(i-1,j-1), &
        z12 => a%log10_L_LTE_low_ortho(i-1,j), &
        z21 => a%log10_L_LTE_low_ortho(i,j-1), &
        z22 => a%log10_L_LTE_low_ortho(i,j), &
        z => tmp1)
        z = calc_four_point_linear_interpol(x, y, x1, x2, y1, y2, z11, z12, z21, z22)
        !associate( &
        !  k1 => (z12-z11) / (y2-y1), &
        !  k2 => (z22-z21) / (y2-y1), &
        !  z0_1 => z11, &
        !  z0_2 => z21)
        !  associate( &
        !    k_k  => (k2-k1)/(x2-x1), &
        !    k_0  => k1, &
        !    k_z0 => (z0_2-z0_1)/(x2-x1), &
        !    dx => x-x1, &
        !    dy => y-y1)
        !    z = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
        !  end associate
        !end associate
      end associate
      associate( &
        x   => log(T), &
        y   => log10N, &
        x1  => log(a%T_low_para(i-1)), &
        x2  => log(a%T_low_para(i)), &
        y1  => a%log10N_low_para(j-1), &
        y2  => a%log10N_low_para(j), &
        z11 => a%log10_L_LTE_low_para(i-1,j-1), &
        z12 => a%log10_L_LTE_low_para(i-1,j), &
        z21 => a%log10_L_LTE_low_para(i,j-1), &
        z22 => a%log10_L_LTE_low_para(i,j), &
        z => tmp2)
        z = calc_four_point_linear_interpol(x, y, x1, x2, y1, y2, z11, z12, z21, z22)
        !associate( &
        !  k1 => (z12-z11) / (y2-y1), &
        !  k2 => (z22-z21) / (y2-y1), &
        !  z0_1 => z11, &
        !  z0_2 => z21)
        !  associate( &
        !    k_k  => (k2-k1)/(x2-x1), &
        !    k_0  => k1, &
        !    k_z0 => (z0_2-z0_1)/(x2-x1), &
        !    dx => x-x1, &
        !    dy => y-y1)
        !    z = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
        !  end associate
        !end associate
      end associate
      get_L_LTE = a%ortho_ratio * tmp1 + a%para_ratio * tmp2
    end if
  end associate
  get_L_LTE = exp(-get_L_LTE*ln10)
end function get_L_LTE


function get_n_12()
  double precision get_n_12
  double precision tmp1, tmp2
  integer i, j
  associate( &
    a => a_Neufeld_cooling_H2O_params, &
    T => a_Neufeld_cooling_H2O_params%T, &
    log10N => a_Neufeld_cooling_H2O_params%log10N)
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
        !associate( &
        !  k1 => (z12-z11) / (y2-y1), &
        !  k2 => (z22-z21) / (y2-y1), &
        !  z0_1 => z11, &
        !  z0_2 => z21)
        !  associate( &
        !    k_k  => (k2-k1)/(x2-x1), &
        !    k_0  => k1, &
        !    k_z0 => (z0_2-z0_1)/(x2-x1), &
        !    dx => x-x1, &
        !    dy => y-y1)
        !    z = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
        !  end associate
        !end associate
      end associate
    else
      associate( &
        x => T, &
        y => a%T_low_ortho, &
        n => n_T_low_ortho, &
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
        y => a%log10N_low_ortho, &
        n => n_log10N_low_ortho, &
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
        x1  => log(a%T_low_ortho(i-1)), &
        x2  => log(a%T_low_ortho(i)), &
        y1  => a%log10N_low_ortho(j-1), &
        y2  => a%log10N_low_ortho(j), &
        z11 => a%log10_n_12_low_ortho(i-1,j-1), &
        z12 => a%log10_n_12_low_ortho(i-1,j), &
        z21 => a%log10_n_12_low_ortho(i,j-1), &
        z22 => a%log10_n_12_low_ortho(i,j), &
        z => tmp1)
        z = calc_four_point_linear_interpol(x, y, x1, x2, y1, y2, z11, z12, z21, z22)
        !associate( &
        !  k1 => (z12-z11) / (y2-y1), &
        !  k2 => (z22-z21) / (y2-y1), &
        !  z0_1 => z11, &
        !  z0_2 => z21)
        !  associate( &
        !    k_k  => (k2-k1)/(x2-x1), &
        !    k_0  => k1, &
        !    k_z0 => (z0_2-z0_1)/(x2-x1), &
        !    dx => x-x1, &
        !    dy => y-y1)
        !    z = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
        !  end associate
        !end associate
      end associate
      associate( &
        x   => log(T), &
        y   => log10N, &
        x1  => log(a%T_low_para(i-1)), &
        x2  => log(a%T_low_para(i)), &
        y1  => a%log10N_low_para(j-1), &
        y2  => a%log10N_low_para(j), &
        z11 => a%log10_n_12_low_para(i-1,j-1), &
        z12 => a%log10_n_12_low_para(i-1,j), &
        z21 => a%log10_n_12_low_para(i,j-1), &
        z22 => a%log10_n_12_low_para(i,j), &
        z => tmp2)
        z = calc_four_point_linear_interpol(x, y, x1, x2, y1, y2, z11, z12, z21, z22)
        !associate( &
        !  k1 => (z12-z11) / (y2-y1), &
        !  k2 => (z22-z21) / (y2-y1), &
        !  z0_1 => z11, &
        !  z0_2 => z21)
        !  associate( &
        !    k_k  => (k2-k1)/(x2-x1), &
        !    k_0  => k1, &
        !    k_z0 => (z0_2-z0_1)/(x2-x1), &
        !    dx => x-x1, &
        !    dy => y-y1)
        !    z = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
        !  end associate
        !end associate
      end associate
      get_n_12 = a%ortho_ratio * tmp1 + a%para_ratio * tmp2
    end if
  end associate
  get_n_12 = exp(-get_n_12*ln10)
end function get_n_12


function get_alpha()
  double precision get_alpha
  double precision tmp1, tmp2
  integer i, j
  associate( &
    a => a_Neufeld_cooling_H2O_params, &
    T => a_Neufeld_cooling_H2O_params%T, &
    log10N => a_Neufeld_cooling_H2O_params%log10N)
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
        z22 => a%alpha_high(i,j), &
        z => get_alpha)
        z = calc_four_point_linear_interpol(x, y, x1, x2, y1, y2, z11, z12, z21, z22)
        !associate( &
        !  k1 => (z12-z11) / (y2-y1), &
        !  k2 => (z22-z21) / (y2-y1), &
        !  z0_1 => z11, &
        !  z0_2 => z21)
        !  associate( &
        !    k_k  => (k2-k1)/(x2-x1), &
        !    k_0  => k1, &
        !    k_z0 => (z0_2-z0_1)/(x2-x1), &
        !    dx => x-x1, &
        !    dy => y-y1)
        !    z = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
        !  end associate
        !end associate
      end associate
    else
      associate( &
        x => T, &
        y => a%T_low_ortho, &
        n => n_T_low_ortho, &
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
        y => a%log10N_low_ortho, &
        n => n_log10N_low_ortho, &
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
        x1  => log(a%T_low_ortho(i-1)), &
        x2  => log(a%T_low_ortho(i)), &
        y1  => a%log10N_low_ortho(j-1), &
        y2  => a%log10N_low_ortho(j), &
        z11 => a%alpha_low_ortho(i-1,j-1), &
        z12 => a%alpha_low_ortho(i-1,j), &
        z21 => a%alpha_low_ortho(i,j-1), &
        z22 => a%alpha_low_ortho(i,j), &
        z => tmp1)
        z = calc_four_point_linear_interpol(x, y, x1, x2, y1, y2, z11, z12, z21, z22)
        !associate( &
        !  k1 => (z12-z11) / (y2-y1), &
        !  k2 => (z22-z21) / (y2-y1), &
        !  z0_1 => z11, &
        !  z0_2 => z21)
        !  associate( &
        !    k_k  => (k2-k1)/(x2-x1), &
        !    k_0  => k1, &
        !    k_z0 => (z0_2-z0_1)/(x2-x1), &
        !    dx => x-x1, &
        !    dy => y-y1)
        !    z = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
        !  end associate
        !end associate
      end associate
      associate( &
        x   => log(T), &
        y   => log10N, &
        x1  => log(a%T_low_para(i-1)), &
        x2  => log(a%T_low_para(i)), &
        y1  => a%log10N_low_para(j-1), &
        y2  => a%log10N_low_para(j), &
        z11 => a%alpha_low_para(i-1,j-1), &
        z12 => a%alpha_low_para(i-1,j), &
        z21 => a%alpha_low_para(i,j-1), &
        z22 => a%alpha_low_para(i,j), &
        z => tmp2)
        z = calc_four_point_linear_interpol(x, y, x1, x2, y1, y2, z11, z12, z21, z22)
        !associate( &
        !  k1 => (z12-z11) / (y2-y1), &
        !  k2 => (z22-z21) / (y2-y1), &
        !  z0_1 => z11, &
        !  z0_2 => z21)
        !  associate( &
        !    k_k  => (k2-k1)/(x2-x1), &
        !    k_0  => k1, &
        !    k_z0 => (z0_2-z0_1)/(x2-x1), &
        !    dx => x-x1, &
        !    dy => y-y1)
        !    z = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
        !  end associate
        !end associate
      end associate
      get_alpha = a%ortho_ratio * tmp1 + a%para_ratio * tmp2
    end if
  end associate
end function get_alpha


function get_L0_vib()
  ! Neufeld 1993, Table 5
  double precision get_L0_vib
  double precision t1
  associate(T => a_Neufeld_cooling_H2O_params%T)
    t1 = exp(-log(T)/3D0)
    get_L0_vib = &
      1.03D-26 * T * exp(-47.5D0 * t1 - 2325D0/T)
  end associate
end function get_L0_vib


function get_L_LTE_vib()
  double precision get_L_LTE_vib
  integer i, j
  associate( &
    a      => a_Neufeld_cooling_H2O_params, &
    T      => a_Neufeld_cooling_H2O_params%T, &
    log10N => a_Neufeld_cooling_H2O_params%log10N)
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
      !associate( &
      !  k1 => (z12-z11) / (y2-y1), &
      !  k2 => (z22-z21) / (y2-y1), &
      !  z0_1 => z11, &
      !  z0_2 => z21)
      !  associate( &
      !    k_k  => (k2-k1)/(x2-x1), &
      !    k_0  => k1, &
      !    k_z0 => (z0_2-z0_1)/(x2-x1), &
      !    dx => x-x1, &
      !    dy => y-y1)
      !    z = (k_k * dx + k_0) * dy + k_z0 * dx + z0_1
      !  end associate
      !end associate
    end associate
    get_L_LTE_vib = exp(-get_L_LTE_vib*ln10 - 2325D0/T)
  end associate
end function get_L_LTE_vib


end module load_Neufeld_cooling_H2O
