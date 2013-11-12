module load_Neufeld_cooling_H2

use trivials

implicit none

integer, parameter, private :: &
  n_T = 22

type, private :: Neufeld_cooling_H2_params
  double precision, dimension(n_T), private :: &
    log10_T_s = (/1.6D0, 1.7D0, 1.8D0, 1.9D0, 1.95D0, 2.0D0, 2.1D0, 2.2D0, &
                  2.3D0, 2.4D0, 2.5D0, 2.6D0, 2.70D0, 2.8D0, 2.9D0, 3.0D0,  &
                  3.1D0, 3.2D0, 3.3D0, 3.4D0, 3.50D0, 3.6D0/)

  double precision, dimension(n_T), private :: &
    log10_L0    = (/24.70D0, 24.69D0, 24.71D0, 24.74D0, 24.76D0, 26.97D0, &
                    26.51D0, 26.11D0, 25.75D0, 25.42D0, 25.12D0, 24.87D0, &
                    24.64D0, 24.43D0, 24.24D0, 24.06D0, 23.87D0, 23.52D0, &
                    23.28D0, 23.01D0, 22.70D0, 22.35D0/)

  double precision, dimension(n_T), private :: &
    log10_L_LTE = (/23.04D0, 23.10D0, 23.18D0, 23.22D0, 23.18D0, 25.39D0, &
                    24.79D0, 24.25D0, 23.77D0, 23.34D0, 22.93D0, 22.53D0, &
                    22.13D0, 21.75D0, 21.38D0, 21.03D0, 20.68D0, 20.35D0, &
                    20.03D0, 19.72D0, 19.43D0, 19.16D0/)

  double precision, dimension(n_T), private :: &
    log10_n_12    = (/1.69D0, 1.60D0, 1.53D0, 1.44D0, 1.38D0, 1.38D0, 1.30D0, &
                      1.23D0, 1.16D0, 1.12D0, 1.08D0, 1.08D0, 1.14D0, 1.24D0, &
                      1.32D0, 1.75D0, 2.27D0, 2.21D0, 2.28D0, 2.34D0, 2.35D0, 2.30D0/)

  double precision, dimension(n_T), private :: &
    alpha_s    = (/0.00D0, 0.00D0, 0.00D0, 0.79D0, 0.75D0, 0.75D0, 0.67D0, &
                   0.58D0, 0.48D0, 0.42D0, 0.39D0, 0.40D0, 0.42D0, 0.43D0, &
                   0.44D0, 0.44D0, 0.55D0, 0.61D0, 0.59D0, 0.56D0, 0.57D0, 0.59D0/)

  ! The user should provide these five quantities.
  ! G is of order 1.
  ! dv/dz is expressed in km s-1 cm-1.
  ! n is in cm-3.
  ! So the N parameter has dimension cm-2 km s-1.
  double precision T

  double precision L_LTE, L0, n_12, alpha
  double precision L_LTE_vib, L0_vib

end type Neufeld_cooling_H2_params

type(Neufeld_cooling_H2_params) a_Neufeld_cooling_H2_params


contains


subroutine get_H2_rot_cool_params
  integer i
  double precision log10T
  double precision, parameter :: t1 = log(10D0)
  !
  log10T = log10(a_Neufeld_cooling_H2_params%T)
  associate( &
    a      => a_Neufeld_cooling_H2_params, &
    T      => a_Neufeld_cooling_H2_params%T)
    ! Find the index.
    associate( &
      x => log10T, &
      y => a%log10_T_s, &
      n => n_T, &
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
    ! Do the interpolation or extrapolations.
    associate( &
      x => a%log10_T_s, &
      xx => log10T)
      associate( &
        y  => a%log10_L0, &
        yy => a%L0)
        associate( &
          k  => (y(i) - y(i-1)) &
              / (x(i) - x(i-1)), &
          dx => xx - x(i-1), &
          y0 => y(i-1))
          yy = k * dx + y0
        end associate
        a%L0 = exp(-a%L0*t1 - 509D0/T)
      end associate
      associate( &
        y  => a%log10_L_LTE, &
        yy => a%L_LTE)
        associate( &
          k  => (y(i) - y(i-1)) &
              / (x(i) - x(i-1)), &
          dx => xx - x(i-1), &
          y0 => y(i-1))
          yy = k * dx + y0
        end associate
        a%L_LTE = exp(-a%L_LTE*t1 - 509D0/T)
      end associate
      associate( &
        y  => a%log10_n_12, &
        yy => a%n_12)
        associate( &
          k  => (y(i) - y(i-1)) &
              / (x(i) - x(i-1)), &
          dx => xx - x(i-1), &
          y0 => y(i-1))
          yy = k * dx + y0
        end associate
        a%n_12 = exp(a%n_12 * t1)
      end associate
      associate( &
        y  => a%alpha_s, &
        yy => a%alpha)
        associate( &
          k  => (y(i) - y(i-1)) &
              / (x(i) - x(i-1)), &
          dx => xx - x(i-1), &
          y0 => y(i-1))
          yy = k * dx + y0
        end associate
      end associate
      if (a%alpha .LT. 0D0) then
        a%alpha = 0D0
      end if
    end associate
  end associate
end subroutine get_H2_rot_cool_params


subroutine get_H2_vib_cool_params
  associate(a => a_Neufeld_cooling_H2_params, T => a_Neufeld_cooling_H2_params%T)
    a%L0 = 1.19D-24 * sqrt(T) * exp(-18100D0/(T+1190D0) - 5897D0/T)
    a%L_LTE_vib = 1.10D-18 * exp(-6744D0/T)
  end associate
end subroutine get_H2_vib_cool_params


end module load_Neufeld_cooling_H2
