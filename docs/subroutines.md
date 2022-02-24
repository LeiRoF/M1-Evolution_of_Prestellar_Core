> Please note that `dp` is defined in the Constat.f90 file. It determine the size of all real variables of the program. By default, it is set to 8.

# density 1D

Compute the density along the radius of a sphere containing a mass distributed in a given way.

$$\rho_{1D} = \frac{Shape}{\int_0^{2 \pi}\int_0^{\pi}\int_{r_{min}}^{r_{max}} Shape * r^2 \sin(\theta) dr d\theta d\phi} * Mass$$

```Fortran
call mean_density_1D(shape, mass, N, r, dr, rho)
```

Inputs:
```Fortran
real (kind=dp), dimension(N) :: shape ! shape of the density distribution
real (kind=dp)               :: mass  ! mass of the system
integer                      :: N     ! array dimension
real (kind=dp), dimension(N) :: r     ! space
real (kind=dp)               :: dr    ! space step
```

Outputs
```Fortran
real (kind=dp), dimension(N) :: rho   ! distribution of the density
```

# derive 🔰2.1

Get the derivated function. A function (and it's derivated form) is represented by an array containing discret values of this function. It call `darive_at` for each element of the array.

```Fortran
call derive(f, dx, N, res)
```

Inputs
```Fortran
real (kind=dp), dimension(N) :: f     ! array
real (kind=dp)               :: dx    ! step of the array
integer                      :: N     ! size of the array
```

Ouputs
```Fortran
real (kind=dp), dimension(N) :: res   ! derivated array
```

# derive at

Get the derivative at a given point using an array containing some discret values of a function.

$$ \frac{\partial F}{\partial r} \approx \frac{F(r+dr) - F(r)}{dr} $$

Numercially, we use the centred differential method:

$$ F'_i  = \frac{F_{i+1} - F_{i-1}}{2dr} $$

Except for the first and last element of the array, where we use respectively right and left differential method. Note that in this subroutine, `x` represent the $i$ index, and `dx` represent $dr$.

```Fortran
call derive_at(f, x, dx, N, res)
```

Inputs
```Fortran
real (kind=dp), dimension(N) :: f     ! array
integer                      :: x     ! position in the array where you want to get the derivative
real (kind=dp)               :: dx    ! step of the array
integer                      :: N     ! size of the array
```

Outputs
```Fortran
real (kind=dp)               :: res   ! value of the derivative
```

# get cs

Get isothermal speed of sound

$$ c_s = \sqrt{\frac{k_B T}{\mu m_H}} $$

```Fortran
call get_cs(T, res)
```

Inputs
```Fortran
real (kind=dp) :: T   ! temprature [K]
```

Outputs
```Fortran
real (kind=dp) :: res ! isothermal speed of sound
```

# grav force 🔰2.1

Get gravitational force

$$ \tag{3} \nabla^2 \Phi = 4 \pi G \rho $$

$$ \rightarrow \frac{\partial \Phi}{\partial r} =  \int 4 \pi G \rho dr $$

```Fortran
call grav_force(rho, N, r, dr, dphi)
```

Inputs
```Fortran
real (kind=dp), dimension(N) :: rho  ! density array
real (kind=dp), dimension(N) :: r    ! space array
integer                      :: N    ! array dimension
real (kind=dp)               :: dr   ! space step
```

Outputs
```Fortran
real (kind=dp)               :: dphi ! gravitational force
```

# integrate on sphere

Integrate a 1D distribution on a sphere

$$ I = \int_0^{2\pi} \int_0^{\pi} \int_{r_{min}}^{r_{max}} D \: dr \: d\theta  \: d\phi$$

```Fortran
call integrate_on_sphere(D, N, r, dr, res)
```

Inputs
```Fortran
real (kind=dp), dimension(N) :: D    ! 1D distribution array to integrate
integer                      :: N    ! array dimension
real (kind=dp), dimension(N) :: r    ! space array
real (kind=dp)               :: dr   ! space step
```

Outputs
```Fortran
real (kind=dp)               :: res ! result
```

# linspace

Get a space linear space array

```Fortran
call linspace(rmin, rmax, N, res, dr)
```

Inputs
```Fortran
integer                      :: N    ! number of discrete space points
real (kind=dp)               :: rmin ! starting position
real (kind=dp)               :: rmax ! stop position
```

Outputs
```Fortran
real (kind=dp), dimension(N) :: res  ! space array
real (kind=dp)               :: dr   ! space step
```

# next rho

Get the density profile in an array at time i+1

$$ \rho_{i+1} \approx \rho_i - \frac{1}{r^2}\bigg(\frac{\partial(r^2\rho v)}{\partial r}\bigg)_i dt$$


```Fortran
call next_rho(rho, v, r, N, dr, dt, res)
```

Inputs
```Fortran
real (kind=dp), dimension(N) :: rho ! array containing densities at time i
real (kind=dp), dimension(N) :: v   ! array containing speeds    at time i
real (kind=dp), dimension(N) :: r   ! array representing space
integer                      :: N   ! dimension of the arrays
real (kind=dp)               :: dt  ! time step
real (kind=dp)               :: dr  ! space step
```

Outputs
```Fortran
real (kind=dp), dimension(N) :: res ! returning array rho        at time i+1
```

# next state 🔰2.1

Simulation evolution of the state at time i (get the state at time i+1). Note that this subroutine just call `next_rho(rho, v, r, N, dr, dt, res)` and `next_v(v, rho, r, N, dr, dt, cs, res)`.

```Fortran
call next_state(rho, v, r, N, dr, dt, cs, new_rho, new_v)
```

Inputs
```Fortran
real (kind=dp), dimension(N) :: rho     ! array containing densities at time i
real (kind=dp), dimension(N) :: v       ! array containing speeds    at time i
real (kind=dp), dimension(N) :: r       ! array representing space
integer                      :: N       ! dimension of the arrays
real (kind=dp)               :: cs      ! isothermal speed of sound
real (kind=dp)               :: dt      ! time step
real (kind=dp)               :: dr      ! space step
```

Outputs
```Fortran
real (kind=dp), dimension(N) :: new_rho ! returning density          at time i+1
real (kind=dp), dimension(N) :: new_v   ! returning speeds           at time i+1
```

# next v

Get the speed profile in an array at time i+1

$$ v_{i+1} \approx v_i + \bigg[-v_i\bigg(\frac{\partial v}{\partial r}\bigg)_i - \frac{c_s^2}{\rho_i}\bigg(\frac{\partial \rho}{\partial r}\bigg)_i - \int 4 \pi G \rho_i dr \bigg] dt $$

```Fortran
call next_v(v, rho, r, N, dr, dt, cs, res)
```

Inputs
```Fortran
real (kind=dp), dimension(N) :: v   ! array containing speeds    at time i
real (kind=dp), dimension(N) :: rho ! array containing densities at time i
real (kind=dp), dimension(N) :: r   ! array representing space
integer                      :: N   ! dimension of the arrays
real (kind=dp)               :: cs  ! isothermal speed of sound
real (kind=dp)               :: dt  ! time step
real (kind=dp)               :: dr  ! space step
```

Outputs
```Fortran
real (kind=dp), dimension(N) :: res ! returning speeds           at time i+1
```

# time free fall

Get the time of free fall for a gas that have an initial density of rho_0

```Fortran
call time_free_fall(rho_0, t)
```

Inputs
```Fortran
real (kind=dp) :: rho_0 ! initial density [g.cm-3]
```

Ouputs
```Fortran
real (kind=dp) :: t     ! time of free fall [s]
```

# total mass 🔰2.2

Get the total mass of the system. Note that this subroutine is equivalent to `integrate_on_sphere(D, N, r, dr, res)`

```Fortran
call total_mass(rho, N, r, dr, res)
```

! Inputs
```Fortran
real (kind=dp), dimension(N) :: rho  ! density array at a given time
integer                      :: N    ! array dimension
real (kind=dp), dimension(N) :: r    ! space array
real (kind=dp)               :: dr   ! space step
```

! Outputs
```Fortran
real (kind=dp)               :: res  ! total mass considered
```