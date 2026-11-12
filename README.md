# MZS elliptical-vortex model

MATLAB code for the MZS moment model of two-dimensional vortex dynamics. Each
vortex is a uniform elliptical patch, described by its centre, its aspect ratio
and its orientation, and the model evolves these under the vortices' mutual
interaction. Three geometries are included:

- **Bounded domains** given by a conformal map of the unit circle: a
  heart-shaped family that includes the unit circle itself.
- **The gap**: the plane with two walls along the real axis, leaving a gap
  between them.
- **The doubly periodic domain** $[0, 2\pi)^2$.

The repository also contains a Fortran contour-dynamics code that simulates a
single vortex patch in the gap domain. Its output is used to test the MZS model.

## Folder layout

| Folder | Contents |
|---|---|
| `MZS/MZS Model - Bounded` | Bounded-domain model (heart, unit circle and gap): main scripts, energy, equations of motion, initial conditions |
| `MZS/MZS Model - Doubly Periodic` | Doubly periodic model: main script, energy, equations of motion, initial conditions |
| `MZS/Conformal Heart` | The heart map $z = F(Z)$, its derivatives, inverse and boundary |
| `MZS/Conformal Gap` | The gap map $z = F(Z)$, its derivatives and inverse |
| `Contour Dynamics of Vortex Patch` | Fortran contour-dynamics code for a vortex patch in the gap domain, and a MATLAB script to plot its output |

Every MATLAB function in `MZS` has a help header (`help <name>` in MATLAB)
giving the formula it computes and what its inputs and outputs are.

## The model

The state of $N$ vortices is a $4N \times 1$ vector `[x; y; X; Y]`:

- $(x_i, y_i)$ is the position of the centre of vortex $i$. In the bounded
  model this is the position $Z_i = x_i + i y_i$ in the unit circle (the
  *image plane*), and the physical position is $F(Z_i)$.
- $X_i + iY_i = \left(\sqrt{\lambda_i} - 1/\sqrt{\lambda_i}\right) e^{2i\phi_i}$
  describes the shape. Here $\lambda_i$ is the aspect ratio and $\phi_i$ the
  orientation, and the semi-axes are $\sqrt{8\varepsilon_i\lambda_i}$ and
  $\sqrt{8\varepsilon_i/\lambda_i}$, so the area is $8\pi\varepsilon_i$.
  Inverting this gives $\lambda = \left(|X| + \sqrt{|X|^2+4}\right)^2/4$ and
  $\phi = \arg(X+iY)/2$.
- Each vortex also has a circulation $\kappa_i$ (`kappas`) and a size
  $\varepsilon_i$ (`eps`), both fixed.

The equations of motion are Hamilton's equations for the energy $H$:

$$
\kappa_i J_i^2 \frac{dx_i}{dt} = \frac{\partial H}{\partial y_i}, \qquad
\kappa_i J_i^2 \frac{dy_i}{dt} = -\frac{\partial H}{\partial x_i}, \qquad
\kappa_i \varepsilon_i \frac{dX_i}{dt} = \frac{\partial H}{\partial Y_i}, \qquad
\kappa_i \varepsilon_i \frac{dY_i}{dt} = -\frac{\partial H}{\partial X_i},
$$

with $J_i = |F'(Z_i)|$ in the bounded model and $J_i = 1$ in the doubly
periodic model.

| | Bounded | Doubly periodic |
|---|---|---|
| Energy $H$ (and its individual terms) | `energy_bounded_MZS.m` | `energy_doubly_MZS.m` |
| Right-hand side of the equations above | `uvUV_bounded_MZS.m` | `uvUV_doubly_MZS.m` |

The terms of $H$ are listed in the headers of the energy functions. A run
conserves $H$, and for the unit circle also the angular momentum
$L = \sum_i \kappa_i \left(|Z_i|^2 + 2\varepsilon_i(|X_i|^2 + 2)\right)$.
The main scripts print the change in each of these over the run and plot the
change in each term of $H$.

## Running the MZS model

Open MATLAB, change into the model's folder and run the script:

```matlab
cd('MZS/MZS Model - Bounded')
main_bounded_MZS
```

The parameters are set at the top of each script.

| Script | Folder | Domain | What it does |
|---|---|---|---|
| `main_bounded_MZS.m` | Bounded | Heart (`q`, `g`); `q = g = 0` is the unit circle | $N$ vortices from random initial positions with a given energy (and, for the unit circle, angular momentum) |
| `main_gap_MZS.m` | Bounded | Gap | One vortex starting above the wall, moving towards the gap |
| `Contour_Dynamics_and_MZS.m` | Bounded | Gap | As `main_gap_MZS.m`, then compares the MZS vortex with the contour-dynamics simulation (see below) |
| `main_doubly_MZS.m` | Doubly Periodic | $[0, 2\pi)^2$ | $N$ vortices, half positive and half negative, from random positions with a given energy |

Each script adds the conformal-map folder it needs to the MATLAB path
(`addpath('../Conformal Heart/')` or `addpath('../Conformal Gap/')`). Both
folders define functions with the same names (`DF1`–`DF4`, `inverse_map`,
`jacob`), and whichever folder was added last is used. Always run the scripts
from their own folder, so that this works and so the model's own functions are
found first.

### Output and restarting

The state, energy and (for the unit circle) angular momentum after every time
step `dt` are appended to a NetCDF file. The file is written in the folder
*above* the script's folder (i.e. in `MZS`), with the script's folder name as
a prefix:

| Script | Output file |
|---|---|
| `main_bounded_MZS.m` | `MZS Model - Bounded - file.nc` |
| `main_gap_MZS.m` | `MZS Model - Bounded - Gap - file.nc` |
| `Contour_Dynamics_and_MZS.m` | `MZS Model - Bounded - Contour & MZS - file.nc` |
| `main_doubly_MZS.m` | `MZS Model - Doubly Periodic - file.nc` |

**If the output file already exists, the script continues the run from its
last saved state** instead of starting a new one. Delete the file to start
afresh, for example after changing the parameters or the number of vortices.

A run stops early once any vortex becomes too elongated:

- aspect ratio 3, i.e. $|X|^2 = 4/3$ (`ElongationEvent.m`), in the bounded
  and doubly periodic runs: Love's criterion for the instability of an
  elliptical vortex;
- $|X|^2 = 25$, an aspect ratio of about 27 (`ElongationEvent2.m`), in the
  gap runs, so that the vortex can be followed for longer.

The plots at the end of each script are switched on and off with
`if 1==1` / `if 1==0`. The movie options save `Vortices.mp4` in the current
folder.

## Contour dynamics of a vortex patch

`Contour Dynamics of Vortex Patch` contains a planar contour-surgery code
(Fortran 90). It follows the boundary of a uniform vortex patch near a wall
along $y = 0$ with a gap at $|x| < 1$, the same geometry as the gap map. It
also allows optional background flows: a uniform stream along the wall, a flow
far from the wall and a flux through the gap. The code solves the
quasi-geostrophic equations, which reduce to the two-dimensional Euler
equations of the MZS model when the Rossby radius is large.

| File | Role |
|---|---|
| `patch.f90` | Program `patch`: asks for the initial patch and run parameters and writes them to `fort.11` |
| `dyn.f90` | Program `dyn`: reads `fort.11`, integrates in time (fourth-order Runge–Kutta) and appends a frame to `fort.11` at every save |
| `dyn_com.f90`, `advect_.f90`, `uv_.f90`, `cubic_.f90`, `renode_.f90`, `surgery_.f90`, `dmarray_.f90`, `dump_.f90`, `tabgen_.f90` | Modules used by `dyn`: shared data, time stepping, velocity (including the wall, the gap and the background flows), contour interpolation, re-noding, contour surgery, array allocation, output and the Green's-function table |
| `makefile` | Builds `dyn` |
| `mkmovie_v7.m` | MATLAB script that plots the frames in `fort.11` and time series of the patch's centroid, orientation, aspect ratio and area |
| `dipole.f90` | Initial conditions for two vortex patches, from the original planar code. It writes `fort.11` in an older layout (no gap or background-flow parameters, no centroid row), so it does not work with this version of `dyn` as it stands. |

### Building and running

1. **Build** both programs in the `Contour Dynamics of Vortex Patch` folder.
   The makefile calls the compiler `f95`; to use a different compiler, set
   `FC90`, e.g. `make dyn FC90=gfortran`.

   ```bash
   make dyn
   ```

   ```bash
   make patch
   ```

2. **Set up the run** with `./patch`, which asks for:
   - the reciprocal of the Rossby radius (use a small value);
   - the patch's half-widths in $x$ and $y$;
   - its vorticity (in units of $2\pi$);
   - its centre;
   - the stream velocity along the wall, the velocity far from the wall, and
     the flux through the gap;
   - the number of Fourier modes used for the flow through the gap;
   - the resolution parameter $\mu$ (at most 0.15);
   - the time step, the number of steps per save, and the number of saves.

   It writes these and the initial contour to `fort.11`, overwriting any
   existing file.
3. **Run** `./dyn`, which appends the saved frames to `fort.11`.
4. **Plot** the output, either on its own with `mkmovie_v7.m` (run from the
   folder containing `fort.11`) or against the MZS model with
   `Contour_Dynamics_and_MZS.m` (copy `fort.11` into `MZS/MZS Model - Bounded`
   first).

### Matching the MZS run in `Contour_Dynamics_and_MZS.m`

For the comparison to be like-for-like, the contour-dynamics run should
describe the same vortex and save at the same times as the MZS run:

- **Same patch.** Use a circular patch of radius $\sqrt{8\varepsilon}$, centred at the same location, with no background
  flows (stream, far-field velocity and flux all zero).
- **Same vorticity.** The MZS circulation $\kappa = 16\pi^2\varepsilon$ gives
  vorticity $\kappa/(8\pi\varepsilon) = 2\pi$. `patch` asks for the vorticity
  in units of $2\pi$, so this is an input of 1. (`dipole.f90` likewise uses 1
  for a vortex of vorticity $2\pi$.)
- **Same save times.** Frame $k$ of `fort.11` is compared with MZS step $k$,
  so the time step times the number of steps per save must equal the MZS `dt`.
- **Long enough.** Use at least as many saves as MZS steps.

`Contour_Dynamics_and_MZS.m` asks for plot settings at the keyboard. It saves
each frame as a PNG in `figfiles/` and can also save a movie, `CD & MZS.mp4`.
It then plots time series of the centroid, aspect ratio and rotation rate from
both models. The layout of `fort.11` is described in a comment in that script.

## Conformal maps

| Map | $z = F(Z)$ | Parameters |
|---|---|---|
| Heart | $\dfrac{a\,Z}{(1 - q^2 Z^2)(1 - i g Z)}$ | `domain_var = [q; g]`. $a$ = `alpha_heart(q, g, pi)` makes the domain's area $\pi$. `q = g = 0` gives the unit circle; `g = 0` gives the Neumann oval. |
| Gap | $\dfrac{2iZ}{Z^2 - 1}$ | None. The walls are $x \le -1$ and $x \ge 1$ on the real axis. The `domain_var` argument is accepted but not used. |

Random initial positions in the heart (`set_ic_bounded_MZS.m`,
`set_ic_mc_MZS.m`) use `boundary.m`, so they only work with the Conformal
Heart folder. `set_ic_boundedL_MZS.m` is for the unit circle only.

## Requirements

- **MZS model and plotting:** MATLAB R2022a or later, as the scripts use
  `Name=Value` arguments (R2021a) and `fontsize` (R2022a). No toolboxes are
  needed: `ode113`, the NetCDF functions and `VideoWriter` are all part of base
  MATLAB.
- **Contour dynamics:** a Fortran 90 compiler (e.g. `gfortran`) and `make`.
