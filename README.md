# MZS elliptical-vortex model

MATLAB code for the MZS moment model of two-dimensional vortex dynamics. Each
vortex is a uniform elliptical patch, described by its centre, its aspect ratio
and its orientation, and the model evolves these under the vortices' mutual
interaction. Three geometries are included:

- **Bounded domains**: given by a conformal map of the unit circle: a
  heart-shaped family that includes the circle itself, as well as the Neumann oval.
- **The gap**: the plane with two walls along the real axis, leaving a gap
  between them.
- **The doubly periodic domain**: $[0, 2\pi)^2$.

Fortran code is also available for Contour Dynamics comparisons.

## MZS folder layout

| Folder | Contents |
|---|---|
| `MZS Model - Bounded` | Bounded-domain model (heart, circle and gap): main scripts, energy, equations of motion, initial conditions |
| `MZS Model - Doubly Periodic` | Doubly periodic model: main script, energy, equations of motion, initial conditions |
| `Conformal Heart` | The heart map $z = F(Z)$, its derivatives, inverse and boundary |
| `Conformal Gap` | The gap map $z = F(Z)$, its derivatives and inverse |

Every function has a help header (`help <name>` in MATLAB) giving the formula
it computes and what its inputs and outputs are.

## The model

The state of $N$ vortices is a $4N \times 1$ vector `[x; y; X; Y]`:

- $(x_i, y_i)$ is the position of the centre of vortex $i$. In the bounded
  model this is the position $Z_i = x_i + i y_i$ in the unit disc (the *image
  plane*), and the physical position is $F(Z_i)$.
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
conserves $H$, and for the unit disc also the angular momentum
$L = \sum_i \kappa_i \left(|Z_i|^2 + 2\varepsilon_i(|X_i|^2 + 2)\right)$.
The main scripts print the change in each of these over the run and plot the
change in each term of $H$.

## Running

Open MATLAB, change into the model's folder and run the script:

```matlab
cd('MZS Model - Bounded')
main_bounded_MZS
```

The parameters are set at the top of each script.

| Script | Folder | Domain | What it does |
|---|---|---|---|
| `main_bounded_MZS.m` | Bounded | Heart (`q`, `g`); `q = g = 0` is the unit disc | $N$ vortices from random initial positions with a given energy (and, for the disc, angular momentum) |
| `main_gap_MZS.m` | Bounded | Gap | One vortex starting above the wall, moving towards the gap |
| `Contour_Dynamics_and_MZS.m` | Bounded | Gap | As `main_gap_MZS.m`, then compares the MZS vortex with a contour-dynamics simulation (see below) |
| `main_doubly_MZS.m` | Doubly Periodic | $[0, 2\pi)^2$ | $N$ vortices, half positive and half negative, from random positions with a given energy |

Each script adds the conformal-map folder it needs to the MATLAB path
(`addpath('../Conformal Heart/')` or `addpath('../Conformal Gap/')`). Both
folders define functions with the same names (`DF1`–`DF4`, `inverse_map`,
`jacob`), and whichever folder was added last is used. Always run the scripts
from their own folder, so that this works and so the model's own functions are
found first.

### Output and restarting

The state, energy and (for the disc) angular momentum after every time step
`dt` are appended to a NetCDF file. The file is written in the folder *above*
the script's folder, with the script's folder name as a prefix:

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
  and doubly periodic runs;
- $|X|^2 = 25$, an aspect ratio of about 27 (`ElongationEvent2.m`), in the
  gap runs.

The plots at the end of each script are switched on and off with
`if 1==1` / `if 1==0`. The movie options save `Vortices.mp4` in the current
folder.

### Contour-dynamics comparison

`Contour_Dynamics_and_MZS.m` reads a contour-dynamics simulation of the same
vortex patch from a text file `fort.11` in the current folder. **The
contour-dynamics code that produces `fort.11` is not part of this
repository.**

The script asks for plot settings at the keyboard. It saves each frame as a
PNG in `figfiles/` and can also save a movie, `CD & MZS.mp4`. It then plots
time series of the centroid, aspect ratio and rotation rate from both models.
It assumes that frame $k$ of `fort.11` is at the same time as MZS step $k$.

## Conformal maps

| Map | $z = F(Z)$ | Parameters |
|---|---|---|
| Heart | $\dfrac{a\,Z}{(1 - q^2 Z^2)(1 - i g Z)}$ | `domain_var = [q; g]`. $a$ = `alpha_heart(q, g, pi)` makes the domain's area $\pi$. `q = g = 0` gives the unit disc. |
| Gap | $\dfrac{2iZ}{Z^2 - 1}$ | None. The walls are $x \le -1$ and $x \ge 1$ on the real axis. The `domain_var` argument is accepted but not used. |

Random initial positions in the heart (`set_ic_bounded_MZS.m`,
`set_ic_mc_MZS.m`) use `boundary.m`, so they only work with the Conformal
Heart folder. `set_ic_boundedL_MZS.m` is for the unit disc only.

## Requirements

- MATLAB R2022a or later. The scripts use `Name=Value` arguments (R2021a)
  and `fontsize` (R2022a).
- No toolboxes. `ode113`, the NetCDF functions and `VideoWriter` are all part
  of base MATLAB.
