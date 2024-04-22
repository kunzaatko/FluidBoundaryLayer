# Physical Background

This repository contains a `C` implementation of the _shooting method_ numerical solution of the
[_Falkner-Skan_](https://en.wikipedia.org/wiki/Falkner–Skan_boundary_layer) 
generalization to the 
[_Blasius boundary layer_](https://en.wikipedia.org/wiki/Blasius_boundary_layer)
differential equation of the form:
```math
y''' + y y'' + \lambda \left( 1 - y'^{2} \right)  = 0, \text{ for } x \in (A, B), \lambda \in \langle 0, 0.5 \rangle
```
and the boundary conditions
```math
 y(A) = y'(A) = 0, y'(B) = 1.
```

This differential equation characterises the flow of a viscous fluid over a flat plate. The Falkner-Skan generalization is in  the fact that the angle of attack of the fluid is not strictly fixed to $0$ as in the Blasius formulation, but forms a wedge with the angle $\lambda \pi$.
<p align="center">
    <img src="static/wedge.png" alt="A wedge over which a flow of fluid is illustrated by two arrows copying the boundary" width="250" />
</p>

The differential equation describes the velocity of flow $u$ around a wedge assuming that the external uniform flow velocity
is $U$. The unknown variable in the equation $y$ represents $y' = y_x = u/U$. Where $x$ is the distance along the normal
vector to the wedge of the position in of the flow velocity $u$. Since a picture is worth a thousand words, here is
a diagram that illustrates the relation of the variables that act in the equation:

<p align="center">
    <img src="static/flow_profile.png" alt="A diagram that illustrates the relation of the variables that act in the equation Falkner-Skan boundary layer differential equation" width="600" />
</p>

It is developed from the _Navier-Stokes_ equations by applying boundary layer approximations, which stem in neglecting
the terms which are small in the boundary layer.
The value of $y$ represents the width of the boundary layer in relation to the distance from the point of first
contact. 
A boundary layer is defined here as the width of the fluid that is affected by the viscosity.

# Contents

In relation to the above, this repo contains:

- `C` implementation (`src/`) with tests (`test/`)
- $\LaTeX$ report of the internals of the numerical algorithm (`report/`)
- `Julia` script for reading the `csv` data and generating a plot of $y$, $y'$ and $y''$ in (`plotting/`)


## Algorithm Implementation 

In the `src/` directory, the `C` code is located. The entry point `src/main.c` holds the parameters that affect the
    starting interval, where the root is to be found by the shooting and it that must be set manually in and must 
    contain the correct value in it, i.e. 

$$y''(A)_\text{true}\in (\texttt{zA2start\\_1}, \texttt{zA2start\\_2}).$$

All the other parameters that are used for the integration of the DE and the parameters to the equation itself are
    stored as `const` values in the file `src/blias.h`.
By default the parameters are set to $(A, B)  = (0, 10)$ and $\lambda = 0.0$ (corresponding to the Blasius boundary 
    layer).
When you are satisfied with setting the parameters, from the root directory you can compile and run the program by

```bash
$ just build-program # or `make program`
$ ./program
```
This creates the file `filename` (`lam_0_0.csv` by default) containing the data.
