# Radial Teukolsky module

## Background

The `RadialTeukolsky` class solves the homogeneous radial Teukolsky equation
$$\left[\Delta^{-s} \frac{d}{dr} \left(\Delta^{s+1} \frac{d }{dr}  \right) + \left(\frac{K^2-2is(r-M)K}{\Delta}+4is\omega r - \lambda_{sjm\omega} \right)\right]R_{sjm\omega} =0$$
given the parameters:
- $s$ : spin-weight of the perturbation
- $j$ : the spheroidal polar mode number
- $m$ : the azimuthal mode number
- $a$ : the dimensionless Kerr spin parameter
- $\omega$ : the frequency

where $\Delta=r^2-2Mr+a^2$, $K=(r^2+a^2)\omega-ma$, and $\lambda_{sjm\omega}$ is the spheroidal eigenvalue (separation constant).


## API

```{eval-rst}
.. automodule:: pybhpt.radial
   :no-index:
   :members:
   :undoc-members:
   :show-inheritance:
```
