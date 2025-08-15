# Radial Teukolsky equation

The radial Teukolsky equation for spin-weight $s$, harmonic numbers $(j,m)$, frequency $\omega$, and $a$ the dimensionless Kerr spin parameter is given by
$$\left[\Delta^{-s} \frac{d}{dr} \left(\Delta^{s+1} \frac{d }{dr}  \right) + \left(\frac{K^2-2is(r-M)K}{\Delta}+4is\omega r - \lambda_{sjm\omega} \right)\right]R_{sjm\omega} = T_{slm\omega},$$
where $\Delta=r^2-2Mr+a^2$, $K=(r^2+a^2)\omega-ma$, $\lambda_{sjm\omega}$ is the spheroidal eigenvalue (separation constant), and $T_{slm\omega}$ is the radial decomposition of the source.

## Homogeneous solutions

For $T_{slm\omega} = 0$, we construct the homogeneous solutions $R^\mathrm{In}_{sjm\omega}$ and $R^\mathrm{Up}_{sjm\omega}$, which correspond to the asymptotic boundary conditions
$$\begin{align}
    R^\mathrm{In}_{sjm\omega} (r \rightarrow r_+) &\sim A^\mathrm{trans}_s \Delta^{-s} e^{-i k r_*},
    \\
    R^\mathrm{In}_{sjm\omega} (r \rightarrow \infty) &\sim A^\mathrm{ref}_s r^{-(2s+1)} e^{i\omega r_*}
    + A^\mathrm{inc}_s r^{-1} e^{-i\omega r_*},
    \\
    R^\mathrm{Up}_{sjm\omega} (r \rightarrow r_+) &\sim B^\mathrm{ref}_s \Delta^{-s} e^{-i k r_*}
    + B^\mathrm{inc}_s e^{i k r_*},
    \\
    R^\mathrm{Up}_{sjm\omega} (r \rightarrow \infty) &\sim 
    B^\mathrm{trans}_s r^{-(2s+1)} e^{i\omega r_*},
\end{align}$$
where $k = \omega - m \omega_+$, $\omega_+ = a/(2Mr_+)$, and the tortoise coordinate is given by the differential relation $dr_{*}/dr = (r^2+a^2)/\Delta$, which can be immediately integrated, leading to
$$r_* = r + \frac{r_+}{\kappa} \ln \frac{r-r_+}{2M} - \frac{r_-}{\kappa} \ln \frac{r-r_-}{2M},$$
where $\kappa = \sqrt{1 - q^2}$ and $q = a/M$.

## Inhomogeneous solutions

Rather than constructing the full inhomogeneous solutions, we can instead construct the so-called *extended homogeneous solutions* for a point-particle source on a bound periodic geodesic,
$$\begin{align}
    \Psi_s &= \Psi_s^+ \Theta(r-r_p) + \Psi_s^- \Theta(r_p-r),
    \\
    \Psi_s^\pm &= \sum_{jmkn}\Psi^\pm_{sjmkn}(r)S_{sj\gamma_{mkn}}(\theta)e^{im\phi}e^{-i\omega_{mkn}t},
\end{align}$$
where we have the mode frequencies $\omega_{mkn} = m\Omega_\phi + k \Omega_\theta + n \Omega_r$, the discrete spheroidicity $\gamma_{mkn} = a\omega_{mkn}$, and the extended homogeneous radial solutions
$$
    \Psi^\pm_{sjmkn}(r) = Z^{\mathrm{Up/In}}_{sjmkn} R^{\mathrm{Up/In}}_{sjmkn}(r),
$$
with Teukolsky amplitudes $Z^{\mathrm{Up/In}}_{sjmkn}$.