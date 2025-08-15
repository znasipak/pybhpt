# Teukolsky point-particle mode module

The `TeukolskyMode` class constructs modes of the so-called extended homogeneous solutions to the radial Teukolsky equation for a point-particle source on a bound periodic geodesic,
$$\begin{align}
    \Psi_s &= \Psi_s^+ \Theta(r-r_p) + \Psi_s^- \Theta(r_p-r),
    \\
    \Psi_s^\pm &= \sum_{jmkn}Z^{\mathrm{Up/In}}_{sjmkn}R^{\mathrm{Up/In}}_{sjmkn}(r)S_{sj\gamma_{mkn}}(\theta)e^{im\phi}e^{-i\omega_{mkn}t},
\end{align}$$
where we have the mode frequencies $\omega_{mkn} = m\Omega_\phi + k \Omega_\theta + n \Omega_r$, the discrete spheroidicity $\gamma_{mkn} = a\omega_{mkn}$.

## API

```{eval-rst}
.. automodule:: pybhpt.teuk
   :no-index:
   :members:
   :undoc-members:
   :show-inheritance:
```
