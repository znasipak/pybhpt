# Flux module

## API

The time-averaged rate of change of the orbital energy $\langle \dot{E}\rangle$, angular momentum $\langle \dot{L}_z\rangle$, and Carter constant $\langle \dot{Q}\rangle$ can be expressed in terms of the Teukolsky amplitudes $Z^\mathrm{Up/In}_{sjmkn}$,
$$ 
\begin{align}
\langle \dot{\mathcal{J}} \rangle & = \langle \dot{\mathcal{J}} \rangle^\mathrm{inf} + \langle \dot{\mathcal{J}} \rangle^\mathrm{hor},
\\
&= \sum_{jmkn} \langle \dot{\mathcal{J}} \rangle^\mathrm{inf/hor}_{jmkn},
\\ 
&= \sum_{jmkn} \alpha_{sjmkn}^{(\mathcal{J})\mathrm{inf/hor}}\left| Z^\mathrm{Up/In}_{sjmkn} \right|^2,
\end{align}$$
where $\mathcal{J} = (E, L_z, Q)$, $\alpha_{sjmkn}^{(\mathcal{J})\mathrm{inf}/\mathrm{hor}}$ are known coefficients, and $s=0, \pm 2$ produce either scalar or gravitational fluxes, respectively. 

The `FluxMode` class takes as input an instances of the Kerr geodesic class and the Teukolsky class for a mode $(s,j,m,k,n)$, and produces the flux contribution for that given mode.

```{eval-rst}
.. automodule:: pybhpt.flux
   :no-index:
   :members:
   :undoc-members:
   :show-inheritance:
```
