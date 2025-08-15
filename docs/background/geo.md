# Geodesics in Kerr

We work in a Kerr background with angular momentum and mass parameters $(J, M)$ and use Boyer-Lindquist coordinates $(t, r, \theta, \phi)$.

Bound periodic timelike geodesics in Kerr spacetime are defined in terms of the turning points of the motion:

- $r_\mathrm{min}$ : minimum Boyer-Lindquist radius
- $r_\mathrm{max}$ : maximum Boyer-Lindquist radius
- $\theta_\mathrm{min}$ : minimum Boyer-Lindquist polar angle
- $\pi-\theta_\mathrm{min}$ : maximum Boyer-Lindquist polar angle

From these we parametrize the geodesic in terms of the generalized Keplerian parameters:

- $a$ : the dimensionless Kerr spin parameter
- $p$ : the dimensionless semilatus rectum
- $e$ : the orbital eccentricty
- $x$ : cosine of the orbital inclination

where $a = J/M^2$, $pM = 2r_\mathrm{max}r_\mathrm{min}/(r_\mathrm{min}+r_\mathrm{max})$ and $e = (r_\mathrm{max}-r_\mathrm{min})/(r_\mathrm{min}+r_\mathrm{max})$. The motion can also be described in terms of the conserved orbital constants:

- $E$ : the specific orbital energy
- $L_z$ : the specific orbital angular momentum
- $Q$ : the Carter constant

which have units $\{1, M, M^2\}$, respectively, along with the mass of the small body $\mu$.

With these conserved quantities, we obtain four first-order ordinary differential equations (ODEs) for $x_p^\mu$, which decouple when parametrized in terms of the Mino(-Carter) time parameter $\lambda$,
$$\begin{align}
    \frac{dt_p}{d\lambda} &= V_{tr}(r_p) + V_{t\theta}(\theta_p),
    \\
    \frac{dr_p}{d\lambda} &= \pm \sqrt{V_r(r_p)},
    \\ 
    \frac{d\theta_p}{d\lambda} &= \pm \sqrt{V_\theta(\theta_p)},
    \\
    \frac{d\phi_p}{d\lambda} &= V_{\phi r}(r_p) + V_{\phi \theta}(\theta_p),
\end{align}$$
where $d\lambda = \Sigma^{-1} d\tau$, and the potential functions are given by
$$\begin{align} 
    V_{r}(r) &= { P^2(r) - \Delta\left(r^2 + {K} \right),} 
    &
    V_{\theta}(\theta) &= {{Q} - {L}_z^2 \cot^2\theta - a^2 (1 -{E}^2)\cos^2\theta,}
    \\
    V_{tr}(r) &= \frac{r^2+a^2}{\Delta}P(r), 
    &
    V_{t\theta}(\theta) &= a{L}_z - a^2 {E} \sin^2\theta, 
    \\
    V_{\phi r}(r) &= \frac{a}{\Delta}P(r), 
    &
    V_{\phi \theta}(\theta) &= {L}_z \csc^2\theta - a {E},
\end{align}$$

with $P(r) = (r^2+a^2){E} - a {L}_z$.

The resulting bound solutions can be separated into terms that are periodic with respect to the Mino time radial and polar frequencies $\Upsilon_r$ and $\Upsilon_\theta$ and terms that grow secularly with the Mino time rates $\Upsilon_t$ and $\Upsilon_\phi$.
Therefore, the fundamental coordinate time frequencies are given by
$$\begin{align}
    \Omega_r &= \frac{\Upsilon_r}{\Upsilon_t},
    &
    \Omega_\theta &= \frac{\Upsilon_\theta}{\Upsilon_t},
    &
    \Omega_\phi &= \frac{\Upsilon_\phi}{\Upsilon_t}.
\end{align}$$

We then represent the radial and polar motion by
$$\begin{align}
    r_p(\lambda) &= \Delta r^{(r)}(\Upsilon_r\lambda) = \Delta r^{(r)}(\Upsilon_r\lambda + 2\pi),
    \\
    \theta_p(\lambda) &= \Delta \theta^{(\theta)}(\Upsilon_\theta\lambda) = \Delta \theta^{(\theta)}(\Upsilon_\theta\lambda + 2\pi),
\end{align}$$
while time and azimuthal angle grow as
$$\begin{align}
    t_p(\lambda) &= \Upsilon_t \lambda + \Delta t^{(r)}(\Upsilon_r\lambda) + \Delta t^{(\theta)}(\Upsilon_\theta\lambda),
    \\ 
    \phi_p(\lambda) &= \Upsilon_\phi \lambda + \Delta \phi^{(r)}(\Upsilon_r\lambda) + \Delta \phi^{(\theta)}(\Upsilon_\theta\lambda),
\end{align}$$
where $\Delta t^{(r)}$, $\Delta \phi^{(r)}$, $\Delta t^{(\theta)}$, and $\Delta \phi^{(\theta)}$ are $2\pi$-periodic odd functions.