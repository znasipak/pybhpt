# Teukolsky equation

In Boyer-Lindquist coordinates $(t,r,\theta,\phi)$, the Teukolsky equation takes the form
$$\begin{align}
    \mathcal{O}^T_s \Psi_s = -4\pi \zeta\bar{\zeta} T_s,
\end{align}$$
where $\zeta=r-ia\cos\theta$, $T_s$ is the (tetrad-projected) source, and the spin-weighted Teukolsky operator is given by,
$$\begin{align}
    \mathcal{O}^T_s &:=-\left(\frac{(r^2+a^2)^2}{\Delta} - a^2 \sin^2\theta\right)\partial_t^2 - \frac{4 M a r}{\Delta} \partial_t \partial_\phi - \left(\frac{a^2}{\Delta} - \frac{1}{\sin^2\theta}\right)\partial_\phi^2
    \\ 
    & \qquad \qquad + \Delta^{-s}\partial_r\left(\Delta^{s+1} \partial_r\right) + \frac{1}{\sin\theta}\partial_\theta\left({\sin\theta}\partial_\theta \right) + 2s\left(\frac{M(r^2-a^2)}{\Delta} - r-ia\cos\theta \right)\partial_t
    \\ \notag
    & \qquad \qquad \qquad \qquad +2 s\left(\frac{a(r-M)}{\Delta}  - \frac{i\cos\theta}{\sin^2\theta} \right)\partial_\phi - s\left(\frac{s\cos^2\theta}{\sin^2\theta} - 1\right),
\end{align}$$
with $\Delta = r^2-2Mr+a^2$, and $\Psi_{0}=\Phi$ for scalar perturbations; $\Psi_{+1}=\phi_0$ and $\Psi_{-1}=\zeta^2\phi_2$ ; and $\Psi_{+2}=\psi_0$ and $\Psi_{-2}=\zeta^4\psi_4$ for gravitational perturbations.

Transforming to the frequency domain, the fields and sources can be mode-decomposed into
$$\begin{align}
    \Psi_s = \int d\omega \sum_{jm} \Psi_{sjm\omega}(r) S_{sjm\gamma}(\theta) e^{im\phi} e^{-i\omega t},
    \\
    T_s = \int d\omega \sum_{jm} T_{sjm\omega}(r) S_{sjm\gamma}(\theta) e^{im\phi} e^{-i\omega t},
\end{align}$$
with $\gamma = a\omega$, leading to separated ODEs
$$\left[\Delta^{-s} \frac{d}{dr} \left(\Delta^{s+1} \frac{d }{dr}  \right) + \left(\frac{K^2-2is(r-M)K}{\Delta}+4is\omega r - \lambda_{sjm\omega} \right)\right]R_{sjm\omega} = T_{slm\omega},$$
and
$$
	\left[\frac{1}{\sin\theta}\frac{d}{d\theta}\left(\sin\theta \frac{d}{d\theta} \right)
	- \left(\gamma^2\sin^2\theta+\frac{(m+s\cos\theta)^2}{\sin^2\theta}
	+2\gamma s\cos\theta-s-2m\gamma-\lambda_{sjm\gamma} \right)\right]S_{sjm\gamma} = 0,
$$
with eigenvalues $\lambda_{sjm\gamma}$ and spheroidicity $\gamma$.

The time-averaged rate of change of the orbital energy $\langle \dot{E}\rangle$, angular momentum $\langle \dot{L}_z\rangle$, and Carter constant $\langle \dot{Q}\rangle$ can then be expressed in terms of the Teukolsky amplitudes $Z^\mathrm{Up/In}_{sjmkn}$,
$$ 
\begin{align}
\langle \dot{\mathcal{J}} \rangle & = \langle \dot{\mathcal{J}} \rangle^\mathrm{inf} + \langle \dot{\mathcal{J}} \rangle^\mathrm{hor},
\\
&= \sum_{jmkn} \langle \dot{\mathcal{J}} \rangle^\mathrm{inf/hor}_{jmkn},
\\ 
&= \sum_{jmkn} \alpha_{sjmkn}^{(\mathcal{J})\mathrm{inf/hor}}\left| Z^\mathrm{Up/In}_{sjmkn} \right|^2,
\end{align}$$
where $\mathcal{J} = (E, L_z, Q)$, $\alpha_{sjmkn}^{(\mathcal{J})\mathrm{inf}/\mathrm{hor}}$ are known coefficients, and $s=0, \pm 2$ produce either scalar or gravitational fluxes, respectively. 