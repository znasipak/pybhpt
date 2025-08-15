# Angular Teukolsky equation

The spin-weighted spheroidal harmonics $S_{sjm\gamma}(\theta)$ satisfy the equation
$$
	\left[\frac{1}{\sin\theta}\frac{d}{d\theta}\left(\sin\theta \frac{d}{d\theta} \right)
	- \left(\gamma^2\sin^2\theta+\frac{(m+s\cos\theta)^2}{\sin^2\theta}
	+2\gamma s\cos\theta-s-2m\gamma-\lambda_{sjm\gamma} \right)\right]S_{sjm\gamma} = 0,
$$
with eigenvalues $\lambda_{sjm\gamma}$ and spheroidicity $\gamma$. For the Teukolsky equation, $\gamma = a\omega$.

The spheroidal harmonics can be expressed as rapidly convergent sums of spin-weighted spherical harmonics,
$$
    S_{sjm\gamma}(\theta)e^{im\phi} = \sum_{\ell = \ell_\mathrm{min}}^\infty b^{\ell}_{sjm\gamma} Y_{s\ell m}(\theta,\phi),
$$
where $\ell_\mathrm{min} = \mathrm{max}[|m|, |s|]$ and the coupling coefficients $b^{\ell}_{sjm\gamma}$ satisfy a five-term recursion relation. Similarly we can represent spin-weighted spherical harmonics as finite series of *scalar* spherical harmonics $Y_{lm}=Y_{0lm}$
$$
    Y_{s\ell m}(\theta)
    = \sin^{-|s|}\theta \sum_{l=|m|}^\infty \mathcal{A}^l_{s\ell m} 
    Y_{lm}(\theta),
$$
where the coefficients are given analytically in terms of the Wigner $3j$ symbol,
$$\begin{align}
	\mathcal{A}^l_{s\ell m} &= 
	(-1)^{m+s(1+\mathrm{sgn}(s))/2} \mathcal{C}_{|s|\ell l}
		\left(
			\begin{array}{ccc}
				|s| & \ell & l
				\\
				0 & m & -m
			\end{array}
		\right)
		\left(
			\begin{array}{ccc}
				|s| & \ell & l
				\\
				s & -s & 0
			\end{array}
		\right),
  \\
  \mathcal{C}_{s\ell l}  &= \sqrt{\frac{4^{s} (s!)^2 (2\ell+1)(2l+1)}{(2s)!}}.
\end{align}$$
The selection rules of the $3j$-symbol mean that the coefficients vanish unless $\ell - |s| \leq l \leq \ell + |s|$.