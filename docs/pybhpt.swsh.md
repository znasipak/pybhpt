# Spin-weighted spheroidal harmonic module

The spin-weighted spheroidal harmonics $S_{sjm\gamma}(\theta)$ satisfy the equation
$$
	\left[\frac{1}{\sin\theta}\frac{d}{d\theta}\left(\sin\theta \frac{d}{d\theta} \right)
	- \left(\gamma^2\sin^2\theta+\frac{(m+s\cos\theta)^2}{\sin^2\theta}
	+2\gamma s\cos\theta-s-2m\gamma-\lambda_{sjm\gamma} \right)\right]S_{sjm\gamma} = 0,
$$
with eigenvalues $\lambda_{sjm\gamma}$ and spheroidicity $\gamma$. For the Teukolsky equation, $\gamma = a\omega$.


## API

```{eval-rst}
.. automodule:: pybhpt.swsh
   :no-index:
   :members:
   :undoc-members:
   :show-inheritance:
```
