from cybhpt_full import RadialTeukolsky as RadialTeukolskyCython
from cybhpt_full import available_methods as available_methods_cython
from cybhpt_full import renormalized_angular_momentum as nu_cython
from cybhpt_full import renormalized_angular_momentum_monodromy as nu_2_cython
from cybhpt_full import hypergeo_2F1 as hypergeo_2F1_cython
import numpy as np

def available_methods():
    """
    Returns a list of available solution methods.
    """
    return available_methods_cython()

def renormalized_angular_momentum(s, j, m, a, omega):
    """
    Computes the renormalized angular momentum for the given parameters.

    Parameters
    ----------
    s : int
        The spin weight of the field.
    j : int
        The spheroidal harmonic mode number.
    m : int
        The azimuthal harmonic mode number.
    a : float
        The black hole spin parameter.
    omega : float
        The frequency of the mode.  

    Returns
    -------
    complex
        The renormalized angular momentum.
    """
    return nu_cython(s, j, m, a, omega)

def renormalized_angular_momentum_monodromy(s, j, m, a, omega, la):
    """
    Computes the renormalized angular momentum using the monodromy method for the given parameters.

    Parameters
    ----------
    s : int
        The spin weight of the field.
    j : int
        The spheroidal harmonic mode number.
    m : int
        The azimuthal harmonic mode number.
    a : float
        The black hole spin parameter.
    omega : complex
        The frequency of the mode.
    la : complex
        The spheroidal eigenvalue.

    Returns
    -------
    complex
        The renormalized angular momentum.
    """
    return nu_2_cython(s, j, m, a, omega, la)

class RadialTeukolsky:
    """A class for solving the homogeneous radial Teukolsky equation.

    Parameters
    ----------
    s : int
        The spin weight of the field
    j : int
        The spheroidal harmonic mode number
    m : int
        The azimuthal harmonic mode number
    a : float
        The black hole spin parameter
    omega : float
        The frequency of the mode
    r : numpy.ndarray
        A numpy array of radial points at which to evaluate the solution

    Attributes
    ----------
    radialpoints : numpy.ndarray
        A numpy array of radial points at which the solution is evaluated.
    base : RadialTeukolskyCython
        The underlying Cython object that performs the computations.
    nsamples : int
        The number of radial points in the radialpoints array.
    
    Properties
    ----------
    blackholespin : float
        The black hole spin parameter.
    spinweight : int
        The spin weight of the field.
    s : int
        Alias for spinweight.
    spheroidalmode : int
        The spheroidal harmonic mode number.
    j : int
        Alias for spheroidalmode.
    azimuthalmode : int
        The azimuthal harmonic mode number.
    m : int
        Alias for azimuthalmode.
    frequency : float
        The frequency of the mode.
    mode_frequency : float
        Alias for frequency.
    omega : float
        Alias for frequency.
    eigenvalue : float
        The spheroidal eigenvalue of the radial Teukolsky equation.
    
    Methods
    -------
    solveboundarycondition(method)
        Solves the boundary condition for the radial Teukolsky equation.
    setboundarycondition(bc, R, Rp, r)
        Sets the boundary condition for the radial Teukolsky equation.
    solve(method="AUTO", bc=None)
        Solves the radial Teukolsky equation.
    flipspinweight()
        Flips the spin weight of the field.
    radialpoint(pos)
        Returns the radial point at the given position.
    boundarypoint(bc)
        Returns the boundary point for the given boundary condition.
    boundarysolution(bc)
        Returns the solution at the boundary for the given boundary condition.
    boundaryderivative(bc)
        Returns the derivative at the boundary for the given boundary condition.
    radialsolution(bc, pos)
        Returns the solution at the radial point for the given boundary condition and position.
    radialderivative(bc, pos)
        Returns the derivative at the radial point for the given boundary condition and position.
    radialderivative2(bc, pos)
        Returns the second derivative at the radial point for the given boundary condition and position.
    radialsolutions(bc)
        Returns the solutions at all radial points for the given boundary condition.
    radialderivatives(bc)
        Returns the derivatives at all radial points for the given boundary condition.
    radialderivatives2(bc)
        Returns the second derivatives at all radial points for the given boundary condition.
    __call__(bc, deriv=0)
        Returns the solutions, first derivatives, or second derivatives at all radial points for the given boundary condition.
        The `deriv` parameter specifies which derivative to return: 0 for solutions,
        1 for first derivatives, and 2 for second derivatives. If `deriv` is not 0, 1, or 2, a ValueError is raised.
    """
    def __init__(self, s, j, m, a, omega, r):
        if a < 0 or a > 1:
            raise ValueError(f"Black hole spin parameter {a} must be in the range [0, 1].")
        if j < np.abs(m):
            raise ValueError(f"Spheroidal harmonic mode number {j} must be greater than or equal to the absolute value of azimuthal harmonic mode number {m}.")
        if np.any(r <= 1 + np.sqrt(1 - a**2)):
            raise ValueError(f"Radial point {r} must be greater than horizon radius r_+ = {1 + np.sqrt(1 - a**2)}.")
        if isinstance(r, list) or (isinstance(r, np.ndarray) and r.ndim > 0):
            self.radialpoints = np.asarray(r)
            self.nsamples = self.radialpoints.shape[0]
        else:
            raise AttributeError("Radial points must be a list or a numpy array.")
        
        if self.nsamples == 0:
            raise ValueError("Radial points array is empty.")
        self.base = RadialTeukolskyCython(a, s, j, m, omega, self.radialpoints)


    @property
    def blackholespin(self):
        return self.base.blackholespin
    
    @property
    def spinweight(self):
        return self.base.spinweight
    
    @property
    def s(self):
        return self.spinweight

    @property
    def spheroidalmode(self):
        return self.base.spheroidalmode

    @property
    def j(self):
        return self.spheroidalmode

    @property
    def azimuthalmode(self):
        return self.base.azimuthalmode
    
    @property
    def m(self):
        return self.azimuthalmode

    @property
    def frequency(self):
        return self.base.frequency
    
    @property
    def mode_frequency(self):
        return self.frequency

    @property
    def omega(self):
        return self.frequency

    @property
    def eigenvalue(self):
        return self.base.eigenvalue
    
    def solveboundarycondition(self, method):
        """Solves the boundary condition for the radial Teukolsky equation.

        Parameters
        ----------
        method : str
            The method to use for solving the boundary condition. Default is "AUTO".
        """
        self.base.solve_bc(method)

    def setboundarycondition(self, bc, R, Rp, r):
        """Sets the boundary condition for the radial Teukolsky equation.

        Parameters
        ----------
        bc : str
            The boundary condition to set. Can be "In" for horizon or "Up" for infinity.
        R : float
            The boundary condition function value at the radial point.
        Rp : float
            The derivative of the boundary condition function at the radial point.
        r : float
            The radial point at which the boundary condition is defined.
        """
        self.base.set_bc(bc, R, Rp, r)

    def solve(self, method = "AUTO", bc = None):
        """Solves the radial Teukolsky equation.

        Parameters
        ----------
        method : str, optional
            The method to use for solving the equation. Default is "AUTO".
        bc : str, optional
            Specifies which homogeneous solutions to compute. If None, both "In" (horizon) and "Up" (infinity) solutions are computed.
            If "In", only the horizon solution is computed. If "Up", only the infinity solution is computed.
        """
        if bc is None:
            self.base.solve(method, "None")
        else:
            self.base.solve(method, bc)

    def flipspinweight(self):
        """Flips the sign of the spin weight of the field."""
        self.base.flip_spinweight()

    def radialpoint(self, pos):
        """Returns the radial point at the given position.
        Parameters
        ----------
        pos : float
            The position at which to evaluate the radial point.
        Returns
        -------
        float
            The radial point at the given position.
        """
        return self.base.radialpoint(pos)
    
    def boundarypoint(self, bc):
        """
        Returns the boundary point for the given boundary condition.
        Parameters
        ----------
        bc : str
            The boundary condition to evaluate. Can be "In" for horizon or "Up" for infinity.
        Returns
        -------
        float
            The boundary point corresponding to the specified boundary condition.
        """
        return self.base.boundarypoint(bc)
    
    def boundarysolution(self, bc):
        """Returns the solution at the boundary for the given boundary condition.
        Parameters
        ---------- 
        bc : str
            The boundary condition to evaluate. Can be "In" for horizon or "Up" for infinity.
        Returns
        -------
        float
            The solution at the boundary corresponding to the specified boundary condition.
        """
        return self.base.boundarysolution(bc)
    
    def boundaryderivative(self, bc):
        """Returns the derivative at the boundary for the given boundary condition.
        Parameters
        ----------
        bc : str
            The boundary condition to evaluate. Can be "In" for horizon or "Up" for infinity.
        Returns
        -------
        float
            The derivative at the boundary corresponding to the specified boundary condition.
        """
        return self.base.boundaryderivative(bc)
    
    def radialsolution(self, bc, pos):
        """Returns the solution at the radial point for the given homogeneous solution and position.
        Parameters
        ----------
        bc : str
            The homogeneous solution to evaluate. Can be "In" for horizon solution or "Up" for infinity solution.
        pos : float
            The position at which to evaluate the solution.
        Returns
        -------
        float
            The solution at the radial point corresponding to the specified boundary condition and position.
        """
        return self.base.solution(bc, pos)
    
    def radialderivative(self, bc, pos):
        """Returns the derivative at the radial point for the given homogeneous solution and position.
        Parameters
        ----------
        bc : str
            The homogeneous solution to evaluate. Can be "In" for horizon solution or "Up" for infinity solution.
        pos : float
            The position at which to evaluate the solution.
        Returns
        -------
        float
            The solution at the radial point corresponding to the specified boundary condition and position.
        """
        return self.base.derivative(bc, pos)
    
    def radialderivative2(self, bc, pos):
        """Returns the second derivative at the radial point for the given homogeneous solution and position.
        Parameters
        ----------
        bc : str
            The homogeneous solution to evaluate. Can be "In" for horizon solution or "Up" for infinity solution.
        pos : float
            The position at which to evaluate the solution.
        Returns
        -------
        float
            The solution at the radial point corresponding to the specified boundary condition and position.
        """
        return self.base.derivative2(bc, pos)

    def radialsolutions(self, bc):
        """Returns a homogeneous solution at all radial points.
        Parameters
        ----------
        bc : str
            The homogeneous solution to evaluate. Can be "In" for horizon solution or "Up" for infinity solution.
        Returns
        -------
        numpy.ndarray
            A numpy array of solutions at all radial points corresponding to the specified boundary condition.
        """
        return np.array([self.base.solution(bc, i) for i in range(self.nsamples)])
    
    def radialderivatives(self, bc):
        """Returns the radial derivative at all radial points.
        Parameters
        ----------
        bc : str
            The homogeneous solution to evaluate. Can be "In" for horizon solution or "Up" for infinity solution.
        Returns
        -------
        numpy.ndarray
            A numpy array of solutions at all radial points corresponding to the specified boundary condition.
        """
        return np.array([self.base.derivative(bc, i) for i in range(self.nsamples)])
    
    def radialderivatives2(self, bc):
        """Returns the second radial derivative at all radial points.
        Parameters
        ----------
        bc : str
            The homogeneous solution to evaluate. Can be "In" for horizon solution or "Up" for infinity solution.
        Returns
        -------
        numpy.ndarray
            A numpy array of solutions at all radial points corresponding to the specified boundary condition.
        """
        return np.array([self.base.derivative2(bc, i) for i in range(self.nsamples)])
    
    def __call__(self, bc, deriv = 0):
        """Returns the solutions, first derivatives, or second derivatives at all radial points for the given boundary condition.

        Parameters
        ----------
        bc : str
            The homogeneous solution to evaluate. Can be "In" for horizon solution or "Up" for infinity solution.
        deriv : int, optional
            Specifies which derivative to return: 0 for solutions, 1 for first derivatives, and 2 for second derivatives.
            Default is 0 (solutions).
        Returns
        -------
        numpy.ndarray
            A numpy array of solutions, first derivatives, or second derivatives at all radial points corresponding to the specified boundary condition.
        Raises
        ------
        ValueError
            If `deriv` is not 0, 1, or 2.
        """ 
        if deriv == 0:
            return self.radialsolutions(bc)
        elif deriv == 1:
            return self.radialderivatives(bc)
        elif deriv == 2:
            return self.radialderivatives2(bc)
        else:
            raise ValueError("RadialTeukolsky only solves up to the second derivative")
        
def hypergeo_2F1(a, b, c, x):
    """
    Gauss hypergeometric function 2F1(a, b; c; x). Note that this function is not very stable across the complex domain.
    
    Parameters
    ----------
    a : complex
        The first parameter of the hypergeometric function.
    b : complex
        The second parameter of the hypergeometric function.
    c : complex
        The third parameter of the hypergeometric function.
    x : complex
        The argument of the hypergeometric function.

    Returns
    -------
    complex
        The value of the hypergeometric function 2F1(a, b; c; x).
    """
    return hypergeo_2F1_cython(a, b, c, x)