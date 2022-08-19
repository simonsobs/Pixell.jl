import PhysicalConstants.CODATA2018 as units
import Unitful: ustrip

@doc raw"""
The derivative of the planck spectrum with respect to temperature, evaluated
at frequencies f and temperature T, in units of Jy/sr/K.

A blackbody has intensity ``I = 2hf^3/c^2/(\exp{hf/kT}-1) = V/(\exp{x}-1)``
with ``V = 2hf^3/c^2``, ``x = hf/kT``.
``dI/dx = -V/(\exp{x}-1)^2 * \exp{x}``
``dI/dT = dI/dx * dx/dT``
``      = 2hf^3/c^2/(\exp{x}-1)^2*\exp{x} * hf/k / T^2``
``      = 2*h^2*f^4/c^2/k/T^2 * \exp{x}/(\exp{x}-1)^2``
``      = 2*x^4 * k^3*T^2/(h^2*c^2) * \exp{x}/(\exp{x}-1)^2``
``      = .... /(4*sinh(x/2)^2)``
"""
function dplanck(f, T=2.72548)
    c = ustrip(units.c_0)
    k = ustrip(units.k_B)
    h = ustrip(units.h)

    x = h*f/(k*T)
    dIdT = 2*x^4 * k^3*T^2/(h^2*c^2) / (4*sinh(x/2)^2) * 1e26
    return dIdT
end