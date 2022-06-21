const h = 6.62606957e-34
const k  = 1.3806488e-23
const T_cmb = 2.72548 # +/- 0.00057
const c  = 299792458.0

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
function dplanck(f, T=T_cmb)
	x = h*f/(k*T)
	dIdT = 2*x^4 * k^3*T^2/(h^2*c^2) / (4*sinh(x/2)^2) * 1e26
	return dIdT
end