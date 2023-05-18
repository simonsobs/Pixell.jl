
```@meta
CurrentModule = Pixell
```

# Maps and Pixels

## Rectangular Pixelizations

Pixell allows you to work with rectangular pixelizations, with an emphasis on equidistant cylindrical projections (ECP) like [plate carrée](https://pro.arcgis.com/en/pro-app/2.8/help/mapping/properties/plate-carree.htm). This pixelization is identified as `CAR` in [WCS](https://www.atnf.csiro.au/people/mcalabre/WCS/) coordinates. Invented in A.D. 100, CAR maps continue to be a useful and efficient way to map the Earth and Universe. 

```@raw html
<img style="padding: 1em" class=center src="assets/Latitude_and_longitude_graticule_on_a_sphere.svg" width="25%"/>
```
**Figure 1**. *A CAR pixel would be one of these small rectangular patches between lines of fixed longitude and latitude. Note the singularities at the poles. CAR pixels are not the same size across the sky.*

One of the most useful aspects of CAR maps is that the representation of the map in computer memory is a simple 2D array. The indices of these pixels are related linearly to the angles on the sky. This package can plot CAR maps as this array, but it's important to note that showing the array in this way distorts both angles and areas on the sphere, since the sphere has curvature.

```@example intro
using Pixell, Plots # hide
shape, wcs = fullsky_geometry(10.0 * Pixell.degree)  # hide
m = Enmap(randn(shape), wcs)  # hide
plot(m)   # hide
```
**Figure 2**. *Each square is 10 degrees on this map of the full sky, now flattened to array form. This visualization distorts the map considerably, particularly at the poles.*

This flattened projection is a poor visualization of the sky for large sky areas, but is very good in the limit of small sky areas centered on the equator. This coincides with regions where the so-called "flat-sky approximation" holds, and one can perform fast fourier transforms instead of spherical harmonic transforms. Let's look at a smaller, high-resolution map centered on the equator:

```@example intro
using Pixell, Plots # hide
import Pixell: degree, arcminute  # hide
box = [4   -4;           # hide
       -2     2] * degree  # hide
shape, wcs = geometry(Pixell.WCS.WCSTransform, box, 10 * arcminute)  # hide
m = Enmap(randn(shape), wcs)  # hide
plot(m)  # hide
```
**Figure 3**. *If you stay within a few degrees of the equator, the actual sky resembles the CAR array representation.*

## Sky to Pixels and Back
The main way you interact with a pixelization is through the mapping of **angles on the sky** to **pixel indices**, as well as the inverse. In Pixell, these operations are called [`sky2pix`](@ref) and [`pix2sky`](@ref). For a person lying on the ground at the equator with their head pointed towards the north celestial pole, the **right ascension** (RA or ``\alpha``) is the angle pointing leftwards on the sky, and the **declination** (DEC or ``\delta``) is the angle away from the equator pointing upwards.


```@raw html
<img style="padding: 1em" class=center src="assets/Ra_and_dec_on_celestial_sphere.png" width="60%"/>
```

For a CAR map `my_map` and pixel `my_map[i,j]` with row index ``i`` and column index ``j``, the resulting right ascension ``\alpha`` and declination ``\delta`` are given by simple linear expressions

```math
\begin{aligned}
    \alpha &= \alpha_0 + \Delta \alpha \, (i - i_0)   \\
    \delta &= \delta_0 + \Delta \delta \, (j - j_0)  
\end{aligned}
```

A CAR map is defined by these constants: a specific reference pixel with index ``i_0``, ``j_0`` that has reference angle ``\alpha_0``, ``\delta_0``, and global pixel angular size specified by ``\Delta \alpha``, ``\Delta \delta``. These are stored in the FITS file header in a system called the World Coordinate System (WCS). WCS stores these numbers under different FITS header keys,

```math
\begin{aligned}
    \texttt{CRVAL} &= (\alpha_0, \delta_0)   \\
    \texttt{CDELT} &= (\Delta \alpha, \Delta \delta)   \\
    \texttt{CRPIX} &= (i_0, j_0).
\end{aligned}
```

The 
