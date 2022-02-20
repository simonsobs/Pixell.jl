# %%
from pixell import enmap, utils, curvedsky
import numpy as np

# add some simple analytic stuff to the map
def gen_simple(imap, powerlaw_exponent):
    for i in range(imap.shape[0]):
        for j in range(imap.shape[1]):
            imap[i,j] = (i * imap.shape[1] + j + 1)**powerlaw_exponent

shape, wcs = enmap.fullsky_geometry(10.0*utils.degree)
imap = enmap.zeros(shape, wcs=wcs)
gen_simple(imap, 2)

alms = curvedsky.map2alm(imap, lmax=18, method="auto")
np.savetxt("data/simple_analytic_sht.txt", np.column_stack([alms.real, alms.imag]), fmt='%.60g')

alms = curvedsky.map2alm(imap[4:-3, 5:-2], lmax=18, method="auto")
np.savetxt("data/simple_analytic_sht_sliced.txt", np.column_stack([alms.real, alms.imag]), fmt='%.60g')

box = np.array([[-5,10],[5,-10]]) * utils.degree
shape, wcs = enmap.geometry(pos=box,res=60 * utils.arcmin,proj='car')
imap = enmap.zeros(shape, wcs=wcs)
gen_simple(imap, 2.5)
alms = curvedsky.map2alm(imap, lmax=100, method="auto")
np.savetxt("data/simple_box_analytic_sht.txt", np.column_stack([alms.real, alms.imag]), fmt='%.60g')
