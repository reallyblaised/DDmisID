"""Models commonly used in fitting.

Credit: Matthew Kenzie
"""

__author__ = "Blaise Delaney"
__email__ = "blaise.delaney at cern.ch"

import numpy as np
from scipy.stats import crystalball, norm, expon

# NOTE: the pdfs (and their composition) generally
# require the appropriate normalisation factor, given
# by the CDF of the distribution, evaluated in the range considered in the fit


def cbl(x, mu, sg, a, n, rng=None):
    cb = crystalball(a, n, mu, sg)
    cbn = 1
    if rng is not None:
        cbn = np.diff(cb.cdf(rng))
    return cb.pdf(x) / cbn


def cbl_cdf(x, mu, sg, a, n, rng=None):
    cb = crystalball(a, n, mu, sg)
    cbn = 1
    if rng is not None:
        cbn = np.diff(cb.cdf(rng))
    return cb.cdf(x) / cbn


def cbr(x, mu, sg, a, n, rng):
    cb = crystalball(a, n, rng[1] - mu + rng[0], sg)
    cbn = np.diff(cb.cdf(rng))
    invx = rng[1] - x + rng[0]
    return cb.pdf(invx) / cbn


def cbr_cdf(x, mu, sg, a, n, rng):
    cb = crystalball(a, n, rng[1] - mu + rng[0], sg)
    cbn = np.diff(cb.cdf(rng))
    invx = rng[1] - x + rng[0]
    return 1 - cb.cdf(invx) / cbn


def dcb(x, f, mul, mur, sgl, sgr, al, ar, nl, nr, rng):
    dcbl = cbl(x, mul, sgl, al, nl, rng)
    dcbr = cbr(x, mur, sgr, ar, nr, rng)
    return f * dcbl + (1 - f) * dcbr


def dcb_cdf(x, f, mul, mur, sgl, sgr, al, ar, nl, nr, rng):
    dcbl = cbl_cdf(x, mul, sgl, al, nl, rng)
    dcbr = cbr_cdf(x, mur, sgr, ar, nr, rng)
    return f * dcbl + (1 - f) * dcbr


def dcbwg(x, f1, f2, mug, mul, mur, sgg, sgl, sgr, al, ar, nl, nr, rng):
    dcbl = cbl(x, mul, sgl, al, nl, rng)
    dcbr = cbr(x, mur, sgr, ar, nr, rng)
    gausp = norm(mug, sgg)
    gausn = np.diff(gausp.cdf(rng))
    gaus = gausp.pdf(x) / gausn
    return f1 * gaus + f2 * dcbl + (1 - f1 - f2) * dcbr


def dcbwg_cdf(x, f1, f2, mug, mul, mur, sgg, sgl, sgr, al, ar, nl, nr, rng):
    dcbl = cbl_cdf(x, mul, sgl, al, nl, rng)
    dcbr = cbr_cdf(x, mur, sgr, ar, nr, rng)
    gausp = norm(mug, sgg)
    gausn = np.diff(gausp.cdf(rng))
    gaus = gausp.cdf(x) / gausn
    return f1 * gaus + f2 * dcbl + (1 - f1 - f2) * dcbr


def dcbw2g(
    x, f1, f2, f3, mug1, mug2, mul, mur, sgg1, sgg2, sgl, sgr, al, ar, nl, nr, rng
):
    dcbl = cbl(x, mul, sgl, al, nl, rng)
    dcbr = cbr(x, mur, sgr, ar, nr, rng)
    gausp1 = norm(mug1, sgg1)
    gausn1 = np.diff(gausp1.cdf(rng))
    gausp2 = norm(mug2, sgg2)
    gausn2 = np.diff(gausp2.cdf(rng))
    gaus1 = gausp1.pdf(x) / gausn1
    gaus2 = gausp2.pdf(x) / gausn2
    return f1 * gaus1 + f2 * gaus2 + f3 * dcbl + (1 - f1 - f2 - f3) * dcbr


def dcbw2g_cdf(
    x, f1, f2, f3, mug1, mug2, mul, mur, sgg1, sgg2, sgl, sgr, al, ar, nl, nr, rng
):
    dcbl = cbl_cdf(x, mul, sgl, al, nl, rng)
    dcbr = cbr_cdf(x, mur, sgr, ar, nr, rng)
    gausp1 = norm(mug1, sgg1)
    gausn1 = np.diff(gausp1.cdf(rng))
    gaus1 = gausp1.cdf(x) / gausn1
    gausp2 = norm(mug2, sgg2)
    gausn2 = np.diff(gausp2.cdf(rng))
    gaus2 = gausp2.cdf(x) / gausn2
    return f1 * gaus1 + f2 * gaus2 + f3 * dcbl + (1 - f1 - f2 - f3) * dcbr


def expon_pdf(x, lb, mrange):
    exp = expon(mrange[0], lb)
    return exp.pdf(x) / np.diff(exp.cdf(mrange))


def expon_cdf(x, lb, mrange):
    exp = expon(mrange[0], lb)
    return exp.cdf(x) / np.diff(exp.cdf(mrange))
