"""
Functions to calculate readout of competition experiments

Used in the preparation of the manuscript "Revisiting labelled ligand affinity
in competition experiments" by Shave et.al.
"""


from mpmath import mpf, sqrt, power, mp, fabs
mp.dps = 500  # Set mpmath to use high accuracy


def calc_amount_p(fraction_bound, l, kdax):
    """ Calculate amount of protein for a given fraction bound and KD"""
    return (-(kdax*fraction_bound) - l*fraction_bound + l*fraction_bound*fraction_bound)/(-1 + fraction_bound)


def calc_kdpi_for_fractionl_bound(p, l, i, kdpl, targetflb):
    p=mpf(p)
    l=mpf(l)
    i=mpf(i)
    kdpl=mpf(kdpl)
    kdpi=mpf(kdpi)
    targetflb=mpf(targetflb)
    return float((kdpl*targetflb*(i - p - i*targetflb + kdpl*targetflb + l*targetflb + p*targetflb - l*power(targetflb, 2)))/((-1 + targetflb)*(-p + kdpl*targetflb + l*targetflb + p*targetflb - l*power(targetflb, 2))).real)

def calc_i_for_fractionl_bound(p,l,kdpl,kdpi,targetflb):
    p=mpf(p)
    l=mpf(l)
    kdpl=mpf(kdpl)
    kdpi=mpf(kdpi)
    targetflb=mpf(targetflb)
    return ((-kdpi + kdpi*targetflb - kdpl*targetflb)*(p - kdpl*targetflb - l*targetflb - p*targetflb + l*power(targetflb,2)))/(kdpl*(-targetflb + power(targetflb,2)))


# 1:1:1 competition - see https://stevenshave.github.io/pybindingcurve/simulate_competition.html
# Readout is PL
def competition_pl(p, l, i, kdpl, kdpi):
    """Calculate PL concentration in competition experiment

    Calculate the protein-ligand complex formed in a competition experiment.
    See https://stevenshave.github.io/pybindingcurve/simulate_competition.html
    The correct solution is chosen based on which KD is larger etc.  A
    correction is applied if the KDs are equal.
    """
    p = mpf(p)
    l = mpf(l)
    i = mpf(i)
    kdpl = mpf(kdpl)
    kdpi = mpf(kdpi)
    if fabs(kdpl - kdpi) < 1e-6:
        kdpi += 1e-6
    if kdpl < kdpi:
        return float((
            -(
                p * kdpi
                - p * kdpl
                + i * kdpl
                + kdpi * kdpl
                - power(kdpl, 2)
                + 2 * kdpi * l
                - kdpl * l
            )
            / (3.0 * (-kdpi + kdpl))
            - (
                power(2, 0.3333333333333333)
                * (
                    -power(
                        p * kdpi
                        - p * kdpl
                        + i * kdpl
                        + kdpi * kdpl
                        - power(kdpl, 2)
                        + 2 * kdpi * l
                        - kdpl * l,
                        2,
                    )
                    + 3
                    * (-kdpi + kdpl)
                    * (
                        -2 * p * kdpi * l
                        + p * kdpl * l
                        - i * kdpl * l
                        - kdpi * kdpl * l
                        - kdpi * power(l, 2)
                    )
                )
            )
            / (
                3.0
                * (-kdpi + kdpl)
                * power(
                    -2 * power(p, 3) * power(kdpi, 3)
                    + 6 * power(p, 3) * power(kdpi, 2) * kdpl
                    - 6 * power(p, 2) * i * power(kdpi, 2) * kdpl
                    - 6 * power(p, 2) * power(kdpi, 3) * kdpl
                    - 6 * power(p, 3) * kdpi * power(kdpl, 2)
                    + 12 * power(p, 2) * i * kdpi * power(kdpl, 2)
                    - 6 * p * power(i, 2) * kdpi * power(kdpl, 2)
                    + 18 * power(p, 2) * power(kdpi, 2) * power(kdpl, 2)
                    - 12 * p * i * power(kdpi, 2) * power(kdpl, 2)
                    - 6 * p * power(kdpi, 3) * power(kdpl, 2)
                    + 2 * power(p, 3) * power(kdpl, 3)
                    - 6 * power(p, 2) * i * power(kdpl, 3)
                    + 6 * p * power(i, 2) * power(kdpl, 3)
                    - 2 * power(i, 3) * power(kdpl, 3)
                    - 18 * power(p, 2) * kdpi * power(kdpl, 3)
                    + 24 * p * i * kdpi * power(kdpl, 3)
                    - 6 * power(i, 2) * kdpi * power(kdpl, 3)
                    + 18 * p * power(kdpi, 2) * power(kdpl, 3)
                    - 6 * i * power(kdpi, 2) * power(kdpl, 3)
                    - 2 * power(kdpi, 3) * power(kdpl, 3)
                    + 6 * power(p, 2) * power(kdpl, 4)
                    - 12 * p * i * power(kdpl, 4)
                    + 6 * power(i, 2) * power(kdpl, 4)
                    - 18 * p * kdpi * power(kdpl, 4)
                    + 12 * i * kdpi * power(kdpl, 4)
                    + 6 * power(kdpi, 2) * power(kdpl, 4)
                    + 6 * p * power(kdpl, 5)
                    - 6 * i * power(kdpl, 5)
                    - 6 * kdpi * power(kdpl, 5)
                    + 2 * power(kdpl, 6)
                    + 6 * power(p, 2) * power(kdpi, 3) * l
                    - 15 * power(p, 2) * power(kdpi, 2) * kdpl * l
                    + 3 * p * i * power(kdpi, 2) * kdpl * l
                    + 3 * p * power(kdpi, 3) * kdpl * l
                    + 12 * power(p, 2) * kdpi * power(kdpl, 2) * l
                    - 9 * p * i * kdpi * power(kdpl, 2) * l
                    - 3 * power(i, 2) * kdpi * power(kdpl, 2) * l
                    - 3 * p * power(kdpi, 2) * power(kdpl, 2) * l
                    - 6 * i * power(kdpi, 2) * power(kdpl, 2) * l
                    - 3 * power(kdpi, 3) * power(kdpl, 2) * l
                    - 3 * power(p, 2) * power(kdpl, 3) * l
                    + 6 * p * i * power(kdpl, 3) * l
                    - 3 * power(i, 2) * power(kdpl, 3) * l
                    - 3 * p * kdpi * power(kdpl, 3) * l
                    + 9 * i * kdpi * power(kdpl, 3) * l
                    + 12 * power(kdpi, 2) * power(kdpl, 3) * l
                    + 3 * p * power(kdpl, 4) * l
                    - 3 * i * power(kdpl, 4) * l
                    - 15 * kdpi * power(kdpl, 4) * l
                    + 6 * power(kdpl, 5) * l
                    - 6 * p * power(kdpi, 3) * power(l, 2)
                    + 12 * p * power(kdpi, 2) * kdpl * power(l, 2)
                    + 3 * i * power(kdpi, 2) * kdpl * power(l, 2)
                    + 3 * power(kdpi, 3) * kdpl * power(l, 2)
                    - 3 * p * kdpi * power(kdpl, 2) * power(l, 2)
                    - 12 * i * kdpi * power(kdpl, 2) * power(l, 2)
                    + 3 * power(kdpi, 2) * power(kdpl, 2) * power(l, 2)
                    - 3 * p * power(kdpl, 3) * power(l, 2)
                    + 3 * i * power(kdpl, 3) * power(l, 2)
                    - 12 * kdpi * power(kdpl, 3) * power(l, 2)
                    + 6 * power(kdpl, 4) * power(l, 2)
                    + 2 * power(kdpi, 3) * power(l, 3)
                    - 3 * power(kdpi, 2) * kdpl * power(l, 3)
                    - 3 * kdpi * power(kdpl, 2) * power(l, 3)
                    + 2 * power(kdpl, 3) * power(l, 3)
                    + sqrt(
                        power(
                            -2 * power(p, 3) * power(kdpi, 3)
                            + 6 * power(p, 3) * power(kdpi, 2) * kdpl
                            - 6 * power(p, 2) * i * power(kdpi, 2) * kdpl
                            - 6 * power(p, 2) * power(kdpi, 3) * kdpl
                            - 6 * power(p, 3) * kdpi * power(kdpl, 2)
                            + 12 * power(p, 2) * i * kdpi * power(kdpl, 2)
                            - 6 * p * power(i, 2) * kdpi * power(kdpl, 2)
                            + 18 * power(p, 2) * power(kdpi, 2) *
                            power(kdpl, 2)
                            - 12 * p * i * power(kdpi, 2) * power(kdpl, 2)
                            - 6 * p * power(kdpi, 3) * power(kdpl, 2)
                            + 2 * power(p, 3) * power(kdpl, 3)
                            - 6 * power(p, 2) * i * power(kdpl, 3)
                            + 6 * p * power(i, 2) * power(kdpl, 3)
                            - 2 * power(i, 3) * power(kdpl, 3)
                            - 18 * power(p, 2) * kdpi * power(kdpl, 3)
                            + 24 * p * i * kdpi * power(kdpl, 3)
                            - 6 * power(i, 2) * kdpi * power(kdpl, 3)
                            + 18 * p * power(kdpi, 2) * power(kdpl, 3)
                            - 6 * i * power(kdpi, 2) * power(kdpl, 3)
                            - 2 * power(kdpi, 3) * power(kdpl, 3)
                            + 6 * power(p, 2) * power(kdpl, 4)
                            - 12 * p * i * power(kdpl, 4)
                            + 6 * power(i, 2) * power(kdpl, 4)
                            - 18 * p * kdpi * power(kdpl, 4)
                            + 12 * i * kdpi * power(kdpl, 4)
                            + 6 * power(kdpi, 2) * power(kdpl, 4)
                            + 6 * p * power(kdpl, 5)
                            - 6 * i * power(kdpl, 5)
                            - 6 * kdpi * power(kdpl, 5)
                            + 2 * power(kdpl, 6)
                            + 6 * power(p, 2) * power(kdpi, 3) * l
                            - 15 * power(p, 2) * power(kdpi, 2) * kdpl * l
                            + 3 * p * i * power(kdpi, 2) * kdpl * l
                            + 3 * p * power(kdpi, 3) * kdpl * l
                            + 12 * power(p, 2) * kdpi * power(kdpl, 2) * l
                            - 9 * p * i * kdpi * power(kdpl, 2) * l
                            - 3 * power(i, 2) * kdpi * power(kdpl, 2) * l
                            - 3 * p * power(kdpi, 2) * power(kdpl, 2) * l
                            - 6 * i * power(kdpi, 2) * power(kdpl, 2) * l
                            - 3 * power(kdpi, 3) * power(kdpl, 2) * l
                            - 3 * power(p, 2) * power(kdpl, 3) * l
                            + 6 * p * i * power(kdpl, 3) * l
                            - 3 * power(i, 2) * power(kdpl, 3) * l
                            - 3 * p * kdpi * power(kdpl, 3) * l
                            + 9 * i * kdpi * power(kdpl, 3) * l
                            + 12 * power(kdpi, 2) * power(kdpl, 3) * l
                            + 3 * p * power(kdpl, 4) * l
                            - 3 * i * power(kdpl, 4) * l
                            - 15 * kdpi * power(kdpl, 4) * l
                            + 6 * power(kdpl, 5) * l
                            - 6 * p * power(kdpi, 3) * power(l, 2)
                            + 12 * p * power(kdpi, 2) * kdpl * power(l, 2)
                            + 3 * i * power(kdpi, 2) * kdpl * power(l, 2)
                            + 3 * power(kdpi, 3) * kdpl * power(l, 2)
                            - 3 * p * kdpi * power(kdpl, 2) * power(l, 2)
                            - 12 * i * kdpi * power(kdpl, 2) * power(l, 2)
                            + 3 * power(kdpi, 2) * power(kdpl, 2) * power(l, 2)
                            - 3 * p * power(kdpl, 3) * power(l, 2)
                            + 3 * i * power(kdpl, 3) * power(l, 2)
                            - 12 * kdpi * power(kdpl, 3) * power(l, 2)
                            + 6 * power(kdpl, 4) * power(l, 2)
                            + 2 * power(kdpi, 3) * power(l, 3)
                            - 3 * power(kdpi, 2) * kdpl * power(l, 3)
                            - 3 * kdpi * power(kdpl, 2) * power(l, 3)
                            + 2 * power(kdpl, 3) * power(l, 3),
                            2,
                        )
                        + 4
                        * power(
                            -power(
                                p * kdpi
                                - p * kdpl
                                + i * kdpl
                                + kdpi * kdpl
                                - power(kdpl, 2)
                                + 2 * kdpi * l
                                - kdpl * l,
                                2,
                            )
                            + 3
                            * (-kdpi + kdpl)
                            * (
                                -2 * p * kdpi * l
                                + p * kdpl * l
                                - i * kdpl * l
                                - kdpi * kdpl * l
                                - kdpi * power(l, 2)
                            ),
                            3,
                        )
                    ),
                    0.3333333333333333,
                )
            )
            + power(
                -2 * power(p, 3) * power(kdpi, 3)
                + 6 * power(p, 3) * power(kdpi, 2) * kdpl
                - 6 * power(p, 2) * i * power(kdpi, 2) * kdpl
                - 6 * power(p, 2) * power(kdpi, 3) * kdpl
                - 6 * power(p, 3) * kdpi * power(kdpl, 2)
                + 12 * power(p, 2) * i * kdpi * power(kdpl, 2)
                - 6 * p * power(i, 2) * kdpi * power(kdpl, 2)
                + 18 * power(p, 2) * power(kdpi, 2) * power(kdpl, 2)
                - 12 * p * i * power(kdpi, 2) * power(kdpl, 2)
                - 6 * p * power(kdpi, 3) * power(kdpl, 2)
                + 2 * power(p, 3) * power(kdpl, 3)
                - 6 * power(p, 2) * i * power(kdpl, 3)
                + 6 * p * power(i, 2) * power(kdpl, 3)
                - 2 * power(i, 3) * power(kdpl, 3)
                - 18 * power(p, 2) * kdpi * power(kdpl, 3)
                + 24 * p * i * kdpi * power(kdpl, 3)
                - 6 * power(i, 2) * kdpi * power(kdpl, 3)
                + 18 * p * power(kdpi, 2) * power(kdpl, 3)
                - 6 * i * power(kdpi, 2) * power(kdpl, 3)
                - 2 * power(kdpi, 3) * power(kdpl, 3)
                + 6 * power(p, 2) * power(kdpl, 4)
                - 12 * p * i * power(kdpl, 4)
                + 6 * power(i, 2) * power(kdpl, 4)
                - 18 * p * kdpi * power(kdpl, 4)
                + 12 * i * kdpi * power(kdpl, 4)
                + 6 * power(kdpi, 2) * power(kdpl, 4)
                + 6 * p * power(kdpl, 5)
                - 6 * i * power(kdpl, 5)
                - 6 * kdpi * power(kdpl, 5)
                + 2 * power(kdpl, 6)
                + 6 * power(p, 2) * power(kdpi, 3) * l
                - 15 * power(p, 2) * power(kdpi, 2) * kdpl * l
                + 3 * p * i * power(kdpi, 2) * kdpl * l
                + 3 * p * power(kdpi, 3) * kdpl * l
                + 12 * power(p, 2) * kdpi * power(kdpl, 2) * l
                - 9 * p * i * kdpi * power(kdpl, 2) * l
                - 3 * power(i, 2) * kdpi * power(kdpl, 2) * l
                - 3 * p * power(kdpi, 2) * power(kdpl, 2) * l
                - 6 * i * power(kdpi, 2) * power(kdpl, 2) * l
                - 3 * power(kdpi, 3) * power(kdpl, 2) * l
                - 3 * power(p, 2) * power(kdpl, 3) * l
                + 6 * p * i * power(kdpl, 3) * l
                - 3 * power(i, 2) * power(kdpl, 3) * l
                - 3 * p * kdpi * power(kdpl, 3) * l
                + 9 * i * kdpi * power(kdpl, 3) * l
                + 12 * power(kdpi, 2) * power(kdpl, 3) * l
                + 3 * p * power(kdpl, 4) * l
                - 3 * i * power(kdpl, 4) * l
                - 15 * kdpi * power(kdpl, 4) * l
                + 6 * power(kdpl, 5) * l
                - 6 * p * power(kdpi, 3) * power(l, 2)
                + 12 * p * power(kdpi, 2) * kdpl * power(l, 2)
                + 3 * i * power(kdpi, 2) * kdpl * power(l, 2)
                + 3 * power(kdpi, 3) * kdpl * power(l, 2)
                - 3 * p * kdpi * power(kdpl, 2) * power(l, 2)
                - 12 * i * kdpi * power(kdpl, 2) * power(l, 2)
                + 3 * power(kdpi, 2) * power(kdpl, 2) * power(l, 2)
                - 3 * p * power(kdpl, 3) * power(l, 2)
                + 3 * i * power(kdpl, 3) * power(l, 2)
                - 12 * kdpi * power(kdpl, 3) * power(l, 2)
                + 6 * power(kdpl, 4) * power(l, 2)
                + 2 * power(kdpi, 3) * power(l, 3)
                - 3 * power(kdpi, 2) * kdpl * power(l, 3)
                - 3 * kdpi * power(kdpl, 2) * power(l, 3)
                + 2 * power(kdpl, 3) * power(l, 3)
                + sqrt(
                    power(
                        -2 * power(p, 3) * power(kdpi, 3)
                        + 6 * power(p, 3) * power(kdpi, 2) * kdpl
                        - 6 * power(p, 2) * i * power(kdpi, 2) * kdpl
                        - 6 * power(p, 2) * power(kdpi, 3) * kdpl
                        - 6 * power(p, 3) * kdpi * power(kdpl, 2)
                        + 12 * power(p, 2) * i * kdpi * power(kdpl, 2)
                        - 6 * p * power(i, 2) * kdpi * power(kdpl, 2)
                        + 18 * power(p, 2) * power(kdpi, 2) * power(kdpl, 2)
                        - 12 * p * i * power(kdpi, 2) * power(kdpl, 2)
                        - 6 * p * power(kdpi, 3) * power(kdpl, 2)
                        + 2 * power(p, 3) * power(kdpl, 3)
                        - 6 * power(p, 2) * i * power(kdpl, 3)
                        + 6 * p * power(i, 2) * power(kdpl, 3)
                        - 2 * power(i, 3) * power(kdpl, 3)
                        - 18 * power(p, 2) * kdpi * power(kdpl, 3)
                        + 24 * p * i * kdpi * power(kdpl, 3)
                        - 6 * power(i, 2) * kdpi * power(kdpl, 3)
                        + 18 * p * power(kdpi, 2) * power(kdpl, 3)
                        - 6 * i * power(kdpi, 2) * power(kdpl, 3)
                        - 2 * power(kdpi, 3) * power(kdpl, 3)
                        + 6 * power(p, 2) * power(kdpl, 4)
                        - 12 * p * i * power(kdpl, 4)
                        + 6 * power(i, 2) * power(kdpl, 4)
                        - 18 * p * kdpi * power(kdpl, 4)
                        + 12 * i * kdpi * power(kdpl, 4)
                        + 6 * power(kdpi, 2) * power(kdpl, 4)
                        + 6 * p * power(kdpl, 5)
                        - 6 * i * power(kdpl, 5)
                        - 6 * kdpi * power(kdpl, 5)
                        + 2 * power(kdpl, 6)
                        + 6 * power(p, 2) * power(kdpi, 3) * l
                        - 15 * power(p, 2) * power(kdpi, 2) * kdpl * l
                        + 3 * p * i * power(kdpi, 2) * kdpl * l
                        + 3 * p * power(kdpi, 3) * kdpl * l
                        + 12 * power(p, 2) * kdpi * power(kdpl, 2) * l
                        - 9 * p * i * kdpi * power(kdpl, 2) * l
                        - 3 * power(i, 2) * kdpi * power(kdpl, 2) * l
                        - 3 * p * power(kdpi, 2) * power(kdpl, 2) * l
                        - 6 * i * power(kdpi, 2) * power(kdpl, 2) * l
                        - 3 * power(kdpi, 3) * power(kdpl, 2) * l
                        - 3 * power(p, 2) * power(kdpl, 3) * l
                        + 6 * p * i * power(kdpl, 3) * l
                        - 3 * power(i, 2) * power(kdpl, 3) * l
                        - 3 * p * kdpi * power(kdpl, 3) * l
                        + 9 * i * kdpi * power(kdpl, 3) * l
                        + 12 * power(kdpi, 2) * power(kdpl, 3) * l
                        + 3 * p * power(kdpl, 4) * l
                        - 3 * i * power(kdpl, 4) * l
                        - 15 * kdpi * power(kdpl, 4) * l
                        + 6 * power(kdpl, 5) * l
                        - 6 * p * power(kdpi, 3) * power(l, 2)
                        + 12 * p * power(kdpi, 2) * kdpl * power(l, 2)
                        + 3 * i * power(kdpi, 2) * kdpl * power(l, 2)
                        + 3 * power(kdpi, 3) * kdpl * power(l, 2)
                        - 3 * p * kdpi * power(kdpl, 2) * power(l, 2)
                        - 12 * i * kdpi * power(kdpl, 2) * power(l, 2)
                        + 3 * power(kdpi, 2) * power(kdpl, 2) * power(l, 2)
                        - 3 * p * power(kdpl, 3) * power(l, 2)
                        + 3 * i * power(kdpl, 3) * power(l, 2)
                        - 12 * kdpi * power(kdpl, 3) * power(l, 2)
                        + 6 * power(kdpl, 4) * power(l, 2)
                        + 2 * power(kdpi, 3) * power(l, 3)
                        - 3 * power(kdpi, 2) * kdpl * power(l, 3)
                        - 3 * kdpi * power(kdpl, 2) * power(l, 3)
                        + 2 * power(kdpl, 3) * power(l, 3),
                        2,
                    )
                    + 4
                    * power(
                        -power(
                            p * kdpi
                            - p * kdpl
                            + i * kdpl
                            + kdpi * kdpl
                            - power(kdpl, 2)
                            + 2 * kdpi * l
                            - kdpl * l,
                            2,
                        )
                        + 3
                        * (-kdpi + kdpl)
                        * (
                            -2 * p * kdpi * l
                            + p * kdpl * l
                            - i * kdpl * l
                            - kdpi * kdpl * l
                            - kdpi * power(l, 2)
                        ),
                        3,
                    )
                ),
                0.3333333333333333,
            )
            / (3.0 * power(2, 0.3333333333333333) * (-kdpi + kdpl))
        ).real)
    else:
        return float((
            -(
                p * kdpi
                - p * kdpl
                + i * kdpl
                + kdpi * kdpl
                - power(kdpl, 2)
                + 2 * kdpi * l
                - kdpl * l
            )
            / (3.0 * (-kdpi + kdpl))
            + (
                (1 - complex(0, 1) * sqrt(3))
                * (
                    -power(
                        p * kdpi
                        - p * kdpl
                        + i * kdpl
                        + kdpi * kdpl
                        - power(kdpl, 2)
                        + 2 * kdpi * l
                        - kdpl * l,
                        2,
                    )
                    + 3
                    * (-kdpi + kdpl)
                    * (
                        -2 * p * kdpi * l
                        + p * kdpl * l
                        - i * kdpl * l
                        - kdpi * kdpl * l
                        - kdpi * power(l, 2)
                    )
                )
            )
            / (
                3.0
                * power(2, 0.6666666666666666)
                * (-kdpi + kdpl)
                * power(
                    -2 * power(p, 3) * power(kdpi, 3)
                    + 6 * power(p, 3) * power(kdpi, 2) * kdpl
                    - 6 * power(p, 2) * i * power(kdpi, 2) * kdpl
                    - 6 * power(p, 2) * power(kdpi, 3) * kdpl
                    - 6 * power(p, 3) * kdpi * power(kdpl, 2)
                    + 12 * power(p, 2) * i * kdpi * power(kdpl, 2)
                    - 6 * p * power(i, 2) * kdpi * power(kdpl, 2)
                    + 18 * power(p, 2) * power(kdpi, 2) * power(kdpl, 2)
                    - 12 * p * i * power(kdpi, 2) * power(kdpl, 2)
                    - 6 * p * power(kdpi, 3) * power(kdpl, 2)
                    + 2 * power(p, 3) * power(kdpl, 3)
                    - 6 * power(p, 2) * i * power(kdpl, 3)
                    + 6 * p * power(i, 2) * power(kdpl, 3)
                    - 2 * power(i, 3) * power(kdpl, 3)
                    - 18 * power(p, 2) * kdpi * power(kdpl, 3)
                    + 24 * p * i * kdpi * power(kdpl, 3)
                    - 6 * power(i, 2) * kdpi * power(kdpl, 3)
                    + 18 * p * power(kdpi, 2) * power(kdpl, 3)
                    - 6 * i * power(kdpi, 2) * power(kdpl, 3)
                    - 2 * power(kdpi, 3) * power(kdpl, 3)
                    + 6 * power(p, 2) * power(kdpl, 4)
                    - 12 * p * i * power(kdpl, 4)
                    + 6 * power(i, 2) * power(kdpl, 4)
                    - 18 * p * kdpi * power(kdpl, 4)
                    + 12 * i * kdpi * power(kdpl, 4)
                    + 6 * power(kdpi, 2) * power(kdpl, 4)
                    + 6 * p * power(kdpl, 5)
                    - 6 * i * power(kdpl, 5)
                    - 6 * kdpi * power(kdpl, 5)
                    + 2 * power(kdpl, 6)
                    + 6 * power(p, 2) * power(kdpi, 3) * l
                    - 15 * power(p, 2) * power(kdpi, 2) * kdpl * l
                    + 3 * p * i * power(kdpi, 2) * kdpl * l
                    + 3 * p * power(kdpi, 3) * kdpl * l
                    + 12 * power(p, 2) * kdpi * power(kdpl, 2) * l
                    - 9 * p * i * kdpi * power(kdpl, 2) * l
                    - 3 * power(i, 2) * kdpi * power(kdpl, 2) * l
                    - 3 * p * power(kdpi, 2) * power(kdpl, 2) * l
                    - 6 * i * power(kdpi, 2) * power(kdpl, 2) * l
                    - 3 * power(kdpi, 3) * power(kdpl, 2) * l
                    - 3 * power(p, 2) * power(kdpl, 3) * l
                    + 6 * p * i * power(kdpl, 3) * l
                    - 3 * power(i, 2) * power(kdpl, 3) * l
                    - 3 * p * kdpi * power(kdpl, 3) * l
                    + 9 * i * kdpi * power(kdpl, 3) * l
                    + 12 * power(kdpi, 2) * power(kdpl, 3) * l
                    + 3 * p * power(kdpl, 4) * l
                    - 3 * i * power(kdpl, 4) * l
                    - 15 * kdpi * power(kdpl, 4) * l
                    + 6 * power(kdpl, 5) * l
                    - 6 * p * power(kdpi, 3) * power(l, 2)
                    + 12 * p * power(kdpi, 2) * kdpl * power(l, 2)
                    + 3 * i * power(kdpi, 2) * kdpl * power(l, 2)
                    + 3 * power(kdpi, 3) * kdpl * power(l, 2)
                    - 3 * p * kdpi * power(kdpl, 2) * power(l, 2)
                    - 12 * i * kdpi * power(kdpl, 2) * power(l, 2)
                    + 3 * power(kdpi, 2) * power(kdpl, 2) * power(l, 2)
                    - 3 * p * power(kdpl, 3) * power(l, 2)
                    + 3 * i * power(kdpl, 3) * power(l, 2)
                    - 12 * kdpi * power(kdpl, 3) * power(l, 2)
                    + 6 * power(kdpl, 4) * power(l, 2)
                    + 2 * power(kdpi, 3) * power(l, 3)
                    - 3 * power(kdpi, 2) * kdpl * power(l, 3)
                    - 3 * kdpi * power(kdpl, 2) * power(l, 3)
                    + 2 * power(kdpl, 3) * power(l, 3)
                    + sqrt(
                        power(
                            -2 * power(p, 3) * power(kdpi, 3)
                            + 6 * power(p, 3) * power(kdpi, 2) * kdpl
                            - 6 * power(p, 2) * i * power(kdpi, 2) * kdpl
                            - 6 * power(p, 2) * power(kdpi, 3) * kdpl
                            - 6 * power(p, 3) * kdpi * power(kdpl, 2)
                            + 12 * power(p, 2) * i * kdpi * power(kdpl, 2)
                            - 6 * p * power(i, 2) * kdpi * power(kdpl, 2)
                            + 18 * power(p, 2) * power(kdpi, 2) *
                            power(kdpl, 2)
                            - 12 * p * i * power(kdpi, 2) * power(kdpl, 2)
                            - 6 * p * power(kdpi, 3) * power(kdpl, 2)
                            + 2 * power(p, 3) * power(kdpl, 3)
                            - 6 * power(p, 2) * i * power(kdpl, 3)
                            + 6 * p * power(i, 2) * power(kdpl, 3)
                            - 2 * power(i, 3) * power(kdpl, 3)
                            - 18 * power(p, 2) * kdpi * power(kdpl, 3)
                            + 24 * p * i * kdpi * power(kdpl, 3)
                            - 6 * power(i, 2) * kdpi * power(kdpl, 3)
                            + 18 * p * power(kdpi, 2) * power(kdpl, 3)
                            - 6 * i * power(kdpi, 2) * power(kdpl, 3)
                            - 2 * power(kdpi, 3) * power(kdpl, 3)
                            + 6 * power(p, 2) * power(kdpl, 4)
                            - 12 * p * i * power(kdpl, 4)
                            + 6 * power(i, 2) * power(kdpl, 4)
                            - 18 * p * kdpi * power(kdpl, 4)
                            + 12 * i * kdpi * power(kdpl, 4)
                            + 6 * power(kdpi, 2) * power(kdpl, 4)
                            + 6 * p * power(kdpl, 5)
                            - 6 * i * power(kdpl, 5)
                            - 6 * kdpi * power(kdpl, 5)
                            + 2 * power(kdpl, 6)
                            + 6 * power(p, 2) * power(kdpi, 3) * l
                            - 15 * power(p, 2) * power(kdpi, 2) * kdpl * l
                            + 3 * p * i * power(kdpi, 2) * kdpl * l
                            + 3 * p * power(kdpi, 3) * kdpl * l
                            + 12 * power(p, 2) * kdpi * power(kdpl, 2) * l
                            - 9 * p * i * kdpi * power(kdpl, 2) * l
                            - 3 * power(i, 2) * kdpi * power(kdpl, 2) * l
                            - 3 * p * power(kdpi, 2) * power(kdpl, 2) * l
                            - 6 * i * power(kdpi, 2) * power(kdpl, 2) * l
                            - 3 * power(kdpi, 3) * power(kdpl, 2) * l
                            - 3 * power(p, 2) * power(kdpl, 3) * l
                            + 6 * p * i * power(kdpl, 3) * l
                            - 3 * power(i, 2) * power(kdpl, 3) * l
                            - 3 * p * kdpi * power(kdpl, 3) * l
                            + 9 * i * kdpi * power(kdpl, 3) * l
                            + 12 * power(kdpi, 2) * power(kdpl, 3) * l
                            + 3 * p * power(kdpl, 4) * l
                            - 3 * i * power(kdpl, 4) * l
                            - 15 * kdpi * power(kdpl, 4) * l
                            + 6 * power(kdpl, 5) * l
                            - 6 * p * power(kdpi, 3) * power(l, 2)
                            + 12 * p * power(kdpi, 2) * kdpl * power(l, 2)
                            + 3 * i * power(kdpi, 2) * kdpl * power(l, 2)
                            + 3 * power(kdpi, 3) * kdpl * power(l, 2)
                            - 3 * p * kdpi * power(kdpl, 2) * power(l, 2)
                            - 12 * i * kdpi * power(kdpl, 2) * power(l, 2)
                            + 3 * power(kdpi, 2) * power(kdpl, 2) * power(l, 2)
                            - 3 * p * power(kdpl, 3) * power(l, 2)
                            + 3 * i * power(kdpl, 3) * power(l, 2)
                            - 12 * kdpi * power(kdpl, 3) * power(l, 2)
                            + 6 * power(kdpl, 4) * power(l, 2)
                            + 2 * power(kdpi, 3) * power(l, 3)
                            - 3 * power(kdpi, 2) * kdpl * power(l, 3)
                            - 3 * kdpi * power(kdpl, 2) * power(l, 3)
                            + 2 * power(kdpl, 3) * power(l, 3),
                            2,
                        )
                        + 4
                        * power(
                            -power(
                                p * kdpi
                                - p * kdpl
                                + i * kdpl
                                + kdpi * kdpl
                                - power(kdpl, 2)
                                + 2 * kdpi * l
                                - kdpl * l,
                                2,
                            )
                            + 3
                            * (-kdpi + kdpl)
                            * (
                                -2 * p * kdpi * l
                                + p * kdpl * l
                                - i * kdpl * l
                                - kdpi * kdpl * l
                                - kdpi * power(l, 2)
                            ),
                            3,
                        )
                    ),
                    0.3333333333333333,
                )
            )
            - (
                (1 + complex(0, 1) * sqrt(3))
                * power(
                    -2 * power(p, 3) * power(kdpi, 3)
                    + 6 * power(p, 3) * power(kdpi, 2) * kdpl
                    - 6 * power(p, 2) * i * power(kdpi, 2) * kdpl
                    - 6 * power(p, 2) * power(kdpi, 3) * kdpl
                    - 6 * power(p, 3) * kdpi * power(kdpl, 2)
                    + 12 * power(p, 2) * i * kdpi * power(kdpl, 2)
                    - 6 * p * power(i, 2) * kdpi * power(kdpl, 2)
                    + 18 * power(p, 2) * power(kdpi, 2) * power(kdpl, 2)
                    - 12 * p * i * power(kdpi, 2) * power(kdpl, 2)
                    - 6 * p * power(kdpi, 3) * power(kdpl, 2)
                    + 2 * power(p, 3) * power(kdpl, 3)
                    - 6 * power(p, 2) * i * power(kdpl, 3)
                    + 6 * p * power(i, 2) * power(kdpl, 3)
                    - 2 * power(i, 3) * power(kdpl, 3)
                    - 18 * power(p, 2) * kdpi * power(kdpl, 3)
                    + 24 * p * i * kdpi * power(kdpl, 3)
                    - 6 * power(i, 2) * kdpi * power(kdpl, 3)
                    + 18 * p * power(kdpi, 2) * power(kdpl, 3)
                    - 6 * i * power(kdpi, 2) * power(kdpl, 3)
                    - 2 * power(kdpi, 3) * power(kdpl, 3)
                    + 6 * power(p, 2) * power(kdpl, 4)
                    - 12 * p * i * power(kdpl, 4)
                    + 6 * power(i, 2) * power(kdpl, 4)
                    - 18 * p * kdpi * power(kdpl, 4)
                    + 12 * i * kdpi * power(kdpl, 4)
                    + 6 * power(kdpi, 2) * power(kdpl, 4)
                    + 6 * p * power(kdpl, 5)
                    - 6 * i * power(kdpl, 5)
                    - 6 * kdpi * power(kdpl, 5)
                    + 2 * power(kdpl, 6)
                    + 6 * power(p, 2) * power(kdpi, 3) * l
                    - 15 * power(p, 2) * power(kdpi, 2) * kdpl * l
                    + 3 * p * i * power(kdpi, 2) * kdpl * l
                    + 3 * p * power(kdpi, 3) * kdpl * l
                    + 12 * power(p, 2) * kdpi * power(kdpl, 2) * l
                    - 9 * p * i * kdpi * power(kdpl, 2) * l
                    - 3 * power(i, 2) * kdpi * power(kdpl, 2) * l
                    - 3 * p * power(kdpi, 2) * power(kdpl, 2) * l
                    - 6 * i * power(kdpi, 2) * power(kdpl, 2) * l
                    - 3 * power(kdpi, 3) * power(kdpl, 2) * l
                    - 3 * power(p, 2) * power(kdpl, 3) * l
                    + 6 * p * i * power(kdpl, 3) * l
                    - 3 * power(i, 2) * power(kdpl, 3) * l
                    - 3 * p * kdpi * power(kdpl, 3) * l
                    + 9 * i * kdpi * power(kdpl, 3) * l
                    + 12 * power(kdpi, 2) * power(kdpl, 3) * l
                    + 3 * p * power(kdpl, 4) * l
                    - 3 * i * power(kdpl, 4) * l
                    - 15 * kdpi * power(kdpl, 4) * l
                    + 6 * power(kdpl, 5) * l
                    - 6 * p * power(kdpi, 3) * power(l, 2)
                    + 12 * p * power(kdpi, 2) * kdpl * power(l, 2)
                    + 3 * i * power(kdpi, 2) * kdpl * power(l, 2)
                    + 3 * power(kdpi, 3) * kdpl * power(l, 2)
                    - 3 * p * kdpi * power(kdpl, 2) * power(l, 2)
                    - 12 * i * kdpi * power(kdpl, 2) * power(l, 2)
                    + 3 * power(kdpi, 2) * power(kdpl, 2) * power(l, 2)
                    - 3 * p * power(kdpl, 3) * power(l, 2)
                    + 3 * i * power(kdpl, 3) * power(l, 2)
                    - 12 * kdpi * power(kdpl, 3) * power(l, 2)
                    + 6 * power(kdpl, 4) * power(l, 2)
                    + 2 * power(kdpi, 3) * power(l, 3)
                    - 3 * power(kdpi, 2) * kdpl * power(l, 3)
                    - 3 * kdpi * power(kdpl, 2) * power(l, 3)
                    + 2 * power(kdpl, 3) * power(l, 3)
                    + sqrt(
                        power(
                            -2 * power(p, 3) * power(kdpi, 3)
                            + 6 * power(p, 3) * power(kdpi, 2) * kdpl
                            - 6 * power(p, 2) * i * power(kdpi, 2) * kdpl
                            - 6 * power(p, 2) * power(kdpi, 3) * kdpl
                            - 6 * power(p, 3) * kdpi * power(kdpl, 2)
                            + 12 * power(p, 2) * i * kdpi * power(kdpl, 2)
                            - 6 * p * power(i, 2) * kdpi * power(kdpl, 2)
                            + 18 * power(p, 2) * power(kdpi, 2) *
                            power(kdpl, 2)
                            - 12 * p * i * power(kdpi, 2) * power(kdpl, 2)
                            - 6 * p * power(kdpi, 3) * power(kdpl, 2)
                            + 2 * power(p, 3) * power(kdpl, 3)
                            - 6 * power(p, 2) * i * power(kdpl, 3)
                            + 6 * p * power(i, 2) * power(kdpl, 3)
                            - 2 * power(i, 3) * power(kdpl, 3)
                            - 18 * power(p, 2) * kdpi * power(kdpl, 3)
                            + 24 * p * i * kdpi * power(kdpl, 3)
                            - 6 * power(i, 2) * kdpi * power(kdpl, 3)
                            + 18 * p * power(kdpi, 2) * power(kdpl, 3)
                            - 6 * i * power(kdpi, 2) * power(kdpl, 3)
                            - 2 * power(kdpi, 3) * power(kdpl, 3)
                            + 6 * power(p, 2) * power(kdpl, 4)
                            - 12 * p * i * power(kdpl, 4)
                            + 6 * power(i, 2) * power(kdpl, 4)
                            - 18 * p * kdpi * power(kdpl, 4)
                            + 12 * i * kdpi * power(kdpl, 4)
                            + 6 * power(kdpi, 2) * power(kdpl, 4)
                            + 6 * p * power(kdpl, 5)
                            - 6 * i * power(kdpl, 5)
                            - 6 * kdpi * power(kdpl, 5)
                            + 2 * power(kdpl, 6)
                            + 6 * power(p, 2) * power(kdpi, 3) * l
                            - 15 * power(p, 2) * power(kdpi, 2) * kdpl * l
                            + 3 * p * i * power(kdpi, 2) * kdpl * l
                            + 3 * p * power(kdpi, 3) * kdpl * l
                            + 12 * power(p, 2) * kdpi * power(kdpl, 2) * l
                            - 9 * p * i * kdpi * power(kdpl, 2) * l
                            - 3 * power(i, 2) * kdpi * power(kdpl, 2) * l
                            - 3 * p * power(kdpi, 2) * power(kdpl, 2) * l
                            - 6 * i * power(kdpi, 2) * power(kdpl, 2) * l
                            - 3 * power(kdpi, 3) * power(kdpl, 2) * l
                            - 3 * power(p, 2) * power(kdpl, 3) * l
                            + 6 * p * i * power(kdpl, 3) * l
                            - 3 * power(i, 2) * power(kdpl, 3) * l
                            - 3 * p * kdpi * power(kdpl, 3) * l
                            + 9 * i * kdpi * power(kdpl, 3) * l
                            + 12 * power(kdpi, 2) * power(kdpl, 3) * l
                            + 3 * p * power(kdpl, 4) * l
                            - 3 * i * power(kdpl, 4) * l
                            - 15 * kdpi * power(kdpl, 4) * l
                            + 6 * power(kdpl, 5) * l
                            - 6 * p * power(kdpi, 3) * power(l, 2)
                            + 12 * p * power(kdpi, 2) * kdpl * power(l, 2)
                            + 3 * i * power(kdpi, 2) * kdpl * power(l, 2)
                            + 3 * power(kdpi, 3) * kdpl * power(l, 2)
                            - 3 * p * kdpi * power(kdpl, 2) * power(l, 2)
                            - 12 * i * kdpi * power(kdpl, 2) * power(l, 2)
                            + 3 * power(kdpi, 2) * power(kdpl, 2) * power(l, 2)
                            - 3 * p * power(kdpl, 3) * power(l, 2)
                            + 3 * i * power(kdpl, 3) * power(l, 2)
                            - 12 * kdpi * power(kdpl, 3) * power(l, 2)
                            + 6 * power(kdpl, 4) * power(l, 2)
                            + 2 * power(kdpi, 3) * power(l, 3)
                            - 3 * power(kdpi, 2) * kdpl * power(l, 3)
                            - 3 * kdpi * power(kdpl, 2) * power(l, 3)
                            + 2 * power(kdpl, 3) * power(l, 3),
                            2,
                        )
                        + 4
                        * power(
                            -power(
                                p * kdpi
                                - p * kdpl
                                + i * kdpl
                                + kdpi * kdpl
                                - power(kdpl, 2)
                                + 2 * kdpi * l
                                - kdpl * l,
                                2,
                            )
                            + 3
                            * (-kdpi + kdpl)
                            * (
                                -2 * p * kdpi * l
                                + p * kdpl * l
                                - i * kdpl * l
                                - kdpi * kdpl * l
                                - kdpi * power(l, 2)
                            ),
                            3,
                        )
                    ),
                    0.3333333333333333,
                )
            )
            / (6.0 * power(2, 0.3333333333333333) * (-kdpi + kdpl))
        ).real)
