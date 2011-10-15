"""Astronomical and physical constants, and unit conversions."""
from __future__ import division
from collections import namedtuple

# Use CODATA 2010 for physics, PDG RPP for astro. constants.  Once I
# have access to Allen's Astrophysical Quantities, PDG can be
# replaced with it.
_refs = dict(nist_codata="NIST CODATA 2010: "
             "http://physics.nist.gov/cuu/Constants/index.html",
             pdg_2011="PDG RPP 2011: http://pdg.lbl.gov/",
             ssd_p="http://ssd.jpl.nasa.gov/?planet_phys_par",
             ssd_c="http://ssd.jpl.nasa.gov/?constants",
             smith="Smith et. al. 1989AJ.....97..265S",
             sofac_20101201_tsc="SOFA Time Scale and Calendar Tools,"
                   " SOFA 20101201, C language, doc version 1.0",
             lieske="Lieske 1979A&A....73..282L")

# Base namedtuple.
Base = namedtuple("Base", "val name unit err src notes")


class AstroConst(Base):
    """Objects of this class behave like floats but have meta data.

    This class has a namedtuple as its base class. The idea is to have
    a type that behaves like a number in mathematical operations, but
    with additional meta data that a user can check when needed. The
    method data includes name, unit, uncertainities, source and notes.
    """
    def __new__(self, val, name=None, unit=None, err=None, src=None,
                 notes=None):
        return Base.__new__(self, val=val, name=name, unit=unit, err=err,
                      src=src, notes=notes)

    def __repr__(self):
        s = Base.__repr__(self)
        return s.replace("Base", "AstroConst")

    def __str__(self):
        return "{0}: {1} {2}".format(self.name, self.val, self.unit)

    def __add__(self, other):
        return self.val + other

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.val - other

    def __rsub__(self, other):
        return self.__sub__(other)

    def __mul__(self, other):
        return self.val * other

    def __rmul__(self, other):
        return self.__mul__(other)

    def __div__(self, other):
        return self.val / other

    def __rdiv__(self, other):
        return self.__div__(other)

    def __truediv__(self, other):
        return float(self.val) / other

    def __rtruediv__(self, other):
        return self.__truediv__(other)

    def __floordiv__(self, other):
        return self.val // other

    def __rfloordiv__(self, other):
        return self.__floordiv__(other)

    def __mod__(self, other):
        return self.val % other

    def __rmod__(self, other):
        return self.__mod__(other)

    def __divmod__(self, other):
        return divmod(self.val, other)

    def __rdivmod__(self, other):
        return self.__divmod__(other)

    def __pow__(self, other, z=None):
        return pow(self.val, other, z)

    def __lt__(self, other):
        return self.val < other

    def __le__(self, other):
        return self.val <= other

    def __eq__(self, other):
        return self.val == other

    def __ne__(self, other):
        return self.val != other

    def __gt__(self, other):
        return self.val > other

    def __ge__(self, other):
        return self.val >= other

    def __neg__(self):
        return -self.val

    def __pos__(self):
        return self.val

    def __abs__(self):
        return abs(self.val)


# Namedtuple for SSD supplied physical properties of planets.
SSDPlanet = namedtuple("SSDPlanet",
                       "erad mrad mass bd srp sop vmag geoalb"
                       " eqgrav escvel")

# Unit conversion.
convert_dict = {
    'length': {
        'm': (1.0, None),
        'cm': (1e+2, None),
        'mm': (1e+3, None),
        'mi': (1e+6, None),
        'nm': (1e+9, None),
        'a': (1e+10, None),
        'info': " m: meter\n cm: centi-meter\n mm: milli-meter\n"
                " mi: micro-meter\n nm: nano-meter\n a: angstrom"
    },
    'energy': {
        # Conversion factors from
        # http://physics.nist.gov/cuu/Constants/energy.html
        'j': (1.0, None),
        'ergs': (1e7, None),
        'ev': (6.24150934e+18, 0.00000014e+18),
        'mev': (6.24150934e+24, 0.00000014e+24),
        'kg': (1.112650056e-17, None),
        # Inverse meter = Joules / (h * c).
        'im': (5.03411701e+24, 0.00000022e+24),
        'icm': (5.03411701e+22, 0.00000022e+22),
        # Hertz = Joules / h.
        'hz': (1.509190311e+33, 0.000000067e+33),
        # Kelvin = Joules / k; k = bolztman constant
        'k': (7.2429716e+22, 0.0000066e+22),
        # Atomic mass unit = Joules / (atomic_mass_unit * c^2)
        'u': (6.70053585e+9, 0.00000030e+9),
        # Hartree = Joules / (Rydberg_infinity * h * c)
        'eh': (2.29371248e+17, 0.00000010e+17),
        'info': " j: Joules\n ergs: ergs\n ev: electron-volt\n"
                " mev: mega electron-volt\n kg: kilogram\n"
                " im: inverse meter\n icm: inverse centi-meter\n"
                " hz: Hertz\n k: Kelvin\n u: atmoic mass unit\n"
                " eh: Hartree\n"
        },
}


def convert_factor(q, from_u, to_u):
    """Return conversion factor between units.

    Parameters
    ----------
    q : str
        Type of physical quantity, for example "length". This must be a
        valid key of the `_convert_dict` dictionary defined in this
        module.
    from_u : str
        The initial unit. This must be a valid key of the
        sub-dictionary of `_convert_dict` corresponding to `q`.
    to_u : str
        The destination unit. This must be valid key of the
        sub-dictionary of `_convert_dict` corresponding to `q`.

    Returns
    -------
    x : float
        Factor for converting from `from_u` to `to_u`. That is quantity
        in `to_u` units can be obtained by multiplying the same in
        `from_u` with this conversion factor.

    Examples
    --------

    """
    m = convert_dict.get(q, None)
    if not m:
        raise ValueError("Unknow quantity %s" % q)
    return m[to_u][0] / m[from_u][0]


def convert_units(val, q, from_u, to_u):
    """Convert given quantity between given units."""
    # This will raise error if q or units are unknown.
    f = convert_factor(q, from_u=from_u, to_u=to_u)
    return val * f


# Planets and Pluto.
mercury = SSDPlanet(
    erad=AstroConst(val=2439.7, err=1.0,
                    unit="km", name="Equatorial radius",
                    notes="Equatorial radius.",
                    src=_refs['ssd_p']),
    mrad=AstroConst(val=2439.7, err=1.0,
                    unit="km", name="Mean radius", notes="Mean radius",
                    src=_refs['ssd_p']),
    mass=AstroConst(val=0.330104e24, err=0.000036e24,
                    unit="kg", name="Mass", notes="Mass.",
                    src=_refs['ssd_p']),
    bd=AstroConst(val=5.427, err=0.007, unit="g cm^-3",
                  name="Bulk density.", notes="Bulk density",
                  src=_refs['ssd_p']),
    srp=AstroConst(val=58.6462, err=None, unit="day",
                   name="Sidereal rotation period",
                   notes="Sidereal roation period.",
                   src=_refs['ssd_p']),
    sop=AstroConst(val=0.2408467, err=None, unit="year",
                   name="Sidereal orbit period",
                   notes="Sidereal orbit period.",
                   src=_refs['ssd_p']),
    vmag=AstroConst(val=-0.60, err=0.10, unit="mag",
                    name="Visual magnitude", notes="Visual magnitude",
                    src=_refs['ssd_p']),
    geoalb=AstroConst(val=0.106, err=None, unit=None,
                      name="Geometric Albedo",
                      notes="Geometric Albedo.",
                      src=_refs['ssd_p']),
    eqgrav=AstroConst(val=3.70, err=None, unit="m s^-2",
                      name="Equatorial gravity",
                      notes="Equatorial gravity.",
                      src=_refs['ssd_p']),
    escvel=AstroConst(val=4.25, err=None, unit="km s^-1",
                      name="Escape velocity", notes="Escape velocity.",
                      src=_refs['ssd_p'])
    )

venus = SSDPlanet(
    erad=AstroConst(val=6051.8, err=1.0,
                    unit="km", name="Equatorial radius",
                    notes="Equatorial radius.",
                    src=_refs['ssd_p']),
    mrad=AstroConst(val=6051.8, err=1.0,
                    unit="km", name="Mean radius", notes="Mean radius",
                    src=_refs['ssd_p']),
    mass=AstroConst(val=4.86732e24, err=0.00049e24,
                    unit="kg", name="Mass", notes="Mass.",
                    src=_refs['ssd_p']),
    bd=AstroConst(val=5.243, err=0.003, unit="g cm^-3",
                  name="Bulk density.", notes="Bulk density",
                  src=_refs['ssd_p']),
    srp=AstroConst(val=-243.018, err=None, unit="day",
                   name="Sidereal rotation period",
                   notes="Sidereal roation period.",
                   src=_refs['ssd_p']),
    sop=AstroConst(val=0.61519726, err=None, unit="year",
                   name="Sidereal orbit period",
                   notes="Sidereal orbit period.",
                   src=_refs['ssd_p']),
    vmag=AstroConst(val=-4.47, err=0.07, unit="mag",
                    name="Visual magnitude", notes="Visual magnitude",
                    src=_refs['ssd_p']),
    geoalb=AstroConst(val=0.65, err=None, unit=None,
                      name="Geometric Albedo",
                      notes="Geometric Albedo.",
                      src=_refs['ssd_p']),
    eqgrav=AstroConst(val=8.87, err=None, unit="m s^-2",
                      name="Equatorial gravity",
                      notes="Equatorial gravity.",
                      src=_refs['ssd_p']),
    escvel=AstroConst(val=10.36, err=None, unit="km s^-1",
                      name="Escape velocity", notes="Escape velocity.",
                      src=_refs['ssd_p'])
    )

earth = SSDPlanet(
    erad=AstroConst(val=6378.14, err=0.01,
                    unit="km", name="Equatorial radius",
                    notes="Equatorial radius.",
                    src=_refs['ssd_p']),
    mrad=AstroConst(val=6371.0, err=0.01,
                    unit="km", name="Mean radius", notes="Mean radius",
                    src=_refs['ssd_p']),
    mass=AstroConst(val=5.97219e24, err=0.00060e24,
                    unit="kg", name="Mass", notes="Mass.",
                    src=_refs['ssd_p']),
    bd=AstroConst(val=5.5134, err=0.0006, unit="g cm^-3",
                  name="Bulk density.", notes="Bulk density",
                  src=_refs['ssd_p']),
    srp=AstroConst(val=0.99726968, err=None, unit="day",
                   name="Sidereal rotation period",
                   notes="Sidereal roation period.",
                   src=_refs['ssd_p']),
    sop=AstroConst(val=1.0000174, err=None, unit="year",
                   name="Sidereal orbit period",
                   notes="Sidereal orbit period.",
                   src=_refs['ssd_p']),
    vmag=AstroConst(val=-3.86, err=None, unit="mag",
                    name="Visual magnitude", notes="Visual magnitude",
                    src=_refs['ssd_p']),
    geoalb=AstroConst(val=0.367, err=None, unit=None,
                      name="Geometric Albedo",
                      notes="Geometric Albedo.",
                      src=_refs['ssd_p']),
    eqgrav=AstroConst(val=9.80, err=None, unit="m s^-2",
                      name="Equatorial gravity",
                      notes="Equatorial gravity.",
                      src=_refs['ssd_p']),
    escvel=AstroConst(val=11.19, err=None, unit="km s^-1",
                      name="Escape velocity", notes="Escape velocity.",
                      src=_refs['ssd_p'])
    )

mars = SSDPlanet(
    erad=AstroConst(val=3396.19, err=0.1,
                    unit="km", name="Equatorial radius",
                    notes="Equatorial radius.",
                    src=_refs['ssd_p']),
    mrad=AstroConst(val=3389.50, err=0.2,
                    unit="km", name="Mean radius", notes="Mean radius",
                    src=_refs['ssd_p']),
    mass=AstroConst(val=0.641693e24, err=0.000064e24,
                    unit="kg", name="Mass", notes="Mass.",
                    src=_refs['ssd_p']),
    bd=AstroConst(val=3.9340, err=0.0008, unit="g cm^-3",
                  name="Bulk density.", notes="Bulk density",
                  src=_refs['ssd_p']),
    srp=AstroConst(val=1.02595676, err=None, unit="day",
                   name="Sidereal rotation period",
                   notes="Sidereal roation period.",
                   src=_refs['ssd_p']),
    sop=AstroConst(val=1.8808476, err=None, unit="year",
                   name="Sidereal orbit period",
                   notes="Sidereal orbit period.",
                   src=_refs['ssd_p']),
    vmag=AstroConst(val=-1.52, err=None, unit="mag",
                    name="Visual magnitude", notes="Visual magnitude",
                    src=_refs['ssd_p']),
    geoalb=AstroConst(val=0.150, err=None, unit=None,
                      name="Geometric Albedo",
                      notes="Geometric Albedo.",
                      src=_refs['ssd_p']),
    eqgrav=AstroConst(val=3.71, err=None, unit="m s^-2",
                      name="Equatorial gravity",
                      notes="Equatorial gravity.",
                      src=_refs['ssd_p']),
    escvel=AstroConst(val=5.03, err=None, unit="km s^-1",
                      name="Escape velocity", notes="Escape velocity.",
                      src=_refs['ssd_p'])
    )

jupiter = SSDPlanet(
    erad=AstroConst(val=71492, err=4,
                    unit="km", name="Equatorial radius",
                    notes="Equatorial radius.",
                    src=_refs['ssd_p']),
    mrad=AstroConst(val=69911, err=6,
                    unit="km", name="Mean radius", notes="Mean radius",
                    src=_refs['ssd_p']),
    mass=AstroConst(val=1898.13e24, err=0.19e24,
                    unit="kg", name="Mass", notes="Mass.",
                    src=_refs['ssd_p']),
    bd=AstroConst(val=1.3262, err=0.0004, unit="g cm^-3",
                  name="Bulk density.", notes="Bulk density",
                  src=_refs['ssd_p']),
    srp=AstroConst(val=0.41354, err=None, unit="day",
                   name="Sidereal rotation period",
                   notes="Sidereal roation period.",
                   src=_refs['ssd_p']),
    sop=AstroConst(val=11.862615, err=None, unit="year",
                   name="Sidereal orbit period",
                   notes="Sidereal orbit period.",
                   src=_refs['ssd_p']),
    vmag=AstroConst(val=-9.40, err=None, unit="mag",
                    name="Visual magnitude", notes="Visual magnitude",
                    src=_refs['ssd_p']),
    geoalb=AstroConst(val=0.52, err=None, unit=None,
                      name="Geometric Albedo",
                      notes="Geometric Albedo.",
                      src=_refs['ssd_p']),
    eqgrav=AstroConst(val=124.79, err=None, unit="m s^-2",
                      name="Equatorial gravity",
                      notes="Equatorial gravity.",
                      src=_refs['ssd_p']),
    escvel=AstroConst(val=60.20, err=None, unit="km s^-1",
                      name="Escape velocity", notes="Escape velocity.",
                      src=_refs['ssd_p'])
    )

saturn = SSDPlanet(
    erad=AstroConst(val=60268, err=4,
                    unit="km", name="Equatorial radius",
                    notes="Equatorial radius.",
                    src=_refs['ssd_p']),
    mrad=AstroConst(val=58232, err=6,
                    unit="km", name="Mean radius", notes="Mean radius",
                    src=_refs['ssd_p']),
    mass=AstroConst(val=568.319e24, err=0.057e24,
                    unit="kg", name="Mass", notes="Mass.",
                    src=_refs['ssd_p']),
    bd=AstroConst(val=0.6871, err=0.0002, unit="g cm^-3",
                  name="Bulk density.", notes="Bulk density",
                  src=_refs['ssd_p']),
    srp=AstroConst(val=0.44401, err=None, unit="day",
                   name="Sidereal rotation period",
                   notes="Sidereal roation period.",
                   src=_refs['ssd_p']),
    sop=AstroConst(val=29.447498, err=None, unit="year",
                   name="Sidereal orbit period",
                   notes="Sidereal orbit period.",
                   src=_refs['ssd_p']),
    vmag=AstroConst(val=-8.88, err=None, unit="mag",
                    name="Visual magnitude", notes="Visual magnitude",
                    src=_refs['ssd_p']),
    geoalb=AstroConst(val=0.47, err=None, unit=None,
                      name="Geometric Albedo",
                      notes="Geometric Albedo.",
                      src=_refs['ssd_p']),
    eqgrav=AstroConst(val=10.44, err=None, unit="m s^-2",
                      name="Equatorial gravity",
                      notes="Equatorial gravity.",
                      src=_refs['ssd_p']),
    escvel=AstroConst(val=36.09, err=None, unit="km s^-1",
                      name="Escape velocity", notes="Escape velocity.",
                      src=_refs['ssd_p'])
    )

uranus = SSDPlanet(
    erad=AstroConst(val=25559, err=4,
                    unit="km", name="Equatorial radius",
                    notes="Equatorial radius.",
                    src=_refs['ssd_p']),
    mrad=AstroConst(val=25362, err=7,
                    unit="km", name="Mean radius", notes="Mean radius",
                    src=_refs['ssd_p']),
    mass=AstroConst(val=86.8103e24, err=0.0087e24,
                    unit="kg", name="Mass", notes="Mass.",
                    src=_refs['ssd_p']),
    bd=AstroConst(val=1.270, err=0.001, unit="g cm^-3",
                  name="Bulk density.", notes="Bulk density",
                  src=_refs['ssd_p']),
    srp=AstroConst(val=-0.71833, err=None, unit="day",
                   name="Sidereal rotation period",
                   notes="Sidereal roation period.",
                   src=_refs['ssd_p']),
    sop=AstroConst(val=84.016846, err=None, unit="year",
                   name="Sidereal orbit period",
                   notes="Sidereal orbit period.",
                   src=_refs['ssd_p']),
    vmag=AstroConst(val=-7.19, err=None, unit="mag",
                    name="Visual magnitude", notes="Visual magnitude",
                    src=_refs['ssd_p']),
    geoalb=AstroConst(val=0.51, err=None, unit=None,
                      name="Geometric Albedo",
                      notes="Geometric Albedo.",
                      src=_refs['ssd_p']),
    eqgrav=AstroConst(val=8.87, err=None, unit="m s^-2",
                      name="Equatorial gravity",
                      notes="Equatorial gravity.",
                      src=_refs['ssd_p']),
    escvel=AstroConst(val=21.38, err=None, unit="km s^-1",
                      name="Escape velocity", notes="Escape velocity.",
                      src=_refs['ssd_p'])
    )

neptune = SSDPlanet(
    erad=AstroConst(val=24764, err=15,
                    unit="km", name="Equatorial radius",
                    notes="Equatorial radius.",
                    src=_refs['ssd_p']),
    mrad=AstroConst(val=24622, err=19,
                    unit="km", name="Mean radius", notes="Mean radius",
                    src=_refs['ssd_p']),
    mass=AstroConst(val=102.410e24, err=0.010e24,
                    unit="kg", name="Mass", notes="Mass.",
                    src=_refs['ssd_p']),
    bd=AstroConst(val=1.638, err=0.004, unit="g cm^-3",
                  name="Bulk density.", notes="Bulk density",
                  src=_refs['ssd_p']),
    srp=AstroConst(val=0.67125, err=None, unit="day",
                   name="Sidereal rotation period",
                   notes="Sidereal roation period.",
                   src=_refs['ssd_p']),
    sop=AstroConst(val=164.79312, err=None, unit="year",
                   name="Sidereal orbit period",
                   notes="Sidereal orbit period.",
                   src=_refs['ssd_p']),
    vmag=AstroConst(val=-6.87, err=None, unit="mag",
                    name="Visual magnitude", notes="Visual magnitude",
                    src=_refs['ssd_p']),
    geoalb=AstroConst(val=0.41, err=None, unit=None,
                      name="Geometric Albedo",
                      notes="Geometric Albedo.",
                      src=_refs['ssd_p']),
    eqgrav=AstroConst(val=11.15, err=None, unit="m s^-2",
                      name="Equatorial gravity",
                      notes="Equatorial gravity.",
                      src=_refs['ssd_p']),
    escvel=AstroConst(val=23.56, err=None, unit="km s^-1",
                      name="Escape velocity", notes="Escape velocity.",
                      src=_refs['ssd_p'])
    )

pluto = SSDPlanet(
    erad=AstroConst(val=1151, err=6,
                    unit="km", name="Equatorial radius",
                    notes="Equatorial radius.",
                    src=_refs['ssd_p']),
    mrad=AstroConst(val=1151, err=6,
                    unit="km", name="Mean radius", notes="Mean radius",
                    src=_refs['ssd_p']),
    mass=AstroConst(val=0.01309e24, err=0.00018e24,
                    unit="kg", name="Mass", notes="Mass.",
                    src=_refs['ssd_p']),
    bd=AstroConst(val=2.05, err=0.04, unit="g cm^-3",
                  name="Bulk density.", notes="Bulk density",
                  src=_refs['ssd_p']),
    srp=AstroConst(val=-6.3872, err=None, unit="day",
                   name="Sidereal rotation period",
                   notes="Sidereal roation period.",
                   src=_refs['ssd_p']),
    sop=AstroConst(val=247.92065, err=None, unit="year",
                   name="Sidereal orbit period",
                   notes="Sidereal orbit period.",
                   src=_refs['ssd_p']),
    vmag=AstroConst(val=-1.0, err=None, unit="mag",
                    name="Visual magnitude", notes="Visual magnitude",
                    src=_refs['ssd_p']),
    geoalb=AstroConst(val=0.3, err=None, unit=None,
                      name="Geometric Albedo",
                      notes="Geometric Albedo.",
                      src=_refs['ssd_p']),
    eqgrav=AstroConst(val=0.066, err=None, unit="m s^-2",
                      name="Equatorial gravity",
                      notes="Equatorial gravity.",
                      src=_refs['ssd_p']),
    escvel=AstroConst(val=1.23, err=None, unit="km s^-1",
                      name="Escape velocity", notes="Escape velocity.",
                      src=_refs['ssd_p'])
    )

# Constants.
# Physical constants.
c = AstroConst(val=299792458.0, err=None,
               unit="m/s", name="Speed of light",
               notes="Speed of light in vacuum.",
               src=_refs['nist_codata'])
light_speed = c
h = AstroConst(val=6.62606957e-34, err=0.00000029e-34,
               unit="J s", name="Planck constant",
               notes="Planck constant in Jules second.",
               src=_refs['nist_codata'])
planck_const = h
h_ev = AstroConst(val=4.135667516e-15, err=0.000000091e-15,
                  unit="ev s", name="Planck constant",
                  notes="Planck constant in electron volts second.",
                  src=_refs['nist_codata'])
h_bar = AstroConst(val=1.054571726e-34, err=0.000000047e-34,
                   unit="J s", name="Reduced Planck constant",
                   notes="h/2pi in Joules second.",
                   src=_refs['nist_codata'])
h_bar_ev = AstroConst(val=6.58211928e-16, err=0.00000015e-16,
                      unit="J s", name="Reduced Planck constant",
                      notes="h/2pi in electron volts second.",
                      src=_refs['nist_codata'])
k = AstroConst(val=1.3806488e-23, err=0.0000013e-23,
               unit="J K^-1", name="Boltzmann constant",
               notes="Boltzmann constant in Jules per kelvin.",
               src=_refs['nist_codata'])
boltzmann_const = k
e_mass = AstroConst(val=9.10938291e-31, err=0.00000040e-31,
                    unit="Kg", name="Electron mass",
                    notes="Electron mass in kilograms.",
                    src=_refs['nist_codata'])
e_mass_u = AstroConst(val=5.4857990946e-4, err=0.0000000022e-4,
                      unit="u", name="Electron mass in u.",
                      notes="Mass of electron in atomic mass units.",
                      src=_refs['nist_codata'])
p_mass = AstroConst(val=1.672621777e-27, err=0.000000074e-27,
                    unit="Kg", name="Proton mass",
                    notes="Proton mass in kilograms.",
                    src=_refs['nist_codata'])
p_mass_u = AstroConst(val=1.007276466812, err=0.000000000090,
                      unit="u", name="Proton mass in u",
                      notes="Proton mass in atomic mass units.",
                      src=_refs['nist_codata'])
sigma = AstroConst(val=5.670373e-8, err=0.000021e-8,
                   unit="W m^-2 K^-4",
                   name="Stefan-Boltzmann constant",
                   notes="Stefan-Boltzmann constant in SI units.",
                   src=_refs['nist_codata'])
wien_freq = AstroConst(val=5.8789254e10, err=0.0000053e10,
                       unit="Hz K^-1",
                       name="Wien frequency displacement law constant",
                       notes="F_max / T = wien_freq.",
                       src=_refs['nist_codata'])
wien_wave = AstroConst(val=2.8977721e-3, err=0.0000026e-3,
                       unit="m K",
                       name="Wien wavelength displacement law constant",
                       notes="Lambda_max * T = wien_wave.",
                       src=_refs['nist_codata'])

# Astronomical constants.
DJ = AstroConst(val=86400.0, err=None,
                unit="s", name="Julian day in seconds",
                notes="Seconds in a day.",
                src=_refs['ssd_c'])
YJ = AstroConst(val=365.25, err=None,
                unit="day", name="Julian year in days",
                notes="Days in a Julian year.",
                src=_refs['ssd_c'])
CJ = AstroConst(val=36525, err=None,
                unit="day", name="Julian century in days",
                notes="Days in a Julian century.",
                src=_refs['ssd_c'])
CB = AstroConst(val=36524.21987817305, err=None,
                unit="day", name="Besselian century in days.",
                notes="Days in a Besselian century (at 1900).",
                src=_refs['smith'])
MJD0 = AstroConst(val=2400000.5, err=None,
                   unit="day", name="MJD zero point",
                   notes="Zero point of Modified Julian Day system.",
                   src=_refs['sofac_20101201_tsc'])
MJD_J2000 = AstroConst(val=51544.5, err=None,
                   unit="day", name="MJD of J2000.0",
                   notes="MJD of January 1 12:00:00 TT.",
                   src=_refs['sofac_20101201_tsc']),
MJD_B1950 = AstroConst(val=33281.92345905, err=None,
                   unit="day", name="MJD of B1950",
                   notes="MJD of Besselian epoch B1950.0.",
                   src=_refs['lieske']),
MJD_B1900 = AstroConst(val=15019.81352, err=None,
                   unit="day", name="MJD of B1900",
                   notes="MJD of Besselian epoch B1900.0.",
                   src=_refs['lieske']),
AU = AstroConst(val=149597870700.0, err=3.0,
                unit="m", name="Astronomical unit",
                notes="1 astronomical unit.",
                src=_refs['pdg_2011'])
PC = AstroConst(val=3.0856776e16, err=None,
                unit="m", name="Parsec",
                notes="1 parsec.",
                src=_refs['pdg_2011'])
parsec = PC
SUN_M = AstroConst(val=1.9884e30, err=0.0002e30,
                   unit="Kg", name="Solar mass",
                   notes="1 solar mass in kilograms.",
                   src=_refs['pdg_2011'])
SUN_RAD = AstroConst(val=6.9551e8, err=0.0004e8,
                     unit="m", name="Solar equatorial radius",
                     notes="Equatorial radius of sun in meters.",
                     src=_refs['pdg_2011'])
SUN_L = AstroConst(val=3.8427e26, err=0.0014e26,
                   unit="J s", name="Solar luminosity",
                   notes="Bolometric luminosity of Sun in Watts.",
                   src=_refs['pdg_2011'])
