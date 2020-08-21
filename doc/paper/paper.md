---
title: 'Pencil Code, a modular MPI code for partial differential equation: multipurpose and multiuser-maintained'
tags:
 - Fortran90
 - fluid dynamics
 - magnetohydrodynamics
 - Python
 - IDL
 - astrophysics
 - radiation
 - inertial particles
 - combustion
authors:
 - name: The Pencil Code Collaboration
   affiliation: 1
 - name: Axel Brandenburg
   affiliation: "1, 2, 3"
   orcid: 0000-0002-7304-021X
 - name: Anders Johansen
   affiliation: 4
 - name: Philippe A. Bourdin
   affiliation: 5
 - name: Wolfgang Dobler
   affiliation: 6
 - name: Wladimir Lyra
   affiliation: 7
 - name: Matthias Rheinhardt
   affiliation: 8
 - name: Frederick Gent
   affiliation: "8, 9"
   orcid: 0000-0002-1331-2260
 - name: Jörn Warnecke
   affiliation: 10
   orcid: 0000-0002-9292-4600
 - name: Nils Erland L. Haugen
   affiliation: 11
   orcid: 0000-0002-9184-8722
 - name: Add Yourself
   affiliation: 12
affiliations:
 - name: Nordita, KTH Royal Institute of Technology and Stockholm University
   index: 1
 - name: Department of Astronomy, Stockholm University
   index: 2
 - name: McWilliams Center for Cosmology & Department of Physics, Carnegie Mellon University
   index: 3
 - name: Lund Observatory, Department of Astronomy and Theoretical Physics, Lund University
   index: 4
 - name: Space Research Institute, Graz
   index: 5
 - name: Bruker, Potsdam
   index: 6
 - name: New Mexico State University
   index: 7
 - name: Astroinformatics, Department of Computer Science, Aalto University
   index: 8
 - name: School of Mathematics, Statistics and Physics, Newcastle University 
   index: 9
 - name: Max Planck Institute for Solar System Research 
   index: 10
 - name: SINTEF Energy Research, Trondheim, Norway
   index: 11
 - name: twelve
   index: 12

date: 20 August 2020
bibliography: paper.bib
---

# Summary

The Pencil Code is a highly modular physics-oriented simulation code
that can be adapted to a wide range of applications.
It is primarily designed to solving partial differential equations
(PDEs) of compressible magnetohydrodynamics [@2010ascl.soft10060B],
but the framework is general and can also be applied to situations not
related to PDEs, for example when just the message passing interface or
input/output strategies of the code are to be used.
The code can also solve for Lagrangian (inertial and noninertial)
particles, their coagulation and condensation, as well as their
interaction with the fluid.

The code is being used for Cartesian, cylindrical, and spherical geometries,
but further extensions are possible.
<!--- AB: comment for Joern, and commented out Yin--Yang
A preliminary implementation of a Yin--Yang mesh is also in place.
--->
One can choose between different time stepping schemes and different
spatial derivative operators.
There are also different diffusion operators to allow for both direct numerical
simulations (DNS) and various types of large-eddy simulations (LES).

# High-level functionality

An idea about the range of available modules can be obtained by inspecting
the examples under pencil-code/samples/, which include Cartesian and
spherical shell convection with Kramers opacities, 
solar atmosphere with radiation and ionization, M dwarfs as a star-in-a-box, 
coronae with heatflux and Boris correction, kinematic dynamos,
helical MHD turbulence, implicit resistivity, interlocked fluxrings,
cylindrical global disks with a dead zone, photoelectric fluids, 
dust in turbulent global disks, cosmic rays, gravitational waves
in the early universe, interstellar supernova-driven turbulence,
inertial particles, potential field boundary condition, sedimentation, sink particles, 
superparticles for condensation-coagulation, nonlinear test fields, ambipolar diffusion,
H2 flame speed tests, reversed field pinch, Boussinesq convection, chiral
dynamos, cylinder deposition with overlayed grid, dynamical alpha effect,
fargo, Kelvin-Helmholtz instability, resolved reacting char particles, shallow water.

Below a list of application that are described in dedicated papers:

* Coagulation and condensation in turbulence [@2008A&A.486.597J; @2017JAMES.9.1116L]
* Radiative transfer [@2006A&A.448.731H; @2014A&A.571A.68B; @2020GApFD.114.162B]
* Chiral magnetic effect in relativistic plasmas [@2018ApJ.858.124S]
* Primordial gravitational waves [@2020GApFD.114.130R]
* Gaseous combustion and detonation [@2011JCoPh.230.1B; @Zhang_etal_2020comb]
* Burning particles, resolved or un-resolved [@2020GApFD.114.58Q; @2017CNF.185a160]
* Flows around immersed solid objects [@2019IJCFD.33.43A; @2020GApFD.114.35A; @2010JFM.661a239]
* Test-field method for turbulent MHD transport [@2010A&A.520A.28R; @2010PhST.142a4028B; @2018A&A.609A.51W]
* Spherical shell dynamos and convection [@2009ApJ.697.923M; @2020GApFD.114.8K]
* Boris correction for coronal physics [@2020GApFD.114.213C]
* Thermal instability and mixing [@2012ApJ.758.48Y]
* Implicit solver for temperature [@2008A&A.484.29G]

# Statement of need and purpose of software

The code provides an easily adaptable tool for solving both standard
equations as well as others, such as the test-field equations.
Significant amounts of runtime diagnostics are available.
Users have access to the latest development version and can ask to
join the circle of developers.
Every revision on GitHub is verified on 9 tests on travis-ci.com.
The current version is also automatically being tested on 59 hourly
tests and on 79 daily tests.
Continuous progress on the code is driven by the research of the
individual developers.
Among the currently 83 developers with check-in permission, there are
currently 16 owners who can give others check-in permission.
Of the developers, 34 have done more than 34 commits.
Further developments and interactions between developers and users is
being promoted through annual user meetings since 2004 and a newsletters
since 2020.
Since 2016, a steering committee of five elected owners reviews the
progress and can take decisions of general concern to the Pencil Code
community.

# Ongoing research using the Pencil Code

Current research topics with the code includes
interstellar and intercluster medium as well as early Universe;
small-scale dynamos and reconnection;
primordial magnetic fields and decaying turbulence;
gravitational waves from turbulent sources;
planet formation and inertial particles;
accretion discs and shear flows;
coronal heating and coronal mass ejections;
helical dynamos, helical turbulence, and catastrophic quenching;
helioseismology;
strongly stratified MHD turbulence and negative effective magnetic pressure instability;
convection in Cartesian domains;
global convection and dynamo simulations;
turbulent transport and test-field methods;
hydrodynamic and MHD instabilities;
chiral MHD;
Hydrodynamic and MHD turbulence;
turbulent gaseous and solid combustion, particle clustering and deposition on solid walls,
front propagation, radiation & ionization.
As of July 2020, 564 papers have been published that acknowledge use of
the Pencil Code [@zenodo.3466444].

# Key references

* Planet and planetesimal formation [@2007Natur.448.1022J; @2009A&A.497.869L]
* Solar dynamo [@2005ApJ.625.539B]
* MHD turbulence and cascades [@2004PhRvE.70a6308H]
* Fully convective star in a box [@2006ApJ.638.336D]
* Dynamo wave in spherical shell convection [@2012ApJ.755L.22K; @2014ApJ.796L.12W]
* Turbulent diffusivity quenching with test fields [@2008ApJ.676.740B; @2014ApJ.795.16K]
* Global simulations of debris disks [@2013Natur.499.184L]
* Supernova-driven turbulence [@2013MNRAS.432.1396G]
* Clustering in in supersonic turbulence [@2019MNRAS.483.5623M]
* Solar corona above active regions [@2011A&A.530A.112B; @2013A&A.555A.123B]
* Stratified shearing box simulations [@2011ApJ.740.18O]

# Acknowledgements

We acknowledge contributions from all submitters and their supporting
funding agencies.
In particular, we mention the Swedish Research Council,
grants 2012-5797, 2013-03992, 2017-03865, and 2019-04234,
the National Science Foundation under the grant AAG-1615100,
the FRINATEK grant 231444 under the Research Council of Norway, SeRC,
the grant "Bottlenecks for particle growth in turbulent aerosols"
from the Knut and Alice Wallenberg Foundation, Dnr.\ KAW 2014.0048,
and the University of Colorado through its support of the
George Ellery Hale visiting faculty appointment,

# References

