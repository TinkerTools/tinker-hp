!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     #################################################################
!     ##                                                             ##
!     ##  module units  --  physical constants and unit conversions  ##
!     ##                                                             ##
!     #################################################################
!
!
!     literature references:
!
!     P. J. Mohr, B. N. Taylor and D. B. Newell, "CODATA Recommended
!     Values of the Fundamental Physical Constants: 2010", Reviews of
!     Modern Physics, 84, 1527-1605 (2012)
!
!     Most values below are taken from 2010 CODATA reference values;
!     available on the web from the National Institute of Standards
!     and Technology at http://physics.nist.gov/constants/
!
!     The conversion from calorie to Joule is the definition of the
!     thermochemical calorie as 1 cal = 4.1840 J from ISO 31-4 (1992)
!
!     The "coulomb" energy conversion factor is found by dimensional
!     analysis of Coulomb's Law, ie, by dividing the square of the
!     elementary charge in Coulombs by 4*pi*eps0*rij, where eps0 is
!     the permittivity of vacuum (the "electric constant"); note that
!     eps0 is typically given in F/m, equivalent to C**2/(J-m)
!
!     The approximate value used for the Debye, 3.33564 x 10-30 C-m,
!     is from IUPAC Compendium of Chemical Technology, 2nd Ed. (1997)
!
!     The value of "prescon" is based on definition of 1 atmosphere
!     as 101325 Pa set by the 10th Conference Generale des Poids et
!     Mesures (1954), where a Pascal (Pa) is equal to a J/m**3
!
!
module units
   implicit none
   real*8 :: avogadro !<Avogadro's number (N) in particles/mole
   real*8 :: lightspd !<speed of light in vacuum (c) in cm/ps
   real*8 :: boltzmann !<Boltzmann constant (kB) in g*Ang**2/ps**2/mole/K
   real*8 :: gasconst !<ideal gas constant (R) in kcal/mole/K
   real*8 :: emass !<mass of an electron in atomic mass units
   real*8 :: planck !<Planck's constant (h) in J-s
   real*8 :: joule !<conversion from calories to joules
   real*8 :: convert !<conversion from kcal to g*Ang**2/ps**2
   real*8 :: bohr !<conversion from Bohrs to Angstroms
   real*8 :: hartree !<conversion from Hartree to kcal/mole
   real*8 :: evolt !<conversion from Hartree to electron-volts
   real*8 :: efreq !<conversion from Hartree to cm-1
   real*8 :: coulomb !<conversion from electron**2/Ang to kcal/mole
   real*8 :: debye !<conversion from electron-Ang to Debyes
   real*8 :: prescon !<conversion from kcal/mole/Ang**3 to Atm
   real*8 :: hbar_planck !<Planck's constant in g*Ang**2/ps**2/mole*ps
   real*8 :: cm1 !<conversion from ps-1 to cm-1
   parameter (avogadro=6.02214129d+23)
   parameter (lightspd=2.99792458d-2)
   parameter (boltzmann=0.831446215d0)
   parameter (gasconst=1.98720415d-3)
   parameter (emass=5.485799095d-4)
   parameter (planck=6.62606957d-34)
   parameter (joule=4.1840d0)
   parameter (convert=4.1840d+2)
   parameter (bohr=0.52917721092d0)
   parameter (hartree=627.5094743d0)
   parameter (evolt=27.21138503d0)
   parameter (efreq=2.194746313708d+5)
   parameter (coulomb=332.063714d0)
   parameter (debye=4.80321d0)
   parameter (prescon=6.85684112d+4)
   parameter (hbar_planck=(planck*1.d11*avogadro)/(2d0*acos(-1.d0)))
   parameter (cm1=1d0/lightspd/(2d0*acos(-1.d0)))
   save
end
