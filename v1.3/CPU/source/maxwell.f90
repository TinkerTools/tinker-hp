!
!     Sorbonne University
!     Washington University in Saint Louis
!     University of Texas at Austin
!
!     ##################################################################
!     ##                                                              ##
!     ##  function maxwell  --  Maxwell-Boltzmann distribution value  ##
!     ##                                                              ##
!     ##################################################################
!
!
!     "maxwell" returns a speed in Angstroms/picosecond randomly
!     selected from a 3-D Maxwell-Boltzmann distribution for the
!     specified particle mass and system temperature
!
!     literature reference:
!
!     P. W. Atkins, "Physical Chemistry, 4th Edition", W. H. Freeman,
!     New York, 1990; see section 24.2 for general discussion
!
!
!> @brief 
!> returns a speed in Angstroms/picosecond randomly
!> selected from a 3-D Maxwell-Boltzmann distribution for the
!> specified particle mass and system temperature
!> @param[in] mass: mass 
!> @param[in] temperature: temperature in K
function maxwell (mass,temper)
   use units
   implicit none
   real*8 maxwell
   real*8 mass,temper
   real*8 rho,beta
   real*8 random,erfinv
   real*8 xspeed,yspeed
   real*8 zspeed
   external erfinv
!
!
!     set normalization factor for cumulative velocity distribution
!
   beta = sqrt(mass / (2.0d0*boltzmann*temper))
!
!     pick a randomly distributed velocity along each of three axes
!
   rho = random ()
   xspeed = erfinv(rho) / beta
   rho = random ()
   yspeed = erfinv(rho) / beta
   rho = random ()
   zspeed = erfinv(rho) / beta
!
!     set the final value of the particle speed in 3-dimensions
!
   maxwell = sqrt(xspeed**2 + yspeed**2 + zspeed**2)
   return
end
