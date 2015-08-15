# pressure_drop.R v0.1
# Copyright 2015 Jim Hargreaves

# Calculates pressure drop for a fluid flowing in a pipe.
# Includes functions for calculating pump power and checking NPSHa.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>

# TODO
# 
# 1. Find a more authoratative source (peer reviewed) for calculating
#    friction factor under transition flow conditions.
#
# 2. Implement sense-checking function on CalcFluidVel. e.g.
#    Warn if fluid velocity will cause water hammer.


# 1. Parameters & Constants

# 1.1 Process Material Properties & Flowrates
  vol.flowrate       <- 3.00     # Volumetric Flowrate (m3 hr-1)
  fluid.dens         <- 2016.61  # Fluid density (kg m-3)
  fluid.visc         <- 40E-3    # Fluid viscosity (Pa s)


# 1.2 Pipe Geometry & Fittings
  # For details refer to A-26 to A-30 of Crane TP-410

  # 1.2.1 Pipe Geometry
    pipe.dia         <- 0.025    # Pipe inside diameter (m)
    pipe.roughness   <- 0.15E-3  # Mean height of roughness (m)
    pipe.len         <- 62.616   # Pipe Length (m)

  # 1.2.2 Quantity of Bends
    lr.bend.qty      <- 17       # Quantity of long-radius bends

  # 1.2.3 Quantity of Valves
    fb.gate.val.qty  <- 1        # Quantity of full-bore gate valves
    
    ball.val.qty     <- 1        # Quantity of ball valves
    gate.val.qty     <- 0        # Qty of gate valves

  # 1.2.4 Contractions & Enlargements
    sud.enlarg.qty   <- 1        # Qty of sudden enlargements (theta <45 deg)
    grad.enlarg.qty  <- 0        # Qty of gradual enlargements (45 <theta <180)
    sud.cont.qty     <- 0        # Qty of sudden contractions (theta <45 deg)
    grad.cont.qty    <- 0        # Qty of gradual contractions (45 <theta <180)

# 1.3 Constants  & First Guesses
  g                  <- 9.81     # Acceleration due to gravity
  darcy.guess        <- 0.00500  # First Guess for Darcy Friction Factor


# 2. Function Definitions

CalcReynoldsNum <- function(fluid.dens, fluid.vel, pipe.dia, fluid.visc) {  
  # Calculates the Reynolds Number for a fluid flowing in a pipe
  # 	
  # Args:
  # fluid.dens - density of fluid (kg m-3)
  # fluid.vel - velocity of fluid (m s-1)
  # pipe.dia - pipe inside diameter (m)
  # fluid.visc - viscosity of fluid (Pa s)
  #
  # Returns:
  # reynolds - Reynolds number of the fluid in the pipe
  
  reynolds.num <- (fluid.dens * fluid.vel * pipe.dia) / fluid.visc
  return(reynolds.num)
}

CalcDarcyFactor <- function(reynolds.num, pipe.dia, pipe.roughness, darcy.guess) {
  # Calculates the Darcy Friction Factor for fluid flow in a pipe
  # Uses Fd = 64 / Re for laminar flow (Re < 2320)
  # Uses a correlation for Transitional Flow ( 2321 < Re < 4000)
  # Uses Colbrook (solved by Newton-Rhapson) for laminar flow (Re > 4000)
  #
  # Args:
  # reynolds.num - The reynolds number of the fluid in a pipe
  # pipe.dia - pipe inside diameter (m)
  # pipe.roughness - mean height of pipe surface roughness (m)
  # darcy.guess - Initial guess at Darcy Factor
  # 
  # Returns:
  # darcy.factor - Darcy friction factor

  if (reynolds.num <= 2320) {
    #Laminar Flow

    darcy.factor <- 64 / reynolds.num
	
  } else if (2320 < reynolds.num & reynolds.num <= 4000) {
      # Transitional Flow
      # N.B. Formula in Morrison's paper originally gave Fanning, not Darcy.
      # TODO(JH): Find a more rigorous (peer-reviewed) correlation.

      darcy.factor <- {(0.0304 * ((3170 / reynolds.num)^0.165)
        / (1 + (3170 / reynolds.num)^7.0)) + (64 / reynolds.num)}

  } else {
    if (4000 < reynolds.num)
      #Turbulent Flow
      
      # Substituting some constants to simplify
      roughness.term <- pipe.roughness / (3.7 * pipe.dia)
      reynolds.term <- 2.51 / reynolds.num

      for (i in 1:50) {

        # Residual form of Colbrook
        f <- {darcy.guess^(-0.5) + 2 *
        log10(roughness.term + reynolds.term * darcy.guess^(-0.5))
        }

        # First Derivative of Colbrook
        f1 <- {-0.5 * darcy.guess^(-1.5) * (1 + (2 * reynolds.term)
        / (log(10) * (roughness.term + reynolds.term * darcy.guess^(-0.5))))
        }
        darcy.guess <- darcy.guess - f/f1
      }
     darcy.factor <- darcy.guess 
  
  }
  return(darcy.factor)
}

CalcFluidVelocity <- function(vol.flowrate, pipe.dia) {
  # Given pipe diameter and flowrate, calculates fluid velocity
  #
  # Args:
  # vol.flowrate - Volumetric Flowrate (m3 h-1)
  # pipe.dia - Pipe inside diameter (m)
  # 
  # Returns:
  # fluid.vel - Linear velocity of the fluid (m s-1)

  fluid.vel <- vol.flowrate / (pi * (pipe.dia / 2) ^ 2 * 3600)
  return(fluid.vel)
}

CalcFricH <- function(fluid.vel, darcy.factor, pipe.len, pipe.dia, k.fit) {
  # Calculates the frictional losses for a fluid flowing in a length of pipe
  #
  # Args:
  # fluid.vel - linear velocity of the fluid (m s-1)
  # darcy.factor - darcy friction factor
  # pipe.len - Pipe length (m)
  # pipe.dia - inside diameter of pipe (m)
  # k.fit - Sum of fitting Resistance cefficients 
  #
  # Returns:
  # fric.h

  fric.h <- {(fluid.vel^2 / (2 * g)) * 
  (darcy.factor * (pipe.len / pipe.dia) + k.fits)
  }
  return(fric.h)
}

#CalcFittingK <- function(fitting.type, d1, d2, theta) {
  # Calculates the k factor of a fitting as per Crane TP-410
  #
  # Args:
  # fitting.type - Type of fitting
  # d1 - Smallest Diameter (m)
  # d2 - Largest Diameter (m)
  # theta - critical angle (degrees)
  #
  # Returns:
  # k.fit - Calculated K factor

#  beta <- (d1/d2)
#  crane.f1 <- (0.8 * sin(theta / 2) * (1 - beta^2)) / beta^4
#  crane.f2 <- {}
#  crane.f3 <- {}
#  crane.f4 <- {}
#  crane.f5 <- {}
#  crane.f6 <- {}
#  crane.f7 <- {}

#}




 


CalcPumpPower <- function(Q, fluid.dens, h, efficiency) {
  # Calculates the power requirement of a centrifugal liquid pump
  #
  # Args:
  # Q - Flowrate (Units?)
  # fluid density - density of fluid
  # h - total differential head
  # efficiency of the pump
  #
  # Returns:
  # abs.power - The estimated absorbed power of the pump.
  
  hyd.power <- (Q * fluid.dens * g * h) / 3.6E6
  abs.power <- hyd.power / efficiency
  return(abs.power)
}


# 3. Executed Statements

fluid.vel <- CalcFluidVelocity(vol.flowrate, pipe.dia)

reynolds.num <- CalcReynoldsNum(fluid.dens, fluid.vel, pipe.dia, fluid.visc)

darcy.factor <- CalcDarcyFactor(reynolds.num, pipe.roughness, darcy.guess)

#fric.h <- CalcFricH(fluid.vel, darcy.factor, pipe.len, pipe.dia, k.fits)

print(fluid.vel)
print(reynolds.num)
print(darcy.factor)


#  Appendix 1 - References

# Crane Valves North America (1999) TP-410M Flow of Fluids through Valves, 
# Fittings and Pipe. Texas: Crane Ltd 
# 
# F.A. Morrison, Data Correlation for Friction Factor in Smooth Pipes. Dept of
# Chemical Engineering, Michigan Technological University. Houston, MI.

