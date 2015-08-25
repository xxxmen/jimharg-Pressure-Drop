# pressure_drop.R v0.2
# Copyright 2015 Jim Hargreaves

# Calculates pressure drop for an incompressible fluid flowing in a pipe.
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
# 1. Formalise output reporting. Different modes? Tables? What options exist?
#
# 2. Way of calculating k.fit to be considered.


# 1. Parameters & Constants

# 1.1 Process Material Properties
  vol.flowrate       <- 3.000      # Volumetric Flowrate (m3 hr-1)
  fluid.dens         <- 1900       # Fluid density (kg m-3)
  fluid.visc         <- 40E-3      # Fluid viscosity (Pa s)
  fluid.vapour.press <- 3169.9     # Absolute vapor pressure of fluid (Pa)

# 1.2 System Properties

  # 1.2.1 Pipe Geometry
    pipe.dia         <- 0.025      # Pipe inside diameter (m)
    pipe.roughness   <- 0.0015E-3  # Mean height of roughness (m)
    pipe.len         <- 62.616     # Pipe Length (m)

  # 1.2.2 System Geometry & Parameters
    suction.stat.h   <- 1.276      # Suction Static Head (m)
    discharge.stat.h <- 3.367      # Discharge Static Head (m)
    suction.tank.p   <- 101.315E3  # absolute pressure in suction tank (Pa)

  # 1.2.3 Pump Parameters
    pump.eff         <- 0.25       # Estimated efficiency of pump


# 1.3 Constants  & First Guesses
  g                  <- 9.81       # Acceleration due to gravity (m s-2)

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
  
  if(!is.numeric(fluid.dens)) 
    stop("One or more elements of fluid.dens are non-numeric.")
  if(!is.numeric(fluid.vel))
    stop("One or more elements of fluid.vel are non-numeric.")
  if(!is.numeric(pipe.dia))
    stop("One or more elements of pipe.dia are non-numeric.")
  if(!is.numeric(fluid.visc))
    stop("One or more elements of fluid.visc are non-numeric.")

  reynolds.num <- (fluid.dens * fluid.vel * pipe.dia) / fluid.visc
  return(reynolds.num)
}

CalcDarcyFactor<- function(reynolds.num, pipe.dia, pipe.roughness, 
                               darcy.guess = 0.005, darcy.iterations = 50) {
  # Calculates the Darcy Friction Factor for fluid flow in a pipe
  # Uses Fd = 64 / Re for laminar flow (Re < 2320)
  # Uses Colbrook (solved by Newton-Rhapson) for non-laminar flow (Re > 2320)
  #
  # Args:
  # reynolds.num - The reynolds number of the fluid in a pipe
  # pipe.dia - pipe inside diameter (m)
  # pipe.roughness - mean height of pipe surface roughness (m)
  # darcy.guess - Initial guess at Darcy Factor (optional, default = 0.005)
  # darcy.iterations - Number of iterations. (optional, default = 50)
  # 
  # Returns:
  # darcy.factor - Darcy friction factor

  if(!is.numeric(reynolds.num)) 
    stop("One or more elements of reynolds.num are non-numeric.")
  if(!is.numeric(pipe.dia))
    stop("One or more elements of pipe.dia are non-numeric.")
  if(!is.numeric(pipe.roughness))
    stop("One or more elements of pipe.roughness are non-numeric.")
  if(!is.numeric(darcy.guess))
    stop("One or more elements of darcy.guess are non-numeric.")
  if(!is.numeric(darcy.iterations))
    stop("One or more elements of darcy.iterations are non-numeric.")

  laminar.darcy <- 64 / reynolds.num
  
  roughness.term <- pipe.roughness / (3.7 * pipe.dia)
  
  reynolds.term <- 2.51 / reynolds.num

  for (i in 1:darcy.iterations) {
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
  turbulent.darcy <- darcy.guess

  darcy.factor <- ifelse(reynolds.num < 2320, laminar.darcy, turbulent.darcy)
  
  return(darcy.factor)
}

CalcFluidVel <- function(vol.flowrate, pipe.dia) {
  # Given pipe diameter and flowrate, calculates fluid velocity
  #
  # Args:
  # vol.flowrate - Volumetric Flowrate (m3 h-1)
  # pipe.dia - Pipe inside diameter (m)
  # 
  # Returns:
  # fluid.vel - Linear velocity of the fluid (m s-1)

  if(!is.numeric(vol.flowrate)) 
    stop("One or more elements of vol.flowrate are non-numeric.")
  if(!is.numeric(pipe.dia))
    stop("One or more elements of pipe.dia are non-numeric.")

  fluid.vel <- vol.flowrate / (pi * (pipe.dia / 2) ^ 2 * 3600)
  return(fluid.vel)
}

CalcFricH <- function(fluid.vel, darcy.factor, pipe.len, pipe.dia, k.fits) {
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

  if(!is.numeric(fluid.vel))
    stop("One or more elements of fluid.vel are non-numeric.")
  if(!is.numeric(darcy.factor))
    stop("One or more elements of darcy.factor are non-numeric.")
  if(!is.numeric(pipe.len))
    stop("One or more elements of pipe.len are non-numeric.")
  if(!is.numeric(pipe.dia))
    stop("One or more elements of pipe.dia are non-numeric.")
  if(!is.numeric(k.fits))
    stop("One or more elements of k.fits are non-numeric.")

  fric.h <- {(fluid.vel^2 / (2 * g)) * 
  (darcy.factor * (pipe.len / pipe.dia) + k.fits)
  }
  return(fric.h)
}

ConvertH2P <- function(head, fluid.dens) {
  # Converts fluid heads (m) to equivalent pressures (Pa)
  # 
  # Args:
  # head - fluid head (m)
  # fluid.dens - density of fluid (kg m-3)
  #
  # Returns:
  # pressure - pressure equivalent to head (Pa)

  if(!is.numeric(head)) 
    stop("One or more elements of head are non-numeric.")
  if(!is.numeric(fluid.dens))
    stop("One or more elements of fluid.dens are non-numeric.")

  pressure <- head * fluid.dens * g
  return(pressure)
}

ConvertP2H <- function(pressure, fluid.dens) {
  # Converts pressures (Pa) to equivalent fluid heads (m)
  # 
  # Args:
  # pressure - pressure equivalent to head (Pa)
  # fluid.dens - density of fluid (kg m-3)
  #
  # Returns:
  # head - fluid head (m)
  if(!is.numeric(pressure)) 
    stop("One or more elements of pressure are non-numeric.")
  if(!is.numeric(fluid.dens))
    stop("One or more elements of fluid.dens are non-numeric.")

  head <- pressure / (fluid.dens * g)
  return(head)
}

CalcPumpPower <- function(vol.flowrate, fluid.dens, total.diff.h, efficiency) {
  # Calculates the power requirement of a centrifugal liquid pump
  #
  # Args:
  # vol.flowrate - Flowrate (m3 h-1)
  # fluid density - density of fluid
  # total.diff.h - total differential head
  # pump.eff - efficiency of the pump
  #
  # Returns:
  # abs.power - The estimated absorbed power of the pump.
  if(!is.numeric(vol.flowrate)) 
    stop("One or more elements of vol.flowrate are non-numeric.")
  if(!is.numeric(fluid.dens))
    stop("One or more elements of fluid.dens are non-numeric.")
  if(!is.numeric(total.diff.h)) 
    stop("One or more elements of total.diff.h are non-numeric.")
  if(!is.numeric(efficiency))
    stop("One or more elements of efficiency are non-numeric.")
 
  hyd.power <- (vol.flowrate * fluid.dens * g * total.diff.h) / 3.6E6
  abs.power <- hyd.power / pump.eff
  return(abs.power)
}


# 3. Executed Statements
# (Customise to suit application)

fluid.vel <- CalcFluidVel(vol.flowrate, pipe.dia)

reynolds.num <- CalcReynoldsNum(fluid.dens, fluid.vel, pipe.dia, fluid.visc)

darcy.factor <- seq(1, length(reynolds.num), 1) 
darcy.factor <- CalcDarcyFactor(reynolds.num, pipe.dia, pipe.roughness)

k.fits <- 584 * darcy.factor + 2.434
# Determined Manually using Crane TP-410

discharge.fric.h <- CalcFricH(fluid.vel, darcy.factor, pipe.len, pipe.dia, 
                              k.fits)
discharge.fric.p <- ConvertH2P(discharge.fric.h, fluid.dens)

suction.fric.h <- 0
suction.fric.p <- 0
# Suction frictional losses assumed zero (short run)

total.diff.h <- (discharge.stat.h - suction.stat.h + discharge.fric.h + 
                 suction.fric.h)
total.diff.p <- ConvertH2P(total.diff.h, fluid.dens)


abs.power <- CalcPumpPower(vol.flowrate, fluid.dens, total.diff.h, pump.eff)

npsha <- (suction.tank.p + ConvertH2P(suction.stat.h, fluid.dens) - 
          suction.fric.p - fluid.vapour.press)

results <- data.frame(vol.flowrate, fluid.dens, fluid.visc, pipe.dia, pipe.len,
                      suction.stat.h, discharge.stat.h, suction.tank.p, 
                      pump.eff, fluid.vel, reynolds.num, darcy.factor,
                      discharge.fric.h, total.diff.h, abs.power, npsha) 

results.summary <- data.frame(signif(fluid.vel, 4), round(reynolds.num, 0), 
                      signif(darcy.factor, 3), round(discharge.fric.h, 3),
                      round(total.diff.h, 3), round(abs.power, 2), 
                      round(ConvertP2H(npsha, fluid.dens),3 ))
colnames(results.summary) <- c(
                       "Fluid Velocity (m s-1)",
                       "Reynolds No.",
                       "Darcy Friction Factor",
                       "discharge friction head (m)",
                       "total differential head (m)",
                       "absorbed power (kW)",
                       "NPSHa (m)")

print(results.summary)

#  Appendix 1 - References

# Crane Valves North America (1999) TP-410M Flow of Fluids through Valves, 
# Fittings and Pipe. Texas: Crane Ltd 
# 
# F.A. Morrison, Data Correlation for Friction Factor in Smooth Pipes. Dept of
# Chemical Engineering, Michigan Technological University. Houston, MI.

