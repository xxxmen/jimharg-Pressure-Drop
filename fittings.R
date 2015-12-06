# Part of the Pressure Drop Calculator Program
# Copyright 2015 Jim Hargreaves

# Calculates Fitting resistance coefficients in terms of the Darcy Friction
# Factor.

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

# I. ABBREVIATIONS
# The following abbreviations are used throughout:

# srb - Short radius bend (centre to end dim = 1d)
# lrb - Long Radius Bend (Centre to end dim = 1.5d)
# tbr - Tee, flow thru' branch
# tst - Tee, flow thru' straight
# btf - Butterfly Valve
# bav - Ball Valve
# 

# Counts - Number of each fittings present
srb.count <- 0
lrb.count <- 0
tbr.count <- 0
tst.count <- 0
btf.count <- 0
bav.count <- 0

# Other variables
btf.size <- 0  # Nominal inside diameter (in mm) of butterfly valves
bav.beta <- 0  # Ratio of ball valves throat dia to pipe dia

# Calculation of K values
srb.k <- srb.count * 20
lrb.k <- lrb.count * 14
tbr.k <- tbr.count * 60
tst.k <- tst.count * 20

CalcBtfK <- function(btf.count, btf.size) {
  # Calculates K resistance coefficient for butterfly valves
  #
  # Args:
  # btf.count - No of butterfly valves
  # btf.size - Nominal ID of butterfly valves in mm
  #
  # Returns:
  # btf.k - resistance coefficient for butterfly valves of that size
  
  
  btf.k.200 <- btf.count * 45
  btf.k.350 <- btf.count * 35
  btf.k.600 <- btf.count * 25

  btf.k <- { ifelse(btf.size <= 200, btf.k.200, 
                    ifelse( btf.size <= 350, btf.k.350, btf.k.600))
  }
  return(btf.k)
}

CalcBavK <- function(bav.count, bav.beta) {
  # Calculates K resistance coefficient for ball valves
  #
  # Args:
  # bav.count - No of ball valves
  # bav.beta - ratio of throat dia to pipe dia
  #
  # Returns:
  # bav.k - resistance coefficient for ball valves of that beta
  
  bav.k1 <- 
  bav.k.fb <- bav.count * 3 # K if full-bore
  bav.k.f5 <- 


#  Appendix 1 - References

# Crane Valves North America (1999) TP-410M Flow of Fluids through Valves, 
# Fittings and Pipe. Texas: Crane Ltd 

