# Oriface Plate.R v0.1
# Part of the Pressure Drop Calculator Program
# Copyright 2015 Jim Hargreaves

# Calculates the pressure drop across an oriface plate.

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

vol.flowrate       <- 408.0      # Volumetric Flowrate (m3 hr-1)
fluid.dens         <- 997.01     # Fluid density (kg m-3)
discharge.coef     <- 0.3        # Dimensionless discharge co-efficient
pipe.dia           <- 0.1485     # Pipe inside diameter (m)
throat.dia         <- 0.0988     # Oriface throat diamter (m)

#app.vel.factor <- sqrt(1 - (throat.dia/pipe.dia)^4)
app.vel.factor <- 1

pipe.area <- (pipe.dia / 2)^2 * pi

pressure.drop <- {(fluid.dens / 2) * 
(((vol.flowrate / 3600) * app.vel.factor) / (pipe.area * discharge.coef))^2
}


