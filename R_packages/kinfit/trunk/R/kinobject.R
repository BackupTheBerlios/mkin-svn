# $Id$

# Copyright (C) 2008-2010 Johannes Ranke
# Contact: mkin-devel@lists.berlios.de

# This file is part of the R package kinfit

# kinfit is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <http://www.gnu.org/licenses/>

kinobject <- function(parent, type, system, 
        layers = NA, sampling_times = NA)
{
        kinobject <- list(parent = parent, 
                type = type, system = system)
        if (!is.na(layers[1])) kinobject$layers = layers
        if (!is.na(sampling_times[1])) {
                kinobject$sampling_times = layers
        }
        return(kinobject)
}
