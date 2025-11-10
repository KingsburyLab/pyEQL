# Partial pressures of reactive gases in atmosphere
ATMOSPHERE = {"CO2": -3.5, "O2": -0.6778}

# The amount, in moles, of an assemblage of pure phases that can react
# reversibly with the aqueous phase.

# See:
#   https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/html/final-38.html
# ".. This number of moles defines the maximum amount of the mineral or gas
#  that can dissolve. It may be possible to dissolve the entire amount without
#  reaching the target saturation index, in which case the solution will have a
#  smaller saturation index for this phase than the target saturation index.".

# We set it to 100 moles to approximate an "infinite" supply mineral/gas.
EQUILIBRIUM_PHASE_AMOUNT = 100
