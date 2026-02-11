# Partial pressures of reactive gases in atmosphere
# Obtained from https://www.noaa.gov/jetstream/atmosphere on 12/15/25
# (N2=78.084%, O2=20.946%, CO2=0.042%)
# "N2": -0.1074379 is removed from ATMOSPHERE as N2 is inert, and want to avoid reacting with O2 to form NO3-.
# See related discussion on the PHREEQC forum:https://phreeqcusers.org/index.php?topic=2371.0
ATMOSPHERE = {"CO2": -3.3767507, "O2": -0.6788989}

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
