# icymoons_thermal_stress
**This repository contains codes to reproduce the results in:
Rudolph, M.L., Manga, M., Walker, M., and Rhoden, A.R. (in review). Cooling crusts create concomittant cryovolcanic conduits.**

Thermal and stress evolution of icy moons, with application to Europa and Enceladus.
Requirements:
- MATLAB 2020a or later due to the use of tiled layouts.
- Fabio Crameri's scientific colormaps
- Parallel computing toolbox

Description of files:
- The main driver program that sweeps parameter space is ```sweep_parameter_space.m```
- The file ```thickening_ice_shell.m``` solves for the stresses in an initially stress-free thickening ice shell with no heat input.
- Helper subroutines are located in the subdirectory ```core/```
- Files in the ```benchmarks/``` subdirectory can be used to reproduce the results in Nimmo (2004).

