Started with v10_z_A6-atmos-addsbackinRerad but made major modifications so that it works for thin lags

Gets rid of extraneous, old parts of code too.

This version combines v11 with Free/Forced Convection basic thermal model to calculate the flat, regional values first
and then calculates and saves surface sublimation at each orbital timestep. Makes sure to reorder arrays to be ordered
in diurnal loops for the surface sublimation calculations.

This version calls several of the scripts to now be in functions, which are supposed to run faster, and will also
make sure nothing is inadvertantly carried over between timesteps.

It also fixes some issues in calcOrbitalParams where there was the possibility of getting imaginary numbers
in the solar azimuth angle which could then propogate into temperatures and all other calculations. It adds
checks to cosi and cos of the azimuth angle to make sure they stay between [-1,1].

As setup now with output, V12 is best used for just calculating free and forced sublimation values throughout
orbital histories. V13 will have a more comprehensive output that compares sublimation through a diffusive lag
to the outputs from runs of this V12 code (or maybe I'll decide to have it calculate the surface sublimation
values so the code is more comprehensive, even though it will make it slower).

New version of sublCalculator and setLayerProperties....

Updates to prevent trying to set really thin pore-filling ice layers. New parameter for thinnest layer. Only uses
for lag thickness though not pore filling because otherwise it has issues with calculating Inf or Nan vapor densities
and then can't interpolate to get z_eq.

May 8, 2018: Added output to tempCalculator (turned it into v13) that uses z_ei_index to print
the temperature convergence at the ice table. Updated calls to tempCalculator in TM_v13 and
findEquilibriumDepth

May 14, 2018: setLayerProperties now uses dry lag skin depth even when pore filling ice all the way up to
surface because otherwise often get layers thicker than the lag/pore filling ice layer(s)

May 14, 2018: Also new version of findEquilbriumDepth (version 13) that deals with cases where
pore filling ice is stable within ~1 mm of the surface for a dt=500 seconds. Ideally would use a smaller
dt (~100 seconds) but code already takes a while to run, and we know pore filling ice should be
stable so now if there are -Inf in the rho_v_mean's (due to temperature instabilities from small layers)
then sets z_pf = s.thinnestLayer and either moves on if groundicecounter > 1 or tries recalculating
using s.thinnestLayer, and it often settles on a z_pf around 1 mm (when s.thinnestLayer = 1 mm). 
Doesn't affect the ice sublimation results since know pore filling ice needs to be stable in this case,
so this gives it a way to move on somewhat quickly without needing to keep reducing the dt

May 14, 2018: Also updates output to include the temperature convergence at the z_ei
...This version as of May 14 is what will be v13_a on the Cluster

As of 10/30/2018 updates for treatment of atmospheric scattering and attenuation
