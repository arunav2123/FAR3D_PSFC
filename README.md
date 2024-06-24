# FAR3D code
The aim of FAR3d code is to offer a fast tool to study the linear MHD and AE
stability of nuclear fusion devices. 
1. FAR3d code solves the linear reduced resistive
MHD equations [1,2,3] adding the Landau damping/growth (wave-particle
resonance effects) [4,5] and geodesic acoustic waves (parallel momentum
response of the thermal plasma) [6].

2. FAR3d includes the moment equations of the
energetic ion density and parallel velocity allowing treatment of linear wave-ion
resonances.

3. The model requires Landau closure relations that can be calibrated by more
complete kinetic models [6]. The simulations are based on an equilibria
calculated by the VMEC code [7], but future version will include other equilibria
options as SIESTA [8] or EFIT [9].

4. From the full MHD equations we can derive a reduced set of three equations that
follow the evolution of the pressure ‘p’, the stream function proportional to the
electrostatic potential ‘Φ’ and the perturbation of the poloidal flux ‘ψ’ retaining
the toroidal angle variation ‘ζ’ (exact three-dimensional equilibrium [1]).

5. The equations for the background or thermal plasma evolution are complemented by
moments of the energetic ion equations to add the effect of the wave-particle
interaction [10], namely the energetic particle density 'nf' and velocity moments
parallel to the magnetic field lines 'v||f'. 

6. The coefficients of the closure relation are selected to match a two-pole approximation of the plasma dispersion function.

## VMEC Setup

```
git submodule update --init --recursive
```
