# FAR3D code
The aim of FAR3d code is to offer a fast tool to study the linear MHD and AE
stability of nuclear fusion devices. 
1. FAR3d code solves the linear reduced resistive
MHD equations adding the Landau damping/growth (wave-particle
resonance effects) and geodesic acoustic waves (parallel momentum
response of the thermal plasma).

2. FAR3d includes the moment equations of the
energetic ion density and parallel velocity allowing treatment of linear wave-ion
resonances.

3. The model requires Landau closure relations that can be calibrated by more
complete kinetic models. The simulations are based on an equilibria
calculated by the VMEC code, but future version will include other equilibria
options as SIESTA or EFIT.

4. From the full MHD equations we can derive a reduced set of three equations that
follow the evolution of the pressure ‘p’, the stream function proportional to the
electrostatic potential ‘Φ’ and the perturbation of the poloidal flux ‘ψ’ retaining
the toroidal angle variation ‘ζ’ (exact three-dimensional equilibrium ).

5. The equations for the background or thermal plasma evolution are complemented by
moments of the energetic ion equations to add the effect of the wave-particle
interaction , namely the energetic particle density 'nf' and velocity moments
parallel to the magnetic field lines 'v||f'. 

6. The coefficients of the closure relation are selected to match a two-pole approximation of the plasma dispersion function.

## Booz_xform setup and compilation

```booz_xform``` is a package for computing Boozer coordinates in toroidal magnetohydrodynamic equilibria, including both stellarators and tokamaks.  ```booz_xform``` requires NetCDF, a C++ compiler, and python3. For many users, the code can be installed by executing
```
pip install booz_xform
```
## VMEC Setup and compilation
1. Simply use git submodule to get ```VMEC``` in ```FAR3D_PSFC``` folder 
```
git submodule update --init --recursive
```
2. ```cd VMEC2000/``` folder, and you could simply run from 
```
pip install .
```
 Sometimes you may get error associated with f90wrap installation. In that case, install numpy before running the above command:
```
pip install numpy 
pip install .
```
3. In this approach, we directly run the setup.py script present in the VMEC root folder. Ensure that the following modules are installed: cmake, scikit-build, ninja, numpy, and f90wrap. These packages could be installed by running the following commands:
```
pip install numpy
pip install cmake scikit-build ninja f90wrap
```
Sometimes you may have to install numpy before trying to install f90wrap.

4. Then, run
```
python setup.py build_ext
python setup.py install
```
5. in succession. At this point, you should be able to import the VMEC module in python. To test this, you can try the following command from the shell:
```
python -c "import vmec; print('success')"
```
