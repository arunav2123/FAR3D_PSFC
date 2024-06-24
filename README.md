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

## VMEC Setup and compilation
1. Get VMEC in FAR3D_PSFC folder
```
git submodule update --init --recursive
```
2. Go to VMEC fodler, and you could simply run from 
```
pip install .
```
 Sometimes you may get error associated with f90wrap installation. In that case, install numpy before running the above command:
```
pip install numpy 
pip install .
```
4. In this approach, we directly run the setup.py script present in the VMEC root folder. Ensure that the following modules are installed: cmake, scikit-build, ninja, numpy, and f90wrap. These packages could be installed by running the following commands:
```
pip install numpy
pip install cmake scikit-build ninja f90wrap
```
Sometimes you may have to install numpy before trying to install f90wrap.

6. Then, run
```
python setup.py build_ext
python setup.py install
```
6. in succession. At this point, you should be able to import the vmec module in python. To test this, you can try the following command from the shell:

```python -c "import vmec; print('success')" ```
