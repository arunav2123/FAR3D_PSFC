# This script is designed for compiling external profiles for FAR3d. This script should work for any CDF, but there are changes that might be needed depending on the variables in your CDF filed. These changes should be made in lines 19-27. Additionally, device specifics should be changed in lines 10, and 70-85. Finally, this script does not include collumns for alphas, but including these or any other profiles follows the same design.
# Leon Nichols, 11/15/2024, lonichol@mit.edu

from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import interp1d
import sys

input_file = 'KSTAR.CDF' 
R0, R_edge = 1.8, 2.255  # Major radius (on axis), and edge (R0+a) (in meters). This is needed for putting profiles in terms of the square root of normalized toroidal flux.
output_filename_template = 'KSTAR_t{time_slice}.txt'

# Loop over time_index values
for time_index in range(0, 27):  # Looping time slices. If you want a single index, simply change to ({index #}, {index #})
    print(f"Processing time_index {time_index}...")

    nc_file = Dataset(input_file, 'r')

    X = nc_file.variables['X'][time_index, :]  # Spatial coordinate. Check with your specific CDF file for what this is. This code assumes it is r/a.
    ne = nc_file.variables['NE'][time_index, :] * 1e6 * 1e-20 # Electron density in x10^20 N/m^3. First multiplying number is to change from cm^3 to m^3. Second is for 10^20.
    ni = nc_file.variables['NI'][time_index, :] * 1e6 * 1e-20 # Ion density
    nb = nc_file.variables['BDENS_D'][time_index, :] * 1e6 * 1e-20 # Beam ion density
    nx = nc_file.variables['NIMP'][time_index, :] * 1e6 * 1e-20 # Impurity density
    te = nc_file.variables['TE'][time_index, :] * 1e-3  # Electron temperature in keV. This was originally in eV.
    ti = nc_file.variables['TI'][time_index, :] * 1e-3  # Ion temperature
    tb = nc_file.variables['EBEAM_D'][time_index, :] * 1e-3  # Beam effective temperature
    q = nc_file.variables['Q'][time_index, :]  # Safety factor

    # Close the NetCDF file
    nc_file.close()

    # Transforming profiles to sqrt of normalized toroidal flux and interpolating to 101 points
    normalized_flux = np.sqrt(R0 / (R0 + (R_edge - R0) * X))
    new_flux = np.linspace(normalized_flux[0], normalized_flux[-1], 101)

    # Defining the fixed first (positional) column for the output (0, 0.01, ..., 1.00)
    rho_column = np.linspace(0, 1, 101)

    # Interpolate each profile
    interpolated_data = {'Rho(norml. sqrt. toroid. flux)': rho_column}
    for label, profile in zip(
            ['q', 'Beam Ion Density(10^20 m^-3)', 'Ion Density(10^20 m^-3)', 'Elec Density(10^20 m^-3)',
             'Impurity Density(10^20 m^-3)', 'Beam Ion Effective Temp(keV)', 'Ion Temp(keV)', 'Electron Temp(keV)'],
            [q, nb, ni, ne, nx, tb, ti, te]):
        interp_func = interp1d(normalized_flux, profile, kind='cubic')
        interpolated_data[label] = interp_func(new_flux)

    kB_keV_to_J = 1.60218e-16  # Boltzmann constant in keV/J

    # Computing derived quantities (Beam Pressure, Thermal Pressure, Equil. Pressure)
    beam_pressure = interpolated_data['Beam Ion Density(10^20 m^-3)'] * 1e20 * interpolated_data['Beam Ion Effective Temp(keV)'] * kB_keV_to_J * 1e-3
    thermal_pressure = interpolated_data['Ion Density(10^20 m^-3)'] * 1e20 * interpolated_data['Ion Temp(keV)'] * kB_keV_to_J * 1e-3 + interpolated_data['Elec Density(10^20 m^-3)'] * 1e20 * interpolated_data['Electron Temp(keV)'] * kB_keV_to_J * 1e-3
    equil_pressure = beam_pressure + 2 * thermal_pressure
    tor_rot = np.zeros_like(rho_column)  # Placeholder, assuming rotation values not provided
    pol_rot = np.zeros_like(rho_column)  # Placeholder, assuming rotation values not provided

    # Adding derived quantities to the data dictionary
    interpolated_data.update({
        'Beam Pressure(kPa)': beam_pressure,
        'Thermal Pressure(kPa)': thermal_pressure,
        'Equil.Pressure(kPa)': equil_pressure,
        'Tor Rot(km/s)': tor_rot,
        'Pol Rot(km/s)': pol_rot
    })

    # Writing the output to the specified text file format. This labels for chosen time slice.
    output_file = output_filename_template.format(time_slice=time_index)

    # All header information included here is for a specific device. This change be changed manually in the output file, but if you are doing multiple time slices it is better to change it here before looping.
    header = """\
PLASMA GEOMETRY
Vacuum Toroidal magnetic field at R=1.8m [Tesla] 
    3
Geometric Center Major radius [m]
    1.8
Minor radius [m]
    0.455
Avg. Elongation
    2
Avg. Top/Bottom Triangularity
    0.482
Main Contaminant Species
    19F
Main Ion Species mass/proton mass
    2
TRYING TO GET TO BETA(0)=0.072 , Rmax=1.8

Rho(norml. sqrt. toroid. flux), q, Beam Ion Density(10^20 m^-3), Ion Density(10^20 m^-3), Elec Density(10^20 m^-3), Impurity Density(10^20 m^-3), Beam Ion Effective Temp(keV), Ion Temp(keV), Electron Temp(keV), Beam Pressure(kPa), Thermal Pressure(kPa), Equil.Pressure(kPa), Tor Rot(km/s), Pol Rot(km/s)
"""
    # This is an alternate way to write the profiles. Here you can make the output file more compact, but slightly less organized.
    #with open(output_file, 'w') as f:
    #    # Write the fixed header
    #    f.write(header )#+ '\n')

    #    for i in range(len(rho_column)):
    #        row = [interpolated_data[key][i] for key in interpolated_data]
    #        f.write(', '.join(f"{value:.7g}" for value in row) + '\n')

    # Open the file for writing
    with open(output_file, 'w') as f:
        # Write the header
        f.write(header )

        # Writing the profile data with specific formatting.
        for i in range(len(rho_column)):
            row = []
            for j, key in enumerate(interpolated_data):
                value = interpolated_data[key][i]
                if j == 0 or j >= len(interpolated_data) - 2:  # First column and last two columns chopped to two decimal places. Change this if including poloidal and toroidal rotation.
                    row.append(f"{value:.2f}")
                else:  # All other columns
                    row.append(f"{value:.6f}")
            f.write(', '.join(row) + '\n')  # Join the formatted row and write it as a line


    print(f"Profiles written to {output_file}")
