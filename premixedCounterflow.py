#!/usr/bin/env python3
"""
Simulate two counter-flow jets of reactants shooting into each other. This
simulation differs from the similar premixed_counterflow_flame.py example as the
latter simulates a jet of reactants shooting into products.

Requires: cantera >= 2.5.0
"""

import cantera as ct
import numpy as np
import sys
import pandas as pd
import os
import time

print('\n*** Computation of premixed counter-flow twin flames  ***\n\n')

# Select the reaction mechanism
mec = 'chemicalMechanism/kee.xml'
gas = ct.Solution(mec)

Tmax_s = []
Sc_s = []
Sd_s = []
strain_s = []
massFlux_s = []

#Set input velocity
axial_velocity = np.linspace(1,5,10)

phi = 1.
folderName ='{:.2f}'.format(phi)

path = './counterFlowResults/CH4/' + folderName
if not os.path.isdir(path):
    os.makedirs(path)
    print("created folder : ", path)
else:
    print(path, " folder already exists.")

print('Path to Save: ' + path) 
time.sleep(6)     

for i in range(0,axial_velocity.size):
    # Create a CH4/Air premixed mixture with equivalence at room
    # temperature and pressure.
    fuel = 'CH4'
    gas.set_equivalence_ratio(phi, fuel, {'O2':1.0, 'N2':3.76})
    gas.TP = 300, ct.one_atm

    # Domain half-width of 2.5 cm, meaning the whole domain is 5 cm wide
    width = 0.025

    # Done with initial conditions
    # Compute the mass flux, as this is what the Flame object requires
    massFlux = gas.density * axial_velocity[i]  # units kg/m2/s
    
    # Create the flame object
    oppFlame = ct.CounterflowTwinPremixedFlame(gas, width=width)
    oppFlame.max_grid_points = 5e4

    # Uncomment the following line to use a Multi-component formulation. Default is
    # mixture-averaged
    #oppFlame.transport_model = 'Multi'
    #oppFlame.soret_enabled=True
    #oppFlame.transport_model = 'UnityLewis'
    oppFlame.transport_model = 'Mix'
         
    oppFlame.reactants.mdot = massFlux
    oppFlame.set_refine_criteria(ratio=2, slope=0.02, curve=0.02, prune=0.00)

    oppFlame.show_solution()
    oppFlame.solve(loglevel = 1, auto=True)
    T_max = np.max(oppFlame.T)

    if T_max < 500:
        print("\n** Flame extinction\ " )
        break
    
    print("Peak temperature: {0:.1f} K".format(T_max))
    print("Mass flux: {0:.4f} Kg/m2s".format(massFlux))

    list_species = ['CH4','O2','CO','CO2',\
                 'H2O','OH','CH2O','H2O2','HO2','HCO']

    #list_species = ['H2','O2','H2O','OH','H2O2','HO2']

    df = pd.DataFrame()
    df['x'] =  oppFlame.grid
    df['rho'] =  oppFlame.density
    df['T'] =  oppFlame.T
    df['velocity'] =  oppFlame.velocity
    for species in list_species:
      df[species] =  oppFlame.Y[gas.species_index(species),:]

    for species in list_species:
      df['wdot' + species] =  oppFlame.net_production_rates[gas.species_index(species),:]*gas.molecular_weights[gas.species_index(species)]  

    for species in list_species:
      df['diff' + species] =  oppFlame.mix_diff_coeffs_mass[gas.species_index(species),:]  

    df['alpha'] =  oppFlame.thermal_conductivity/(oppFlame.cp_mass*oppFlame.density)
    df['k'] =  oppFlame.thermal_conductivity
    df['Qdot'] =  abs(oppFlame.heat_release_rate)

    #df['Z_C'] =  oppFlame.elemental_mass_fraction('C')
    df['Z_O'] =  oppFlame.elemental_mass_fraction('O')
    df['Z_H'] =  oppFlame.elemental_mass_fraction('H')
    df['Z_N'] =  oppFlame.elemental_mass_fraction('N')
    
    fileName = '{:.3f}'.format( axial_velocity[i] )    
    df.to_csv(path + '/' + fileName, index = False)
     
