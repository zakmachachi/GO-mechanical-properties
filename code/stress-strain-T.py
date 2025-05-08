import ase
from ase.io import read, write
import numpy as np
from mace.calculators import MACECalculator
import re
from ase import units
from ase.md.bussi import Bussi
from ase.optimize import LBFGS
from ase.filters import FrechetCellFilter
from ase.md import MDLogger
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, ZeroRotation, Stationary
from mace.calculators import MACECalculator 

# z axis is the zigzag direction whilst x axis is the armchair direction 

calculator = MACECalculator(
    model_paths="MACE_model_swa.model", # Replace with path to GO-MACE-23
    device="cuda",
    default_dtype="float32",
    enable_cueq=True
    )
def run(config:ase.Atoms,  simulation_name : str, direction = 'z'): 
    config.set_pbc([True,True,True]) 
    
    match = re.search(r'_(\d+)K', simulation_name)
    if match:
        T = int(match.group(1))
    else: T = 300
    # by defualt T is set to be 300K, change here if needed

    config.calc = calculator
    T_init = T # K
    Eqr_arg = {
    'taut' :100 * 0.5 * units.fs, # ps
    'timestep' :0.5 * units.fs , # fs time step
    }

    strain_tensor = np.array([
    [1.0, 0, 0],
    [0, 1.0, 0],
    [0, 0, 1.0]])

    opt = LBFGS(FrechetCellFilter(config,constant_volume=True))
    opt.run(fmax=0.0001,  steps=10000)
    stress_GPa = config.get_stress() * 160.2176621
    config.info['stress_GPa'] = stress_GPa
    config.info['strain'] = 0  
    write(f"result_loading_{direction}_{simulation_name}.xyz", config, format='extxyz')
    cell = config.get_cell()

    # Set the momenta corresponding to T=300 K
    MaxwellBoltzmannDistribution(config, temperature_K=T_init)
    Stationary(config)
    ZeroRotation(config)

    # run NPT to relax the GO 
    dyn = Bussi(config, temperature_K=T_init, **Eqr_arg)
    def write_frame():
        dyn.atoms.write(f"result_loading_{direction}_{simulation_name}.xyz", append=True)

    dyn.attach(write_frame, interval=100)
    dyn.attach(
    MDLogger(dyn, config, f"result_loading_{direction}_{simulation_name}.txt", stress=True, peratom=True, mode="a"), interval=100
    )
    dyn.run(10000)

    if direction =='z': 
        dim = 2
    elif direction =='x': 
        dim = 0
    for i in np.linspace(0.01, 0.51, 1000):
        eqr_config = read(f"result_loading_{direction}_{simulation_name}.xyz")
        eqr_config.set_calculator(calculator)

        strain_tensor[dim,dim] = 1 + i # step wise increment of the strain
        eqr_config.set_cell(np.dot(cell, strain_tensor), scale_atoms=True)

        load = Bussi(eqr_config, temperature_K=T_init, **Eqr_arg)
        def write_frame_stress():
            load.atoms.write(f"result_loading_{direction}_{simulation_name}.xyz", append=True)
        load.attach(write_frame_stress, interval=100)
        load.attach(
        MDLogger(load, eqr_config, f"result_loading_{direction}_{simulation_name}.txt", stress=True, peratom=True, mode="a"), interval=100
        )
        load.run(2000) 
        stress_GPa = eqr_config.get_stress() * 160.2176621
        eqr_config.info['stress_GPa'] = stress_GPa
        eqr_config.info['strain'] = i 
    return 

go = read('GO.xyz')

run(go, 'go_300K','x')
run(go, 'go_300K','z')
