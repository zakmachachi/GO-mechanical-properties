import ase
from ase.io import read, write
import numpy as np
from mace.calculators import MACECalculator
import re
from ase import units
from ase.optimize import LBFGS
from ase.filters import FrechetCellFilter
from ase.md import MDLogger
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from mace.calculators import MACECalculator 

# z axis is the zigzag direction whilst x axis is the armchair direction 

calculator = MACECalculator(
    model_paths="MACE_model_swa.model", # replace with path to GO-MACE-23
    device="cuda",
    default_dtype="float32",
    enable_cueq=False
    )
def run(config:ase.Atoms,  simulation_name : str, direction = 'z'): 
    config.set_pbc([True,True,True]) 

    config.calc = calculator

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

    if direction =='z': 
        dim = 2
    elif direction =='x': 
        dim = 0

    for i in np.linspace(0.01, 0.51, 1000):
        eqr_config = read(f"result_loading_{direction}_{simulation_name}.xyz", index=-1)
        eqr_config.set_calculator(calculator)
        strain_tensor[dim,dim] = 1 + i # step wise increment of the strain
        eqr_config.set_cell(np.dot(cell, strain_tensor), scale_atoms=True)
        load = LBFGS(eqr_config)
        load.run(fmax=0.001, steps=1000)
        stress_GPa = eqr_config.get_stress() * 160.2176621
        eqr_config.info['stress_GPa'] = stress_GPa
        eqr_config.info['strain'] = i 
        write(f"result_loading_{direction}_{simulation_name}.xyz", eqr_config, write_results=True, write_info=True, append=True)
    return 

go = read('rGO.xyz')

run(go, 'rGO_geom','x')
run(go, 'rGO_geom','z')
