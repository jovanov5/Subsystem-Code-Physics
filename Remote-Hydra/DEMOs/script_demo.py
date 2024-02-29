#!/usr/bin/env python3

import matplotlib
matplotlib.use('Agg')
import numpy as np
import sys
import matplotlib.pyplot as plt
import store_bogoliubov
import bogoliubov_space as bogoliubov_space_final

from full_system_diagonalize import Full_system



#CODE store_data()

total_width= float(sys.argv[1])
chemical_potential=float(sys.argv[2])
delta=float(sys.argv[3])
T=float(sys.argv[4])
number=700
num_basis=110
N=200
L_y=150



number_electron=0
solver=Full_system(num_basis,N,L_y,total_width,number,delta,chemical_potential,T=T,number_electrons=number_electron)
filling=chemical_potential
if solver.num_electron_initial!=0:
 solver.bisection_chemical_potential()
solver.run_full_system()
store_data(solver,filling)