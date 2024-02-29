import numpy as np
 
total_widths=[40]
Apf_widths=[0]
 
 
chemical_potent=[1.0]
 
delta=[0.1]
 
Ts=np.geomspace(0.0001,0.02,35)
 
with open('jobs.sh','w') as f:
   for L in total_widths:
           for alpha in chemical_potent:
               for d in delta:
                   for T in Ts:
                       text='addqueue -c "15 min" -m 2 /usr/bin/python3 /mnt/users/dperkovic/analytic_helper_calculations/bogoliubov/run_bdg.py '+str(L)+' '+str(alpha)+' '+str(d)+' '+str(T)+'\n'
                       f.write(text)