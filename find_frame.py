#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:

# ----------------------------------------
# PREAMBLE:

import numpy as np
import MDAnalysis
from MDAnalysis.analysis.align import *
import sys

data_file = sys.argv[1]        # the RMSD data file with the avg structure as ref 001.100.wt_ssrna_atp_1.exemplar_structure_rmsd.dat
index = int(sys.argv[2])       # the column in the data file you want to use i.e. 2 for full protein 
pdb_file = sys.argv[3]         # ../traj/truncated.pdb
ref_pdb_file = sys.argv[4]     # the average structure pdb associated with the RMSD data file
system = sys.argv[5]           # the system of interest i.e. wt_ssrna_atp_1
start = int(sys.argv[6])       # the start frame i.e. if you want to look at 071.100, then the start frame would be 355000
end = int(sys.argv[7])         # the end frame i.e. if you want to look at 071.100, then the end frame would be 500000
start_traj = int(sys.argv[8])  # if you are looking at 071.100, the start traj would be 071
end_traj = int(sys.argv[9])    # if you are looking at 071.100, the end traj would be 100

#alignment = 'protein and name CA and (resid 68:73 91:96 117:120 135:138 171:175 197:199 211:214)'
alignment = 'protein and name CA and (resid 20:25 50:55 73:75 90:94 112:116 142:147 165:169 190:194 214:218 236:240 253:258 303:307)'
def summary():
    with open('%03d.%03d.%s.find_frame.summary' %(start_traj,end_traj,system),'w') as f:
        f.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
        f.write('To recreate this analysis, run this line:\n')
        for i in range(len(sys.argv)):
            f.write('%s ' %(sys.argv[i]))
        f.write('\n\n')
        f.write('output is written to:\n')
        f.write('%03d.%03d.%s.frame_%s.pdb\n' %(start_traj,end_traj,system,index_min))

data = np.loadtxt(data_file)
index_min = np.argmin(data[start:end,index]) + start
min_value = np.amin(data[start:end,index])

print index_min, min_value

traj_num = int(index_min/5000)+1	# trajectory index starts at 1; 5000 frames per trajectory...

frame_num = index_min%5000		# frame indexing (for MDAnalysis) starts at 0

print traj_num, frame_num

u = MDAnalysis.Universe(pdb_file, '../traj/Truncated/production.%s/production.%s.dcd' %(traj_num, traj_num))
ref = MDAnalysis.Universe(ref_pdb_file)

u_system = u.select_atoms('all')
ref_system = ref.select_atoms('all')
u.trajectory[frame_num]

align_u = u.select_atoms(alignment)
align_ref = ref.select_atoms(alignment)

u_system.translate(-align_u.center_of_mass())
ref_system.translate(-align_ref.center_of_mass())

R, rmsd = rotation_matrix(align_u.positions,align_ref.positions)
u_system.rotate(R)

u_system.write('%03d.%03d.%s.frame_%s.pdb' %(start_traj,end_traj,system,index_min))

summary()
