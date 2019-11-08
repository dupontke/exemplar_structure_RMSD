# exemplar_structure_RMSD

First run the RMSD calculation with the Average structure as the reference structure:
./rmsd.ref.py rmsd.config

Then run the find_frame.py to obtain the exemplar structure:
./find_frame.py 021.100.wt_ssrna_atp_1.exemplar_structure_rmsd.dat 2 ../traj/truncated.pdb ../Average_structure_calcs/021.100.average_structure.pdb wt_ssrna_atp_1 105000 500000 021 100