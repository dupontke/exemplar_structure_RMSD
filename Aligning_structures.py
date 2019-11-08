#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:

# ----------------------------------------
# PREAMBLE:

import MDAnalysis
from MDAnalysis.analysis.align import *
import sys
from distance_functions import *

# ----------------------------------------
# VARIABLE DECLARATION

selection = 'all'

# ZIKV 5gjb
reference_pdb = '/home/rbdavid/Projects/Molecular_Machines/Helicase_ZIKV/original_crystal_structures/5gjb.pdb'
reference_alignment_landmark = 'bynum 3480 3499 3522 3544 3564'
# structure_list
structure_list = [
                ['/home/rbdavid/Projects/Molecular_Machines/Helicase_DENV/original_crystal_structures/2JLV_chainA.pdb','bynum 3623 3642 3665 3687 3707','2JLV_chainA_above_aligned_5gjb.pdb'],
                ['/home/rbdavid/Projects/Molecular_Machines/Helicase_DENV/original_crystal_structures/2JLV_chainA.pdb','bynum 3645 3665 3687 3707 3727','2JLV_chainA_below_aligned_5gjb.pdb'],
                ['/home/rbdavid/Projects/Molecular_Machines/Helicase_ZIKV/original_crystal_structures/5gjb.pdb','bynum 3480 3499 3522 3544 3564','5GJB_aligned_5gjb.pdb'],
                ['/home/rbdavid/Projects/Molecular_Machines/Helicase_ZIKV/original_crystal_structures/5mfx.pdb','bynum 3362 3381 3404 3426 3446','5MFX_aligned_5gjb.pdb'],
                ['/home/rbdavid/Projects/Molecular_Machines/Helicase_HCV/denv_hcv_alignment/1a1v_aligned_to_2jlv_A.pdb','bynum 3908 3925 3945 3965 3985','1A1V_aligned_5gjb.pdb'],
                ['/home/rbdavid/Projects/Molecular_Machines/Helicase_HCV/denv_hcv_alignment/2f55_chainA_aligned_to_2jlv_A.pdb','bynum 99 115 134 153 172','2F55_chainA_aligned_5gjb.pdb'],
                ['/home/rbdavid/Projects/Molecular_Machines/Helicase_HCV/denv_hcv_alignment/3kqh_chainA_aligned_to_2jlv_A.pdb','bynum 3297 3315 3336 3357 3378','2KQH_chainA_aligned_5gjb.pdb'],
                ['/home/rbdavid/Projects/Molecular_Machines/Helicase_HCV/denv_hcv_alignment/3kqh_chainB_aligned_to_2jlv_A.pdb','bynum 3297 3315 3336 3357 3378','2KQH_chainB_aligned_5gjb.pdb'],
                ['/home/rbdavid/Projects/Molecular_Machines/Helicase_HCV/denv_hcv_alignment/3kqk_chainA_aligned_to_2jlv_A.pdb','bynum 3296 3313 3333 3353 3373','2KQK_chainA_aligned_5gjb.pdb']]

# ----------------------------------------
# FUNCTIONS:

# ----------------------------------------
# MAIN:
# ----------------------------------------
# PREP THE REFERENCE UNIVERSE
ref = MDAnalysis.Universe(reference_pdb)
ref_align = ref.select_atoms(reference_alignment_landmark)
print ref_align, ref_align.n_atoms
ref_sel = ref.select_atoms(selection)
ref_sel.translate(-ref_align.center_of_mass())
pos0 = ref_align.positions

ref_sel.write('reference_structure.pdb')

print RMSD(ref_align.positions,ref_align.positions,ref_align.n_atoms)

for i in structure_list:
        # ----------------------------------------
        # PREP THE UNIVERSE TO BE ALIGNED TO THE REFERENCE
        u = MDAnalysis.Universe(i[0])
        u_align = u.select_atoms(i[1])
        u_sel = u.select_atoms(selection)
        
        print i[2], u_align, u_align.n_atoms
        
        u_sel.translate(-u_align.center_of_mass())
        pos1 = u_align.positions
        R, rmsd = rotation_matrix(pos1,pos0)
        u_sel.rotate(R)
        u_sel.write(i[2])
        print RMSD(u_align.positions,ref_align.positions,u_align.n_atoms)

