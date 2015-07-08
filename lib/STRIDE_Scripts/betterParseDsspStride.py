# -*- coding: utf-8 -*-

import os
import numpy as np

from prody.atomic import ATOMIC_FIELDS
from prody.atomic import AtomGroup
from prody.utilities import gunzip, which, PLATFORM

def parseDSSP_better(dssp, ag, parseall=True, biomol=False):
    """Parse DSSP data from file *dssp* into :class:`~.AtomGroup` instance
    *ag*.  DSSP output file must be in the new format used from July 1995
    and onwards.  When *dssp* file is parsed, following attributes are added
    to *ag*:

    * *dssp_resnum*: DSSP's sequential residue number, starting at the first
      residue actually in the data set and including chain breaks; this number
      is used to refer to residues throughout.

    * *dssp_acc*: number of water molecules in contact with this residue \*10.
      or residue water exposed surface in Angstrom^2.

    * *dssp_kappa*: virtual bond angle (bend angle) defined by the three Cα
      atoms of residues I-2,I,I+2.  Used to define bend (structure code 'S').

    * *dssp_alpha*: virtual torsion angle (dihedral angle) defined by the four
      Cα atoms of residues I-1,I,I+1,I+2.Used to define chirality (structure
      code '+' or '-').

    * *dssp_phi* and *dssp_psi*: IUPAC peptide backbone torsion angles

    The following attributes are parsed when ``parseall=True`` is passed:

    * *dssp_bp1*, *dssp_bp2*, and *dssp_sheet_label*: residue number of first
      and second bridge partner followed by one letter sheet label

    * *dssp_tco*: cosine of angle between C=O of residue I and C=O of residue
      I-1.  For α-helices, TCO is near +1, for β-sheets TCO is near -1.  Not
      used for structure definition.

    * *dssp_NH_O_1_index*, *dssp_NH_O_1_energy*, etc.: hydrogen bonds; e.g.
      -3,-1.4 means: if this residue is residue i then N-H of I is h-bonded to
      C=O of I-3 with an electrostatic H-bond energy of -1.4 kcal/mol.  There
      are two columns for each type of H-bond, to allow for bifurcated H-bonds.

    See http://swift.cmbi.ru.nl/gv/dssp/DSSP_3.html for details."""

    if not os.path.isfile(dssp):
        raise IOError('{0} is not a valid file path'.format(dssp))
    if not isinstance(ag, AtomGroup):
        raise TypeError('ag argument must be an AtomGroup instance')

    dssp_file = open(dssp, 'r')
    dssp_output = dssp_file.readlines()

    n_atoms = ag.numAtoms()
    NUMBER = np.zeros(n_atoms, int)
    SHEETLABEL = np.zeros(n_atoms, np.array(['a']).dtype.char + '1')
    ACC = np.zeros(n_atoms, float)
    KAPPA = np.zeros(n_atoms, float)
    ALPHA = np.zeros(n_atoms, float)
    PHI = np.zeros(n_atoms, float)
    PSI = np.zeros(n_atoms, float)
    SECSTR = np.zeros(n_atoms, dtype=ATOMIC_FIELDS['secondary'].dtype)

    if parseall:
        BP1 = np.zeros(n_atoms, int)
        BP2 = np.zeros(n_atoms, int)
        NH_O_1 = np.zeros(n_atoms, int)
        NH_O_1_nrg = np.zeros(n_atoms, float)
        O_HN_1 = np.zeros(n_atoms, int)
        O_HN_1_nrg = np.zeros(n_atoms, float)
        NH_O_2 = np.zeros(n_atoms, int)
        NH_O_2_nrg = np.zeros(n_atoms, float)
        O_HN_2 = np.zeros(n_atoms, int)
        O_HN_2_nrg = np.zeros(n_atoms, float)
        TCO = np.zeros(n_atoms, float)

    segnm = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'*20)
    old_res_num = 0

    for i in range(len(dssp_output)):
        if dssp_output[i].startswith('  #  RESIDUE'):
            start = i

    start+=1
    chain_tracker = dssp_output[start][11]

    for i in range(start, len(dssp_output)):
        line = dssp_output[i]

        if line[13] == '!':
            continue

        new_res_num = int(line[5:10])
        new_chain_id = line[11]

        if biomol:
            if int(line[:5]) == 1:
                segidx = 0
            elif new_chain_id == chain_tracker and new_res_num < old_res_num:
                segidx += 1
            res = ag[(segnm[segidx], new_chain_id, int(line[5:10]), line[10].strip())]
        else:
            res = ag[(line[11], int(line[5:10]), line[10].strip())]

        if res is None:
            continue

        indices = res.getIndices()
        SECSTR[indices] = line[16].strip()
        NUMBER[indices] = int(line[:5])
        SHEETLABEL[indices] = line[33].strip()
        ACC[indices] = int(line[35:38])
        KAPPA[indices] = float(line[91:97])
        ALPHA[indices] = float(line[97:103])
        PHI[indices] = float(line[103:109])
        PSI[indices] = float(line[109:115])

        if parseall:
            BP1[indices] = int(line[25:29])
            BP2[indices] = int(line[29:33])
            NH_O_1[indices] = int(line[38:45])
            NH_O_1_nrg[indices] = float(line[46:50])
            O_HN_1[indices] = int(line[50:56])
            O_HN_1_nrg[indices] = float(line[57:61])
            NH_O_2[indices] = int(line[61:67])
            NH_O_2_nrg[indices] = float(line[68:72])
            O_HN_2[indices] = int(line[72:78])
            O_HN_2_nrg[indices] = float(line[79:83])
            TCO[indices] = float(line[85:91])

        old_res_num = int(line[5:10])

    ag.setData('dssp_secstr', SECSTR)
    ag.setData('dssp_resnum', NUMBER)
    ag.setData('dssp_sheet_label', SHEETLABEL)
    ag.setData('dssp_acc', ACC)
    ag.setData('dssp_kappa', KAPPA)
    ag.setData('dssp_alpha', ALPHA)
    ag.setData('dssp_phi', PHI)
    ag.setData('dssp_psi', PSI)

    if parseall:
        ag.setData('dssp_bp1', BP1)
        ag.setData('dssp_bp2', BP2)
        ag.setData('dssp_NH_O_1_index', NH_O_1)
        ag.setData('dssp_NH_O_1_energy', NH_O_1_nrg)
        ag.setData('dssp_O_NH_1_index', O_HN_1)
        ag.setData('dssp_O_NH_1_energy', O_HN_1_nrg)
        ag.setData('dssp_NH_O_2_index', NH_O_2)
        ag.setData('dssp_NH_O_2_energy', NH_O_2_nrg)
        ag.setData('dssp_O_NH_2_index', O_HN_2)
        ag.setData('dssp_O_NH_2_energy', O_HN_2_nrg)
        ag.setData('dssp_tco', TCO)

    return ag

def parseSTRIDE_better(stride, ag, biomol=False):
    """Parse STRIDE output from file *stride* into :class:`~.AtomGroup`
    instance *ag*.  STRIDE output file must be in the new format used
    from July 1995 and onwards.  When *stride* file is parsed, following
    attributes are added to *ag*:

    * *stride_resnum*: STRIDE's sequential residue number, starting at the
      first residue actually in the data set.

    * *stride_phi*, *stride_psi*: peptide backbone torsion angles phi and psi

    * *stride_area*: residue solvent accessible area"""

    if not os.path.isfile(stride):
        raise IOError('{0} is not a valid file path'.format(stride))
    if not isinstance(ag, AtomGroup):
        raise TypeError('ag argument must be an AtomGroup instance')

    stride_file = open(stride, 'r')
    stride_output = stride_file.readlines()

    n_atoms = ag.numAtoms()
    NUMBER = np.zeros(n_atoms, int)
    AREA = np.zeros(n_atoms, float)
    PHI = np.zeros(n_atoms, float)
    PSI = np.zeros(n_atoms, float)
    SECSTR = np.zeros(n_atoms, dtype=ATOMIC_FIELDS['secondary'].dtype)

    segnm = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ'*20)
    old_res_num = 0

    for i in range(len(stride_output)):
        line = stride_output[i]
        if line.startswith("ASG"):
            start = i
            break

    old_chain_id = stride_output[start][9]
    old_stride_num = int(stride_output[start][16:20])

    for i in range(start, len(stride_output)):
        line = stride_output[i]

        alt_loc = ''
        try:
            new_res_num = int(line[10:15])
        except:
            new_res_num = int(line[10:14])
            alt_loc = line[14]
        new_chain_id = line[9]
        new_stride_num = int(line[16:20])

        if biomol:
            if new_stride_num == 1:
                segidx = 0
            elif new_res_num < old_res_num and new_chain_id == old_chain_id:
                segidx += 1
            elif new_res_num < old_res_num and new_chain_id != old_chain_id and new_stride_num != old_stride_num:
                segidx = 0
            res = ag[(segnm[segidx], new_chain_id, new_res_num, alt_loc)]
        else:
            res = ag[(new_chain_id, new_res_num, alt_loc)]
        if res is None:
            continue
        indices = res.getIndices()
        SECSTR[indices] = line[24].strip()
        NUMBER[indices] = int(line[16:20])
        PHI[indices] = float(line[42:49])
        PSI[indices] = float(line[52:59])
        AREA[indices] = float(line[64:69])

        if alt_loc:
            old_res_num = int(line[10:14])
        else:
            old_res_num = int(line[10:15])
        old_chain_id = line[9]
        old_stride_num = int(line[16:20])

    ag.setData('stride_secstr', SECSTR)
    ag.setData('stride_resnum', NUMBER)
    ag.setData('stride_phi', PHI)
    ag.setData('stride_psi', PSI)
    ag.setData('stride_area', AREA)
    return ag