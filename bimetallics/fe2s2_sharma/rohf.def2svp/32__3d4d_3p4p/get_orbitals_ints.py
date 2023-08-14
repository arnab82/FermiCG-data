import pyscf
import pyscf.tools

from orbitalpartitioning import *


molecule = """
Fe  5.48 1.15 -8.03
S   4.05 -0.61 -8.75
S   7.49 0.42 -9.04
Fe  6.04 -1.22 -9.63
S   5.47 1.25 -5.58
S   4.63 3.28 -8.77
S   5.75 -1.50 -12.05
S   6.86 -3.41 -8.86
C   5.51 4.45 -7.51
H   6.49 4.83 -7.92
H   4.87 5.33 -7.25
H   5.72 3.84 -6.59
C   3.60 1.70 -5.54
H   3.01 0.80 -5.82
H   3.28 2.06 -4.52
H   3.42 2.48 -6.31
C   5.21 -4.22 -9.46
H   5.10 -4.01 -10.55
H   5.21 -5.32 -9.26
H   4.37 -3.72 -8.93
C   7.63 -1.85 -12.24
H   7.90 -2.06 -13.31
H   8.20 -0.96 -11.86
H   7.89 -2.72 -11.59
 """

pymol = pyscf.gto.Mole(
        atom    =   molecule,
        symmetry=   True,
        spin    =   10, # number of unpaired electrons
        charge  =   -2)
#
#                   ROHF/sto-3g
#
pymol.basis = 'sto-3g'
pymol.build(False, False)
mf = pyscf.scf.ROHF(pymol).density_fit()
mf.verbose = 4
mf.conv_tol = 1e-8
mf.conv_tol_grad = 1e-5
mf.chkfile = "scf.fchk"
mf.run()

#
#                   ROHF/def2-svp
#
pymol.basis = 'def2-svp'
pymol.build(False, False)
mf = pyscf.scf.ROHF(pymol).density_fit().newton()
mf.verbose = 4
mf.conv_tol = 1e-8
mf.conv_tol_grad = 1e-5
mf.chkfile = "scf.fchk"
mf.init_guess = "chkfile"
mf.run()


pyscf.tools.molden.from_mo(mf.mol, "Csing.molden", mf.mo_coeff[:,mf.mo_occ==1])
##################################################################

#
#   Get data
F = mf.get_fock()
C = mf.mo_coeff
S = mf.get_ovlp()
Cdocc = mf.mo_coeff[:,mf.mo_occ==2]
Csing = mf.mo_coeff[:,mf.mo_occ==1]
Cvirt = mf.mo_coeff[:,mf.mo_occ==0]
ndocc = Cdocc.shape[1]
nsing = Csing.shape[1]
nvirt = Cvirt.shape[1]


# Find AO's corresponding to atoms
full = []
frag1 = []
frag2 = []
frag3 = []
frag4 = []
for ao_idx,ao in enumerate(mf.mol.ao_labels(fmt=False)):
    if ao[0] == 0:
        if ao[2] in ("3d","4d"):
            frag1.append(ao_idx)
            full.append(ao_idx)
    elif ao[0] == 1:
        if ao[2] in ("3p","4p"):
            frag2.append(ao_idx)
            full.append(ao_idx)
    elif ao[0] == 2:
        if ao[2] in ("3p","4p"):
            frag3.append(ao_idx)
            full.append(ao_idx)
    elif ao[0] == 3:
        if ao[2] in ("3d","4d"):
            frag4.append(ao_idx)
            full.append(ao_idx)


frags = [frag1, frag2, frag3, frag4]
print(frags)



#
#   Get full active space
(Oact, Sact, Vact), (Cenv, Cerr, _) = svd_subspace_partitioning_nonorth((Cdocc, Csing, Cvirt), full, S)
assert(Cerr.shape[1] == 0)


#
#   Split active space into fragments
init_fspace = []
clusters = []
Cfrags = []
orb_index = 1


for fi,f in enumerate(frags):
    print()
    print(" Fragment: ", f)
    (Of, Sf, Vf), (_, _, _) = svd_subspace_partitioning_nonorth((Oact, Sact, Vact), f, S)

    Cfrags.append(np.hstack((Of, Sf, Vf)))
    ndocc_f = Of.shape[1]
    init_fspace.append((ndocc_f+Sf.shape[1], ndocc_f))
    nmof = Of.shape[1] + Sf.shape[1] + Vf.shape[1]
    clusters.append(list(range(orb_index, orb_index+nmof)))
    orb_index += nmof


# Orthogonalize Fragment orbitals
Cfrags = sym_ortho(Cfrags, S)

# Pseudo canonicalize fragments
Cfrags = canonicalize(Cfrags, F)


Cact = np.hstack(Cfrags)

print(" sing vals of active space: ", np.linalg.svd(Cact.T @ S @ Cact)[1])
# Write Molden files for visualization
pyscf.tools.molden.from_mo(mf.mol, "Cact.molden", Cact)
pyscf.tools.molden.from_mo(mf.mol, "Cenv.molden", Cenv)
for i in range(len(frags)):
    di = Cfrags[i] @ Cfrags[i].T
    pyscf.tools.cubegen.density(mf.mol, 'fragden_{:02d}.cube'.format(i+1), di)


print(" init_fspace = ", init_fspace)
print(" clusters    = ", clusters)




#
#   Make integrals
#
d1_embed = 2 * Cenv @ Cenv.T

h0 = pyscf.gto.mole.energy_nuc(mf.mol)
h  = pyscf.scf.hf.get_hcore(mf.mol)
j, k = pyscf.scf.hf.get_jk(mf.mol, d1_embed, hermi=1)

h0 += np.trace(d1_embed @ ( h + .5*j - .25*k))

h = Cact.T @ h @ Cact;
j = Cact.T @ j @ Cact;
k = Cact.T @ k @ Cact;
nact = h.shape[0]

h2 = pyscf.ao2mo.kernel(pymol, Cact, aosym="s4", compact=False)
h2.shape = (nact, nact, nact, nact)
# The use of d1_embed only really makes sense if it has zero electrons in the
# active space. Let's warn the user if that's not true

S = pymol.intor("int1e_ovlp_sph")
n_act = np.trace(S @ d1_embed @ S @ Cact @ Cact.T)
if abs(n_act) > 1e-8 == False:
    print(n_act)
    error(" I found embedded electrons in the active space?!")

h1 = h + j - .5*k;

np.save("ints_h0", h0)
np.save("ints_h1", h1)
np.save("ints_h2", h2)
np.save("mo_coeffs", Cact)
np.save("overlap_mat", S)

Pa = mf.make_rdm1()[0]
Pb = mf.make_rdm1()[1]
np.save("Pa", Cact.T @ S @ Pa @ S @ Cact)
np.save("Pb", Cact.T @ S @ Pb @ S @ Cact)
