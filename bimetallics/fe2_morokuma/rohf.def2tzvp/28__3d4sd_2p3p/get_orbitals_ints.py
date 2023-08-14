import pyscf
import pyscf.tools

from orbitalpartitioning import *


molecule = """
 Fe 1.67785607 0.00052233 0.06475932
 O 0.00000000 0.00000000 -0.47099074
 Fe -1.67785607 -0.00052233 0.06475932
 Cl 1.87002704 -1.09796437 1.99091682
 Cl 2.93244917 -0.98210488 -1.47467288
 Cl 2.37160936 2.07954091 -0.50446591
 Cl -1.87002704 1.09796437 1.99091682
 Cl -2.93244917 0.98210488 -1.47467288
 Cl -2.37160936 -2.07954091 -0.50446591
 """

basis = "def2-tzvp"
pymol = pyscf.gto.Mole(
        atom    =   molecule,
        symmetry=   True,
        spin    =   10, # number of unpaired electrons
        charge  =   -2,
        basis   =   basis)


pymol.build()

print("symmetry: ",pymol.topgroup)
mf = pyscf.scf.ROHF(pymol)
mf.verbose = 4
mf.conv_tol = 1e-8
mf.conv_tol_grad = 1e-5
mf.chkfile = "scf.fchk"
mf.init_guess = "chkfile"
dm = np.load("rohf_dm.npy")
mf.run(dm)
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
for ao_idx,ao in enumerate(mf.mol.ao_labels(fmt=False)):
    if ao[0] == 0:
        if ao[2] in ("3d","4d","4s"):
            frag1.append(ao_idx)
            full.append(ao_idx)
    elif ao[0] == 1:
        if ao[2] in ("2p","3p"):
            frag2.append(ao_idx)
            full.append(ao_idx)
    elif ao[0] == 2:
        if ao[2] in ("3d","4d","4s"):
            frag3.append(ao_idx)
            full.append(ao_idx)


frags = [frag1, frag2, frag3]
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
