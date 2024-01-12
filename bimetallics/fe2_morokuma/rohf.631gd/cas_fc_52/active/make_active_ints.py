import numpy as np
import pyscf
import pyscf.tools


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

basis = "6-31g*"
pymol = pyscf.gto.Mole(
        atom    =   molecule,
        symmetry=   True,
        spin    =   10, # number of unpaired electrons
        charge  =   -2,
        basis   =   basis).build()


mf = pyscf.scf.ROHF(pymol)
mf.verbose = 4
mf.conv_tol = 1e-8
mf.conv_tol_grad = 1e-5
mf.chkfile = "scf.fchk"
mf.init_guess = "chkfile"
mf.run()

Ccmf = np.load("Ccmf.npy")
pyscf.tools.molden.from_mo(pymol, "Ccmf.molden", Ccmf)



Ccore = Ccmf[:,0:15]
Cact = Ccmf[:,15:67]
pyscf.tools.molden.from_mo(pymol, "Cact.molden", Cact)


#   Make integrals
#
d1_embed = 2 * Ccore @ Ccore.T

S = pymol.intor("int1e_ovlp_sph")
print(" Number of electrons in core: ", np.trace(d1_embed@S))

h0 = pyscf.gto.mole.energy_nuc(mf.mol)
h  = pyscf.scf.hf.get_hcore(mf.mol)
j, k = pyscf.scf.hf.get_jk(mf.mol, d1_embed, hermi=1)

print(h0)
h0 += np.trace(d1_embed @ ( h + .5*j - .25*k))
print(h0)

h = Cact.T @ h @ Cact;
j = Cact.T @ j @ Cact;
k = Cact.T @ k @ Cact;
nact = h.shape[0]

h2 = pyscf.ao2mo.kernel(pymol, Cact, aosym="s4", compact=False)
h2.shape = (nact, nact, nact, nact)
# The use of d1_embed only really makes sense if it has zero electrons in the
# active space. Let's warn the user if that's not true

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
