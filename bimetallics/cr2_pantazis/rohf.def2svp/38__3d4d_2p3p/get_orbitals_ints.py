import pyscf
import pyscf.tools

from orbitalpartitioning import *


molecule = """
 Cr -1.32077675789660 0.00004567211811 -0.00007096467317
 Cr 1.32077324283384 0.00005158303881 -0.00007157859457
 O -0.00000104868858 -0.16583237580723 1.45467585745211
 O -0.00000489764937 1.34277292975137 -0.58371627199903
 O 0.00000067419445 -1.17682896850597 -0.87101026120097
 H 0.00001790208267 0.50127981586150 2.15993328303163
 H 0.00056213619050 1.61869136067886 -1.51448427772977
 H -0.00044104651121 -2.12078532888951 -0.64413286649204
 N -2.64979639176962 -1.44569378219198 0.71141991748241
 C -1.96723877683042 -2.53150336334491 1.49741316302721
 C -3.24368037350626 -2.03126564893063 -0.53390524471276
 C -3.70971767634131 -0.79482652965216 1.56203025179475
 H -1.50611109500905 -2.09494287212245 2.38227714666532
 H -1.20003332072637 -2.99746593446751 0.87727778775072
 H -2.69709429184401 -3.29417670409302 1.78680599374014
 H -2.45713087989536 -2.63082378918971 -1.00368039677010
 H -4.06474101618296 -2.71373529840676 -0.28162892858883
 H -4.62553838458745 -0.71938852428422 0.97193419906511
 H -3.94180643056544 -1.43233893352354 2.41956824499303
 N -2.64980256138524 1.33901556828090 0.89630056028452
 N -2.64980038336441 0.10677232402498 -1.60776835378432
 C -3.70971800240785 -0.95531484258297 -1.46940407035238
 C -3.24368562009935 0.55331550766570 2.02608271830195
 N 2.64980067794514 -1.44568193229971 0.71141868600201
 N 2.64979447720107 1.33902743110622 0.89629932793762
 N 2.64979663544791 0.10678418068707 -1.60776957870374
 C -1.96725035879470 2.56261235577474 1.44364406188149
 C -3.70972375060080 1.75022808457439 -0.09267425723559
 C -1.96724550946031 -0.03101115292470 -2.94110526383651
 C -3.24368848245842 1.47803937363856 -1.49222464129291
 H -4.62558984511764 -0.48221330653117 -1.10887816237189
 H -3.94175413156726 -1.37940439028017 -2.45019713433286
 H -2.45710023893015 0.44599637215694 2.78010888133082
 H -4.06480735461771 1.11280848002045 2.49108824409826
 C 1.96724725554052 -2.53149455336756 1.49741224378802
 C 3.24368569356644 -2.03125112915353 -0.53390676339270
 C 3.70971844014110 -0.79480993568548 1.56202852714936
 C 1.96723601210309 2.56262114404852 1.44364312933477
 C 3.24368056139230 0.55333001361486 2.02608120711149
 C 3.70971236083504 1.75024468234002 -0.09267599934625
 C 1.96724081992101 -0.03100234186851 -2.94110620541212
 C 3.24367767216167 1.47805387752771 -1.49222616709369
 C 3.70971807894818 -0.95529823271254 -1.46940580237183
 H -1.50375290990065 3.10925190551415 0.62354092184389
 H -1.20171313853662 2.25877752450171 2.15909815420654
 H -2.69767375679224 3.19578710394873 1.95706724009970
 H -4.62543371684035 1.20116961382526 0.13681738047320
 H -3.94200796666106 2.81155859771218 0.03059305603763
 H -1.50276864801509 -1.01409155231681 -3.00420649677802
H -1.20309418638945 0.74187206739252 -3.03572600131324
 H -2.69797575471462 0.09546145526857 -3.74614669451067
 H -2.45703697100368 2.18462404933411 -1.77634058940058
 H -4.06475178550869 1.60090730318433 -2.20935862370879
 H 1.50614862542804 -2.09494317731569 2.38229660706312
 H 1.20002265002298 -2.99741709798907 0.87727577268110
 H 2.69710252932602 -3.29417828809696 1.78677807515430
 H 2.45714267657494 -2.63082369046460 -1.00366709050502
 H 4.06476961090010 -2.71369418770397 -0.28163224842656
 H 4.62555565455524 -0.71940605121393 0.97195330268582
 H 3.94178307498783 -1.43230113185123 2.41958979056064
 H 1.50357479424781 3.10918222005415 0.62358103033223
 H 1.20182998798639 2.25879076334749 2.15924302589083
 H 2.69768680938841 3.19588061502689 1.95692445402491
 H 2.45708429727129 0.44597962062585 2.78009130737457
 H 4.06478059679567 1.11284184879546 2.49110060500510
 H 4.62541221159885 1.20116314525263 0.13680039348052
 H 3.94201224677884 2.81156990751106 0.03060321408914
 H 1.50267872785295 -1.01404489669261 -3.00416424100455
 H 1.20314783285142 0.74193593536788 -3.03579980829302
 H 2.69798268801952 0.09538099720936 -3.74615088499019
 H 2.45703083735398 2.18464798306617 -1.77634272668040
 H 4.06473833125488 1.60093248495915 -2.20936143690374
 H 4.62558786531673 -0.48219041570644 -1.10888361173207
 H 3.94174955515341 -1.37938388964014 -2.45020164669110
 """

pymol = pyscf.gto.Mole(
        atom    =   molecule,
        symmetry=   True,
        spin    =   6, # number of unpaired electrons
        charge  =   3)
#
#                   ROHF/sto-3g
#
pymol.basis = 'sto-3g'
pymol.build()
mf = pyscf.scf.ROHF(pymol).density_fit()
mf.verbose = 4
mf.conv_tol = 1e-8
mf.conv_tol_grad = 1e-5
mf.chkfile = "scf.fchk"
mf.run()

pyscf.tools.molden.from_mo(mf.mol, "Csing_sto3g.molden", mf.mo_coeff[:,mf.mo_occ==1])

#
#                   ROHF/def2-svp
#
pymol.basis = 'def2-svp'
pymol.build()
mf = pyscf.scf.ROHF(pymol).density_fit()
mf.verbose = 4
mf.conv_tol = 1e-8
mf.conv_tol_grad = 1e-5
mf.chkfile = "scf.fchk"
mf.init_guess = "chkfile"
mf.run()


pyscf.tools.molden.from_mo(mf.mol, "Csing.molden", mf.mo_coeff[:,mf.mo_occ==1])
##################################################################
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
frag5 = []
for ao_idx,ao in enumerate(mf.mol.ao_labels(fmt=False)):
    if ao[0] == 0:
        if ao[2] in ("3d","4d"):
            frag1.append(ao_idx)
            full.append(ao_idx)
    elif ao[0] == 1:
        if ao[2] in ("3d","4d"):
            frag2.append(ao_idx)
            full.append(ao_idx)
    elif ao[0] == 2:
        if ao[2] in ("2p","3p"):
            frag3.append(ao_idx)
            full.append(ao_idx)
    elif ao[0] == 3:
        if ao[2] in ("2p","3p"):
            frag4.append(ao_idx)
            full.append(ao_idx)
    elif ao[0] == 4:
        if ao[2] in ("2p","3p"):
            frag5.append(ao_idx)
            full.append(ao_idx)


frags = [frag1, frag2, frag3, frag4, frag5]
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
