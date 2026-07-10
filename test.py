#from Qlapoti_main.sage_implementation.applications.PRISM.quaternions.ideals import *
#from Qlapoti_main.sage_implementation.applications.PRISM.precomputations_PRSIM import *
#from sqisign_prism_v2_main import params
import sage.all as sage
from POKE_PKE_main.montgomery_isogenies.kummer_line import *
#from POKE_PKE_main.montgomery_isogenies.kummer_isogeny import *
from POKE_PKE_main.montgomery_isogenies.isogenies_x_only import isogeny_from_scalar_x_only
from sqisign_prism_v2_main.quaternions import *
from sqisign_prism_v2_main.qlapoti import IdealToIsogeny
#from sqisign_prism_v2_main.theta.theta_isogenies.product_isogeny import *
from sqisign_prism_v2_main.theta.theta_structures.couple_point import *
from sqisign_prism_v2_main.hd import Dim2Iso
from sqisign_prism_v2_main.ec import *

p=2**512 * 3**5-1

Fp = GF(p)
Fp2, Fp2_i = GF(p**2, name="i", modulus=[1, 0, 1]).objgen()
E0 = EllipticCurve(Fp2, [1, 0])

E0.set_order((p + 1) ** 2)
print("E0=",E0)
P8,Q8=TorsionBasis(E0)
params.set_spider_params(1)
#print("P8=",P8)
print("E0",params.E0)
#print("P8",P8)
#print("Q8",Q8)