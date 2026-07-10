from Qlapoti_main.sage_implementation.ideal_to_isogeny.quaternion_helpers.helpers import *
from Qlapoti_main.sage_implementation.ideal_to_isogeny.theta_isogenies.product_isogeny import *
from Qlapoti_main.sage_implementation.applications.PRISM.precomputations_PRSIM import *
from Qlapoti_main.sage_implementation.ideal_to_isogeny.ideal_to_isogeny_qlapoti import *
from Qlapoti_main.sage_implementation.applications.PRISM.params import *
from Qlapoti_main.sage_implementation.applications.PRISM.PRISM_sign import *
#from POKE_PKE_main.utilities.supersingular import *
#from POKE_PKE_main.montgomery_isogenies.kummer_isogeny import *

p=5*2**248 - 1
e=248

B = QuaternionAlgebra(-1, -p) 
i, j, k = B.gens()
O0 = B.quaternion_order([1, i, (i + j) / 2, (1 + k) / 2])
F = GF(p**2, name="i", modulus=[1, 0, 1])
E0 = EllipticCurve(F, [1, 0])
E0.set_order((p + 1) ** 2)
#print(E0)

#p2=p**2
#N=p2.next_prime()
#print("N= ",N)
N = next_prime(randint(p, p**2))

O_N = O_mod_N(O0, N)
J=O_N.random_ideal()
mu1, mu2, d1, d2, θ, N1, I, betaij, e =qlapoti(J, e, odd_norm_output=True)

print("mu1=", mu1)
print("mu2= ",mu2)
print("theta= ",θ)
print("theta=mu1*mu2.conjugate/N", θ==mu2*mu1.conjugate()/N1)
print("N1/N=",round(N1**5/N,5))
