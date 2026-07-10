import hashlib
import logging
import time
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()
formatter = logging.Formatter('%(name)s [%(levelname)s] %(message)s')
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)

from sage.all import *
proof.all(False)

#try:
from sqisign_prism_v2_main import params
#except ImportError:
 #   import params

from sqisign_prism_v2_main.quaternions import *
from sqisign_prism_v2_main.qlapoti import IdealToIsogeny
from sqisign_prism_v2_main.theta.theta_structures.couple_point import *
from sqisign_prism_v2_main.hd import Dim2Iso
from sqisign_prism_v2_main.ec import *


class parameters():
    '''
    pp for SPIDER sigma protocol '''
    def __init__(self,q,level=1,r=1):
        assert level in [1,3,5]
        params.set_spider_params(level)
        self.r=r
        self.q=q
        self.a=params.a
        self.d=q*(2**(self.a/2)-q)
        assert self.d>0, "d must be positive"
        self.p=params.p
        self.E0=params.E0
        self.O0=params.O0
        
        self.lb=params.lb
        self.rb=params.rb
        Fq=self.E0.base_ring()
        Fp4 = Fq.extension(2, 'j')
        self.Fp4 = Fp4
        
          

    def PRNG(self, seed):
        s_u=[]
        set_random_seed(seed)
        for _ in range(self.r):
            s_u.append(ZZ.random_element(1,2**self.lb))
        return s_u
          

def spider_sample(pp):
    '''
    return a HD representation of phi:E--->E' of degree d=q(2^a-q).
    Note that we return (phi(P),phi(Q))= phi(E[2^(2a)]) instead of phi(E[2^(a)]).
    Tho comput the isogeny dimon, one should use (R=[2^a]P, S=[2^a]Q) and (phi(R),phi(S))
    '''
    p=pp.p
    E0=pp.E0
    O0=pp.O0
    d=pp.d
    q=pp.q
    lb=500
    P0=params.P0
    Q0=params.Q0
    
    a=pp.a/2
    N=next_prime(2**(4*lb))
    I0=RandomIdealGivenNorm(N, prime=True)
    if not I0:
        raise ValueError(f"Failed to find a random ideal of norm {N}")
    E,phi_P0,phi_Q0=IdealToIsogeny(I0)
    P,Q=TorsionBasis(E)
    M=ChangeOfBasis((phi_P0,phi_Q0),(P,Q), e=None) 
    assert EvalMatrix(M, (phi_P0, phi_Q0)) == (P, Q)
    J=RandomIdealGivenNorm(d, prime=False)
    if not J:
        raise ValueError(f"Failed to find a random ideal of norm {d}")
    I=Pushforward(J,I0)
    I1=I0*I
    E_prim,P_prim0,Q_prim0=IdealToIsogeny(I1)
    (phi_P,phi_Q)=EvalMatrix(M, basis=(P_prim0, Q_prim0))
    phi=E_prim,(phi_P,phi_Q)  
    return E,phi

def sigma_Rd_Commit(pp,E,phi,seed=None):
    a=pp.a/2   # This is because the actual a used in the code is 2a from the paper, for simplicity. We should fix this in the future.
    if seed is None:
        seed = randint(0, 2**a - 1)
    r=pp.r  
    assert r==1, "Currently only supports r=1" #TODO: Handle the case r>1 in the future
    s_u=pp.PRNG(seed)
    P,Q=TorsionBasis(E)
    
    Pu=(2**a)*P
    Qu=(2**a)*Q
    Ru=Pu+s_u[0]*Qu
    psi0=E.isogeny(Ru,algorithm="factored")
    E1=psi0.codomain().montgomery_model()
    state=(pp,seed,E,P,Q,a,phi,s_u,psi0,E1)
    return E1,state

def sigma_Rd_Response(state,ch):
    (pp,seed,E,P,Q,a,phi,s_u,psi0,E1)=state
    if ch==0:
        return seed
    else:
        E_prim,(phi_P,phi_Q)=phi
        psi0_dual=psi0.dual()
        P1,Q1=TorsionBasis(E1)
        iso=E1.isomorphism_to(psi0.codomain())
        psi0_dual_P1=psi0_dual(iso(P1))
        psi0_dual_Q1=psi0_dual(iso(Q1))
        A=ChangeOfBasis((P,Q),(psi0_dual_P1,psi0_dual_Q1), e=None)   # (psi0_dual_P1,psi0_dual_Q1)= A*(P,Q)
        assert EvalMatrix(A, (P, Q)) == (psi0_dual_P1, psi0_dual_Q1)
        (phi_psi0_dual_P1, phi_psi0_dual_Q1)= EvalMatrix(A, (phi_P,phi_Q))
        Ru_prim=2**a*phi_P+s_u[0]*2**a*phi_Q
        psi0_prim=E_prim.isogeny(Ru_prim,algorithm="factored")
        E1_prim=psi0_prim.codomain().montgomery_model()
        iso1=psi0_prim.codomain().isomorphism_to(E1_prim)
        psi1_P1=iso1(psi0_prim(phi_psi0_dual_P1))
        psi1_Q1=iso1(psi0_prim(phi_psi0_dual_Q1))

        return E1, E1_prim, (psi1_P1, psi1_Q1)
        
def sigma_Rd_Verify(pp,E,com,ch,resp):
    a=pp.a/2
    E1=com
    if ch==0:
        seed=resp
        s_u=pp.PRNG(seed)
        P,Q=TorsionBasis(E)
        Pu,Qu=2**a*P, 2**a*Q
        Ru=Pu+s_u[0]*Qu
        psi0=E.isogeny(Ru,algorithm="factored")
        return psi0.codomain().is_isomorphic(E1) 
    else:
        q=pp.q
        a=pp.a/2   # This is because the actual a used in the code is 2a from the paper.
        E1, E1_prim, (psi1_P1, psi1_Q1)=resp
        P1,Q1=TorsionBasis(E1)
        P1_u,Q1_u=2**a*P1, 2**a*Q1
        K=(CouplePoint(q*P1_u,psi1_P1),CouplePoint(q*Q1_u,psi1_Q1))
        PSI1=Dim2Iso(K,a)
        R1=PSI1(CouplePoint(P1_u,E1_prim(0)))[0]
        S1=PSI1(CouplePoint(Q1_u,E1_prim(0)))[0]
        w1= pari.ellweilpairing(E1, P1_u, Q1_u, 2**a)
        pair=pari.ellweilpairing(R1.curve(), R1, S1, 2**a)
        w1_q=w1**q
        w1_q_inv=w1_q**(-1)
        return pair in [w1_q,w1_q_inv] 

'''
if __name__ == "__main__":
    
    q=next_prime(2**249)
    #a=320
    pp=parameters(q,level=1,r=1)
    d=pp.d
    a=pp.a/2
    Fp2=params.Fp2
    #print("\n",sample(pp))
    E,phi=spider_sample(pp)
    E_prim,(phi_P,phi_Q)=phi
    P,Q=TorsionBasis(E)
    R=2**a*P
    S=2**a*Q
    phi_R=2**a*phi_P
    phi_S=2**a*phi_Q
    Ka=(CouplePoint(q*R,phi_R),CouplePoint(q*S,phi_S)) 
    PHI=Dim2Iso(Ka,a)

    KL_E=KummerLine(E)
    KP_P, KQ_P = KL_E(P), KL_E(Q)
    #print("\n KP_P=",KP_P.XZ(), "\n KQ_P=",KQ_P.XZ())
    T=CouplePoint(P,E_prim(0))
    print("\n PHI(T)=",PHI(T))
    #print("\n E_prim=",E_prim)  
    seed=randint(0, 2**a - 1)
    #seed=1226161722578291778072544579321671146702749075272823018708845811560162129143444254133558177112908189115150677810277390656415346694509704997388011667993
    #t1=time.time()
    E1,state=sigma_Rd_Commit(pp,E,phi,seed)
    #print("P0.order()=",factor(params.P0.order()), "\n Q0.order()=",factor(params.Q0.order()))
    resp=sigma_Rd_Response(state,1)
    
    b=sigma_Rd_Verify(pp,E,E1,1,resp)

    print("\n verif=",b)
'''

