import sys
import time
from sqisign_prism_v2_main import params
#except ImportError:
 #   import params

from sqisign_prism_v2_main.quaternions import *
from sqisign_prism_v2_main.qlapoti import IdealToIsogeny
from sqisign_prism_v2_main.theta.theta_structures.couple_point import *
from sqisign_prism_v2_main.hd import Dim2Iso
from sqisign_prism_v2_main.ec import *
from sqisign_prism_v2_main.prism import hash_to_prime
#from SPIDER_sigma import parameters,spider_sample,sigma_Rd_Commit,sigma_Rd_Response,sigma_Rd_Verify
from SPIDER_signature import SPIDER_pp,SPIDER_sign,SPIDER_vrify

class parameters_IBS(SPIDER_pp):
    '''
    pp for SPIDER IBS protocol '''
    
    def __init__(self):
        super().__init__(q=1,level=1,r=1)

    def H_prism(self,id,mpk,r=None):
        a=self.a/2
        return hash_to_prime(id,mpk,r=r,a=a)

def MkeyGen(pp):
    print("SPIDER IBS master key generation process")
    p=pp.p
    E0=pp.E0
    O0=pp.O0
    d=pp.d
    q=pp.q
    lb=pp.lb 
    P0=params.P0
    Q0=params.Q0

    a=pp.a/_sage_const_2 
    N=next_prime(p**2)
    I0=RandomIdealGivenNorm(N, prime=True)
    
    if not I0:
        raise ValueError(f"Failed to find a random ideal of norm {N}")
    E,phi_P0,phi_Q0=IdealToIsogeny(I0)
    E.set_order((p+1 )**2 , num_checks=0 )
    P,Q=TorsionBasis(E)
    M=ChangeOfBasis((phi_P0,phi_Q0),(P,Q), e=None) 
    assert EvalMatrix(M, (phi_P0, phi_Q0)) == (P, Q)
    mpk=E
    msk=M,I0
    return mpk,msk

def KeyDer(pp,mpk,msk,id):
    print("SPIDER IBS key derivation process")
    "Derivation of the user key from the master secret key and the identity"
    p=pp.p
    E0=pp.E0
    O0=pp.O0
    lb=pp.lb 
    P0=params.P0
    Q0=params.Q0
    q,ct=pp.H_prism(id,mpk)
    pp.set_q(q)
    d=pp.d
    M,I0=msk
    J=RandomIdealGivenNorm(d, prime=False)
    if not J:
        raise ValueError(f"Failed to find a random ideal of norm {d}")
    I=Pushforward(J,I0)
    I1=I0*I
    E_prim,P_prim0,Q_prim0=IdealToIsogeny(I1)
    (phi_P,phi_Q)=EvalMatrix(M, basis=(P_prim0, Q_prim0))
    E_prim.set_order((p+1 )**2 , num_checks=0 )
    P1,Q1=TorsionBasis(E_prim)
    M1=ChangeOfBasis((P1,Q1),(phi_P,phi_Q), e=None) 
    assert EvalMatrix(M1, (P1, Q1)) == (phi_P, phi_Q)
    
    
    usk=E_prim,M1,ct
    return usk

def sign_IBS(pp,mpk,id,usk,msg):
    print("SPIDER IBS signing process")
    E_prim,M1,ct=usk
    q,ct=pp.H_prism(id,mpk,r=ct)
    if pp.q!=q: pp.set_q(q)
    
    phi_P,phi_Q=EvalMatrix(M1, basis=TorsionBasis(E_prim))
    phi=E_prim,(phi_P,phi_Q)
    sig=SPIDER_sign(pp,mpk,phi,msg),ct
    return  sig

def verify_IBS(pp,mpk,id,msg,sig):
    print("SPIDER IBS verification process")
    "Verifying a signature with the master public key and the identity"
    q,_=pp.H_prism(id,mpk,r=sig[1])
    if pp.q!=q: pp.set_q(q)
    v=SPIDER_vrify(pp,mpk,msg,sig[0])
    return v


#==================================================================================
#========================= TESTING =========================================================
#==================================================================================
pp=parameters_IBS()
t0=time.time()
mpk,msk=MkeyGen(pp)
t1=time.time()
print(f"Master key generation time: {t1-t0:.3f} s")

id="Gustave Tchoffo"
msg="Hello world"
t2=time.time()
usk=KeyDer(pp,mpk,msk,id)
t3=time.time()
sig=sign_IBS(pp,mpk,id,usk,msg)
t4=time.time()
ver=verify_IBS(pp,mpk,id,msg,sig)
t5=time.time()
print("ver=",ver,'\n')
print(f"User key derivation time: {t3-t2:.3f} s")
print(f"Signature time: {t4-t3:.3f} s")
print(f"Verification time: {t5-t4:.3f} s\n")
print(f'mpk size: {sys.getsizeof(mpk)} B')
print(f'msk size: {(sys.getsizeof(msk)+ sum(sys.getsizeof(s) for s in msk))} B')
print(f'Usk size: {(sys.getsizeof(usk)+ sum(sys.getsizeof(s) for s in usk))} B')
print(f'signature size: {(sys.getsizeof(sig[0])+ sum(sys.getsizeof(s) for s in sig[0])+sys.getsizeof(sig[1]))/1024.0: .2f} KB')