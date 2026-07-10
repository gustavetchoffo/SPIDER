from sage.all_cmdline import *   # import sage library

_sage_const_1 = Integer(1); _sage_const_3 = Integer(3); _sage_const_5 = Integer(5); _sage_const_2 = Integer(2); _sage_const_0 = Integer(0); _sage_const_210 = Integer(210); _sage_const_33 = Integer(33); _sage_const_500 = Integer(500); _sage_const_4 = Integer(4)
import hashlib
import logging
import time
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()
formatter = logging.Formatter('%(name)s [%(levelname)s] %(message)s')
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)

from SPIDER_sigma import *
from generalities import SeedTree,Node_seed,recover_leaves,parse_hashs_t_w

class SPIDER_pp():
    '''
    pp for SPIDER signature '''
    def __init__(self,q,level=1,r=1):
        self.pp_sigma=parameters(q,level,r)
        if level==1:
            self.t=210
            self.w=33
    
    def H(self,data):
        #msg=data.encode()
        ch=parse_hashs_t_w(data, 2, self.t, self.w)
        ch=[ZZ(mod(ch[i]+1,2)) for i in range(self.t)]
        return ch



def SPIDER_sign(pp,pk,sk,msg):
    print("SPIDER signing process")
    E=pk
    phi=sk
    t=pp.t
    lb=pp.pp_sigma.lb
    size=Integer(2**lb)
    val_seed_root=ZZ.random_element(size)
    seed_root=Node_seed(val_seed_root,parent=None,left_child=None,right_child=None,h=0,i=0,nb_leaves=t)
    seed_tree=SeedTree(seed_root,t,lb)
    seed_leaves=seed_tree.leaves()
        
    seeds=[int.from_bytes(leaf.value) for leaf in seed_leaves]
    #seeds=[int.from_bytes(val) for val in seed_bytes]
    #com=[]
    state=[]
    data=E.j_invariant().to_bytes()
    for i in range(t):
        com_i=sigma_Rd_Commit(pp.pp_sigma,E,phi,seeds[i])
        #com.append(com_i[0])
        state.append(com_i[1])
        data+=com_i[0].j_invariant().to_bytes()      #preferable to use the A- coeff?
    msg=msg.encode()
    data+=msg
    ch=pp.H(data)
    seed_internal=seed_tree.releaseSeed(ch,0)
    resp={}
    for i in range(t):
        if ch[i]==1:
            resp[i]=sigma_Rd_Response(state[i],1)
    return ch,seed_internal,resp

def SPIDER_vrify(pp,pk,msg,sign):
    print("SPIDER verification process")
    E=pk
    t=pp.t
    pp_sigma=pp.pp_sigma
    lb=pp.pp_sigma.lb
    ch,seed_internal,resp=sign
    seeds_0=recover_leaves(seed_internal,ch,0,lb)
    ct=0
    for seed in seeds_0:
        for i in range(ct,t):
            if ch[i]==0:
                resp[i]=seed
                ct=i+1
                break
    assert len(resp)==t
    data=E.j_invariant().to_bytes()
    for i in range(t):
        if ch[i]==0:
            seed=int.from_bytes(resp[i])
            s_u=pp_sigma.PRNG(seed)
            P,Q=TorsionBasis(E)
            Pu,Qu=2**a*P, 2**a*Q
            Ru=Pu+s_u[0]*Qu
            psi0=E.isogeny(Ru,algorithm="factored")
            j1=psi0.codomain().j_invariant()
            data+=j1.to_bytes()
        else:
            E1,_,_=resp[i]
            assert sigma_Rd_Verify(pp_sigma,E,E1,1,resp[i])
            j1=E1.j_invariant()
            data+=j1.to_bytes()
    msg=msg.encode()
    data+=msg
    return pp.H(data)==ch



if __name__ == "__main__":
    
    q=next_prime(2**249)
    #a=320
    pp=SPIDER_pp(q,level=1,r=1)
    pp_sigma=pp.pp_sigma
    d=pp_sigma.d
    a=pp_sigma.a/2
    msg="hallo world"
    #Fp2=params.Fp2
    #print("\n",sample(pp))
    t1=time.time()
    E,phi=spider_sample(pp_sigma)
    t2=time.time()
    
    sign=SPIDER_sign(pp,E,phi,msg)
    t3=time.time()
    v=SPIDER_vrify(pp,E,msg,sign)
    t4=time.time()
    print(f'key generation time: {t2-t1:.3f} s')
    print(f'Signature time: {t3-t2:.3f} s')
    print(f'verification time: {t4-t3:.3f} s')
    print("\n verif=",v)
