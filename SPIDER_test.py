import sys
import logging
import time
logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()

from SPIDER_signature import *


if __name__ == "__main__":
    
    q=next_prime(2**249)
    #a=320
    #t,w=247,30
    t=Integer(210);w=Integer(33)
    #t, w=222, 32
    pp=SPIDER_pp(q,level=1,r=1,t=t,w=w)
    d=pp.d
    a=pp.a/2
    msg="hallo world"
    #Fp2=params.Fp2
    #print("\n",sample(pp))
    t1=time.time()
    E,phi=spider_sample(pp)
    t2=time.time()
    print(f'SPIDER parameters:\n lambda={pp.lb}\n p={factor(pp.p+1)}-1 \n a={a}\n t={pp.t}\n w={pp.w}')
    t3=time.time()
    sign=SPIDER_sign(pp,E,phi,msg)
    t4=time.time()
    v=SPIDER_vrify(pp,E,msg,sign)
    t5=time.time()
    print(f'key generation time: {t2-t1:.3f} s')
    print(f'Signature time: {t4-t3:.3f} s')
    print(f'verification time: {t5-t4:.3f} s')
    
    sign_bytes=sys.getsizeof(sign) + sum(sys.getsizeof(s) for s in sign)
    usk_bytes=sys.getsizeof(phi) + sum(sys.getsizeof(s) for s in phi)
    print("\n verif=",v,'\n signature size=',sign_bytes,'B\n','Mkey size=',sys.getsizeof(E),'B\n','B\n','Usk size=',usk_bytes,'B\n' )
