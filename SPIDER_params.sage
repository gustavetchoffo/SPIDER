from sqisign_prism_v2_main import params, constants
from sqisign_prism_v2_main.ec import *
from two_isogenies.Theta_SageMath.utilities.supersingular import *
#SPIDER parameters

def set_spider_params(lvl):
    """
    Return parameters for SPIDER at the given level
    """

    global p, f, lb, a, rb, D_mix
    global FP_ENC_BYTES

    if lvl == 1:
        #p = ZZ(27*2**500 - 1)
        p = 11 * 2**257 * 3**163 - 1
        f = ZZ(257)
        a = ZZ(257)
        b=ZZ(163)
        rb = ZZ(32)
        lb = ZZ(128)
        
# Parameters for level 3 and 5 are still to be defined.
    elif lvl == 3:
        p = ZZ(65*2**376 - 1)
        f = ZZ(376)
        a = ZZ(256)
        rb = ZZ(20)
        lb = ZZ(376)

    elif lvl == 5:
        p = ZZ(27*2**500 - 1)
        f = ZZ(500)
        a = ZZ(500)
        rb = ZZ(32)
        lb = ZZ(128)

    else:
        raise ValueError(f"level {lvl} not recognized")

    D_mix = next_prime(2**(4*lb))
    FP_ENC_BYTES = 8 * floor((log(p, 2) + 63) / 64)

    # Quaternions
    global B, O0
    global QUAT_prime_cofactor
    B = QuaternionAlgebra(-1, -p)
    _i, _j, _k = B.gens()
    O0 = B.maximal_order(order_basis=(B(1), _i, (_i+_j)/2, (1-_k)/2))
    QUAT_prime_cofactor = next_prime(2**ceil(log(p, 2)))

    # Elliptic Curves
    global Fp, Fp2, Fp2_i
    global E0, P0, Q0
    Fp = GF(p)
    Fp2, Fp2_i = GF(p**2, name="i", modulus=[1, 0, 1]).objgen()
    E0 = EllipticCurve(Fp2, [1, 0])
    E0.set_order((p + 1) ** 2)

    # Hardcoded SQIsign implementation points and endomorphisms actions
    global iota, frob
    global mat_1, mat_i, mat_ij2, mat_1k2
    mat_1 = matrix(Zmod(2**f), 2, [1, 0, 0, 1])
    P0,Q0=torsion_basis(E0, 2**a)
    print("P0=",P0, "\n Q0=",Q0)

set_spider_params(1)
'''
    if lvl == 1:
        P0 = E0(constants.Px1, constants.Py1)
        Q0 = E0(constants.Qx1, constants.Qy1)
        mat_i = matrix(Zmod(2**f), 2, constants.mat_i_1)
        mat_ij2 = matrix(Zmod(2**f), 2, constants.mat_ij2_1)
        mat_1k2 = matrix(Zmod(2**f), 2, constants.mat_1k2_1)

# Parameters for level 3 and 5 are still to be defined.
    elif lvl == 3:
        P0 = E0(constants.Px3, constants.Py3)
        Q0 = E0(constants.Qx3, constants.Qy3)
        mat_i = matrix(Zmod(2**f), 2, constants.mat_i_3)
        mat_ij2 = matrix(Zmod(2**f), 2, constants.mat_ij2_3)
        mat_1k2 = matrix(Zmod(2**f), 2, constants.mat_1k2_3)

    elif lvl == 5:
        P0 = E0(constants.Px5, constants.Py5)
        Q0 = E0(constants.Qx5, constants.Qy5)
        mat_i = matrix(Zmod(2**f), 2, constants.mat_i_5)
        mat_ij2 = matrix(Zmod(2**f), 2, constants.mat_ij2_5)
        mat_1k2 = matrix(Zmod(2**f), 2, constants.mat_1k2_5)
'''