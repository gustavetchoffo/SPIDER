 '''
    e=pp.e
    d=pp.d
    q=pp.q
    precomps=pp.prism_param
    P0,Q0=precomps.basis
    I0=precomps.random_ideal_of_given_norm(N, prime=True)
    phi0=IdealToIsogenyQlapotis(I0, e, precomps)
    E,P,Q=phi0.images()                #Qlapoti
    J=precomps.random_ideal_of_given_norm(d, prime=False)
    Ir = pushforward(J,I0)
    
    K = I0 * Ir
    #print(K.left_order(),"\n K.left_order()==O0: ",K.left_order()==O0)
    #Iaux=I0.intersection(J)
    #print("\n K==Iaux",K==Iaux,Ir.norm()==d)
    phi_phi0=IdealToIsogenyQlapotis(K, e, precomps)
    Eaux,Paux,Qaux=phi_phi0.images()
    print("\n Norm of K=d.N: ",K.norm()==N*d)

    ker_Phi=((N*P0,Paux),(N*Q0,Qaux))
    #Phi = Dim2Isogeny(ker_Phi, a)

    prism=PRISM_sign(role='signer', precomps=None, level=5, qlapoti=True,SPIDER=True)
    msg='sample SPIDER'
    prism.keygen()
    E=prism.pk_curve
    I0,M=prism.sk
    Eprim,(Pr,Qr),q=prism.sign(msg.encode())
    Phi=Eprim,(Pr,Qr)
    return E,Phi,q
    '''
def commit_1(pp,E,Phi,q):
    pp=pp
    q=q
    
    Eprim,(Pr,Qr)=Phi
    a=pp.a
    d=q*(2**(a)-q)
    P,Q,w=torsion_basis_with_pairing_2e(E, 2**a)
    ker_Phi=((q*P,Pr),(q*Q,Qr))
    Phi = Dim2Isogeny(ker_Phi, sign_torsion)
    r=pp.r
    dcom=2**a
    dcom_prev=2**(a-1)
    
    s=[]
    for _ in range(r):
        s.append(ZZ.random_element(1,dcom))

    E1=E
    Qu=E1(compute_point_order_2e(E1, a, x_start=0))
    #print("Qu= ",s[0]*Qu)
    L1=KummerLine(E1)
    for u in range(1,r+1):
        Pu,_=compute_linearly_independent_point_with_pairing_2e(E1, Qu, a, x_start=0)
        #print("Pu= ",Pu)
        Ru=Pu+s[u-1]*Qu
        if dcom_prev*Ru==E1((0,0)):
            if s[u-1]!=0:
                Ru=Pu
        xRu=L1(Ru[0])
        print("xRu.order()==dcom ",Ru.order()==dcom)
        phi_com=KummerLineIsogeny(L1,xRu,dcom)
        xQu=L1(Qu[0])
        x=phi_com(xQu)
        E1=phi_com.codomain().curve()
        L1=KummerLine(E1)
        Qu=x.curve_point()
    return s,E1
    