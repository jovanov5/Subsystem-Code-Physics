""" Measurement
    - methods rely on the correct definiton of each `system` object
"""

function entanglement_negativity(state::MixedDestabilizer,
                                system)
    """this calculates the TEN of the state on a torus, using the usual quadpartition. 
    state            : state to be evolved. It is taken by reference (hence the ! in the name by convention).
    system           : system object that has knows everyhting!
    """

    S = system.S
    P_A = system.P_A
    P_B = system.P_B
    P_C = system.P_C

    T = stab_to_gf2(stabilizerview(state))'

    K = T'*(P_A)*S*(P_A)*T
    K = K.%2
    K = K.!=0
    e_A = rank(K)/2
    K = T'*(P_B)*S*(P_B)*T
    K = K.%2
    K = K.!=0
    e_B = rank(K)/2
    K = T'*(P_C)*S*(P_C)*T
    K = K.%2
    K = K.!=0
    e_C = rank(K)/2
    K = T'*(P_A+P_B)*S*(P_A+P_B)*T
    K = K.%2
    K = K.!=0
    e_AB = rank(K)/2
    K = T'*(P_A+P_C)*S*(P_A+P_C)*T
    K = K.%2
    K = K.!=0
    e_AC = rank(K)/2
    K = T'*(P_B+P_C)*S*(P_B+P_C)*T
    K = K.%2
    K = K.!=0
    e_BC = rank(K)/2
    K = T'*(P_A+P_B+P_C)*S*(P_A+P_B+P_C)*T
    K = K.%2
    K = K.!=0
    e_ABC = rank(K)/2

    e_topo = -e_A-e_B-e_C+e_AB+e_AC+e_BC-e_ABC
    return e_topo
end

function entanglement_entropy_topo(state::MixedDestabilizer,
                            system)
    """this calculates the TEE of the state on a torus, using the usual quadpartition. 
    state            : state to be evolved. It is taken by reference (hence the ! in the name by convention).
    n                : total number of qubits.
    """

    r_A = system.r_A
    r_B = system.r_B
    r_C = system.r_C
    r_AB = system.r_AB
    r_BC = system.r_BC
    r_ABC = system.r_ABC

    e_A = entanglement_entropy(state,     # state to compute entropy for
                r_A,       # subsystem as list of indices. Specifying full subsystem gives the van-neumann entropy
                Val(:rref) # algorithm to use (see documentation)
                )
    e_B = entanglement_entropy(state,     # state to compute entropy for
                r_B,       # subsystem as list of indices. Specifying full subsystem gives the van-neumann entropy
                Val(:rref) # algorithm to use (see documentation)
                )
    e_C = entanglement_entropy(state,     # state to compute entropy for
                r_C,       # subsystem as list of indices. Specifying full subsystem gives the van-neumann entropy
                Val(:rref) # algorithm to use (see documentation)
                )
    e_AB = entanglement_entropy(state,     # state to compute entropy for
                r_AB,       # subsystem as list of indices. Specifying full subsystem gives the van-neumann entropy
                Val(:rref) # algorithm to use (see documentation)
                )
    e_BC = entanglement_entropy(state,     # state to compute entropy for
                r_BC,       # subsystem as list of indices. Specifying full subsystem gives the van-neumann entropy
                Val(:rref) # algorithm to use (see documentation)
                )
    e_ABC = entanglement_entropy(state,     # state to compute entropy for
                r_ABC,       # subsystem as list of indices. Specifying full subsystem gives the van-neumann entropy
                Val(:rref) # algorithm to use (see documentation)
                )
    e_AC = entanglement_entropy(state,     # state to compute entropy for
                [r_A; r_C],       # subsystem as list of indices. Specifying full subsystem gives the van-neumann entropy
                Val(:rref) # algorithm to use (see documentation)
                )

    e_topo = -e_A-e_B-e_C+e_AB+e_AC+e_BC-e_ABC
    # print(e_topo)
    return e_topo
end

function entanglement_entropy_cut(state::MixedDestabilizer, system, subdiv_array)
    """this calculates the entanglement entropy along a cut of the state on a torus, using the anular subdivisions.
    state            : state to be examined.
    system           : system object that has knows everything!
    n_sub            : number of subdivisions.
    """

    L = system.L
    k = system.k
    ee_array = []
    for sub_size in subdiv_array
        subsys_rande = sub_size*L*k
        e = entanglement_entropy(state,     # state to compute entropy for
                1:subsys_rande,       # subsystem as list of indices. Specifying full subsystem gives the van-neumann entropy
                Val(:rref) # algorithm to use (see documentation)
                )
        push!(ee_array, e)
    end
    return ee_array
end