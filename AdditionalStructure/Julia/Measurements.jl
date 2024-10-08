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
    subdiv_array     : where to cut the system.
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

function entanglement_entropy_half_cut(state::MixedDestabilizer, system)
    """this calculates the entanglement entropy along a cut of the state on a torus, using the anular subdivisions, at the middle.
    state            : state to be examined.
    system           : system object that has knows everything!
    """

    L = system.L
    sub_size = Int(round(L/2))
    k = system.k
    subsys_rande = sub_size*L*k
    e = entanglement_entropy(state,     # state to compute entropy for
            1:subsys_rande,       # subsystem as list of indices. Specifying full subsystem gives the van-neumann entropy
            Val(:rref) # algorithm to use (see documentation)
            )
    return e
end

function general_zassenhausen_correlator(state::MixedDestabilizer, system, rep_funtion::Function, deformator)
    """this calculates the represented correlator between two sites (0, 0) and (all, all) of the state on a torus.
    state            : state to be examined.
    system           : system object that has knows everything!
    rep_funtion      : function that returns the representative string operator.
    deformator       : deformator to the representative
    """

    L = system.L
    k = system.k
    nbits = system.nbits
    cell_num = system.cell_num

    corr_function = zeros(Int, L, L)

    for r_x  = 0:L-1
        for r_y = 0:L-1
            if r_x == 0 && r_y == 0
                corr_function[r_x+1, r_y+1] = 1
            else
                # rep_function is an input
                rep = rep_funtion((r_x, r_y), system)
                # size_of_bff = 0
                # zassenhausen_alg also in HelperTools.jl
                dim_source, _ = zassenhausen_alg(rep, deformator, stabilizerview(state), nbits)
                dim_norm, _ = zassenhausen_alg(deformator, stabilizerview(state), nbits)
                corr_function[r_x+1, r_y+1] = Int(round(dim_source - dim_norm))
            end
        end
    end

    return corr_function

end

function general_zassenhausen_correlator_fast(state::MixedDestabilizer, system, rep_funtion::Function, def_function_fast::Function)
    """this calculates the represented correlator between two sites (0, 0) and (all, all) of the state on a torus.
    state            : state to be examined.
    system           : system object that has knows everything!
    rep_funtion      : function that returns the representative string operator.
    deformator       : deformator to the representative
    """

    L = system.L
    k = system.k
    nbits = system.nbits
    cell_num = system.cell_num

    corr_function = zeros(Int, L, L)

    for r_x  = 0:L-1
        for r_y = 0:L-1
            if r_x == 0 && r_y == 0
                corr_function[r_x+1, r_y+1] = 1
            else
                # Removes non-commuting elements from the general (deformator, representative) group so we do not need to doible the system size.
                deformator = def_function_fast((r_x, r_y), system)
                # rep_function is an input
                rep = rep_funtion((r_x, r_y), system)
                # size_of_bff = 0
                # zassenhausen_alg also in HelperTools.jl
                dim_source, _ = zassenhausen_alg_fast(rep, deformator, stabilizerview(state), nbits)
                dim_norm, _ = zassenhausen_alg_fast(deformator, stabilizerview(state), nbits)
                corr_function[r_x+1, r_y+1] = Int(round(dim_source - dim_norm))
            end
        end
    end

    return corr_function

end

function particular_zassenhausen_correlator(r_2::Tuple{Int, Int}, state::MixedDestabilizer, system, rep_funtion::Function, deformator)
    """this calculates the represented correlator between two sites (0, 0) and (i, j) of the state on a torus.
    r_2              : coordinate of the second site.
    state            : state to be examined.
    system           : system object that has knows everything!
    rep_funtion      : function that returns the representative string operator.
    deformator       : deformator to the representative
    """

    L = system.L
    k = system.k
    nbits = system.nbits
    cell_num = system.cell_num

    r_x = r_2[1]
    r_y = r_2[2]

    if r_x == 0 && r_y == 0
        return 1
    else
        rep = rep_funtion((r_x, r_y), system)
        # size_of_bff = 0
        dim_source, _ = zassenhausen_alg(rep, deformator, stabilizerview(state), nbits)
        dim_norm, _ = zassenhausen_alg(deformator, stabilizerview(state), nbits)
        return Int(round(dim_source - dim_norm))
    end
end

function particular_zassenhausen_correlator_fast(r_2::Tuple{Int, Int}, state::MixedDestabilizer, system, rep_funtion::Function, def_function_fast::Function)
    """this calculates the represented correlator between two sites (0, 0) and (i, j) of the state on a torus.
    r_2              : coordinate of the second site.
    state            : state to be examined.
    system           : system object that has knows everything!
    rep_funtion      : function that returns the representative string operator.
    deformator       : deformator to the representative
    """

    L = system.L
    k = system.k
    nbits = system.nbits
    cell_num = system.cell_num

    r_x = r_2[1]
    r_y = r_2[2]

    if r_x == 0 && r_y == 0
        return 1
    else
        deformator = def_function_fast((r_x, r_y), system)
        rep = rep_funtion((r_x, r_y), system)
        # size_of_bff = 0
        dim_source, _ = zassenhausen_alg_fast(rep, deformator, stabilizerview(state), nbits)
        dim_norm, _ = zassenhausen_alg_fast(deformator, stabilizerview(state), nbits)
        return Int(round(dim_source - dim_norm))
    end
end

function get_e_reprentative(r_2::Tuple{Int, Int}, system)
    """ This gets you a PauliOperator of one string connecting two e-type defects between sites (0,0) and r_2.
    """

    nbits = system.nbits
    L = system.L
    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)

    for cell_index = 0:r_2[1]-1
        Z_op_list[2*cell_index + 1 + 1] = true
    end

    for cell_index = 0:r_2[2]-1
        Z_op_list[2*(L*cell_index + r_2[1]) + 1] = true
    end

    return PauliOperator(0x0, X_op_list, Z_op_list)
end

function get_m_reprentative(r_2::Tuple{Int, Int}, system)
    """ This gets you a PauliOperator of one string connecting two m-type defects between sites (0,0) and r_2.
    """

    nbits = system.nbits
    L = system.L
    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)

    for cell_index = 0:r_2[1]-1
        X_op_list[2*(cell_index + 1) + 1] = true
    end

    for cell_index = 0:r_2[2]-1
        X_op_list[2*(L*(cell_index + 1) + r_2[1]) + 1 + 1] = true
    end

    return PauliOperator(0x0, X_op_list, Z_op_list)
end


function get_f_reprentative(r_2::Tuple{Int, Int}, system)
    """ This gets you a PauliOperator of one string connecting two f-type defects between sites (0,0) and r_2.
    """

    nbits = system.nbits
    L = system.L
    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)

    for cell_index = 0:r_2[1]-1
        Z_op_list[2*cell_index + 1 + 1] = true
        X_op_list[2*(cell_index + 1) + 1] = true
    end

    for cell_index = 0:r_2[2]-1
        Z_op_list[2*(L*cell_index + r_2[1]) + 1] = true
        X_op_list[2*(L*(cell_index + 1) + r_2[1]) + 1 + 1] = true
    end

    return PauliOperator(0x0, X_op_list, Z_op_list)
end

function get_e_deformator(system)
    nbits = system.nbits
    L = system.L
    e_deformators = []

    for cell_index = 0:L*L-1
        push!(e_deformators, tc_stab(Plaquette, cell_index, system))
    end

    # The noncontractible loops are actually very important deformators.

    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)
    for cell_index = 0:L-1
        Z_op_list[2*cell_index+1+1] = true
    end
    Z_loop_1 = PauliOperator(0x0, X_op_list, Z_op_list)
    push!(e_deformators, Z_loop_1) # Up Loop

    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)
    for cell_index = 0:L-1
        Z_op_list[2*L*cell_index+1] = true
    end
    Z_loop_2 = PauliOperator(0x0, X_op_list, Z_op_list)
    push!(e_deformators, Z_loop_2) # Right Loop

    return e_deformators

end

function get_e_deformator_fast(r_2::Tuple{Int, Int}, system) # no changes needed, This group is abelian!
    nbits = system.nbits
    L = system.L
    e_deformators = []

    for cell_index = 0:L*L-1
        push!(e_deformators, tc_stab(Plaquette, cell_index, system))
    end

    # The noncontractible loops are actually very important deformators.

    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)
    for cell_index = 0:L-1
        Z_op_list[2*cell_index+1+1] = true
    end
    Z_loop_1 = PauliOperator(0x0, X_op_list, Z_op_list)
    push!(e_deformators, Z_loop_1) # Up Loop

    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)
    for cell_index = 0:L-1
        Z_op_list[2*L*cell_index+1] = true
    end
    Z_loop_2 = PauliOperator(0x0, X_op_list, Z_op_list)
    push!(e_deformators, Z_loop_2) # Right Loop

    return e_deformators

end

function get_m_deformator(system)
    nbits = system.nbits
    L = system.L
    m_deformators = []

    for cell_index = 0:L*L-1
        push!(m_deformators, tc_stab(Star, cell_index, system))
    end

    # The noncontractible loops are actually very important deformators.

    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)
    for cell_index = 0:L-1
        X_op_list[2*cell_index+1] = true
    end
    X_loop_1 = PauliOperator(0x0, X_op_list, Z_op_list)
    push!(m_deformators, X_loop_1) # Up Loop

    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)
    for cell_index = 0:L-1
        X_op_list[2*L*cell_index+1+1] = true
    end
    X_loop_2 = PauliOperator(0x0, X_op_list, Z_op_list)
    push!(m_deformators, X_loop_2) # Right Loop

    return m_deformators

end

function get_f_deformator(system)

    nbits = system.nbits
    L = system.L
    f_deformators = []

    for cell_index = 0:L*L-1
        push!(f_deformators, tc_stab(Star, cell_index, system))
    end

    for cell_index = 0:L*L-1
        push!(f_deformators, tc_stab(Plaquette, cell_index, system))
    end

    # The noncontractible loops are actually very important deformators.

    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)
    for cell_index = 0:L-1
        Z_op_list[2*cell_index+1+1] = true
    end
    Z_loop_1 = PauliOperator(0x0, X_op_list, Z_op_list)
    push!(f_deformators, Z_loop_1) # Up Loop

    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)
    for cell_index = 0:L-1
        Z_op_list[2*L*cell_index+1] = true
    end
    Z_loop_2 = PauliOperator(0x0, X_op_list, Z_op_list)
    push!(f_deformators, Z_loop_2) # Right Loop

    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)
    for cell_index = 0:L-1
        X_op_list[2*cell_index+1] = true
    end
    X_loop_1 = PauliOperator(0x0, X_op_list, Z_op_list)
    push!(f_deformators, X_loop_1) # Up Loop

    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)
    for cell_index = 0:L-1
        X_op_list[2*L*cell_index+1+1] = true
    end
    X_loop_2 = PauliOperator(0x0, X_op_list, Z_op_list)
    push!(f_deformators, X_loop_2) # Right Loop

    return f_deformators

end

function get_f_deformator_fast(r_2::Tuple{Int, Int}, system)

    nbits = system.nbits
    L = system.L
    f_deformators = []

    problematic_index = r_2[1] + r_2[2]*L

    for cell_index = 1:L*L-1
        if cell_index != problematic_index
            push!(f_deformators, tc_stab(Star, cell_index, system))
        end
    end

    for cell_index = 0:L*L-1
        if cell_index != problematic_index
            push!(f_deformators, tc_stab(Plaquette, cell_index, system))
        end
    end

    push!(f_deformators, tc_stab(Plaquette, problematic_index, system) * tc_stab(Star, problematic_index, system))
    push!(f_deformators, tc_stab(Plaquette, 0, system) * tc_stab(Star, 0, system))

    # The noncontractible loops are actually very important deformators.

    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)
    for cell_index = 0:L-1
        Z_op_list[2*cell_index+1+1+L] = true # I want the Z to be to the left to X so that proper fermion symmetry deforms this string!
    end
    Z_loop_1 = PauliOperator(0x0, X_op_list, Z_op_list)
    # push!(f_deformators, Z_loop_1) # Up Loop

    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)
    for cell_index = 0:L-1
        Z_op_list[2*L*cell_index+1+2] = true # Simaraly, I want Z on top of Xs in the horizontal loop!
        # The only issue is that the deformation power of this group is not adequate, but let us see how it goes!!! ;)
    end
    Z_loop_2 = PauliOperator(0x0, X_op_list, Z_op_list) 
    # push!(f_deformators, Z_loop_2) # Right Loop

    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)
    for cell_index = 0:L-1
        X_op_list[2*cell_index+1] = true
    end
    X_loop_1 = PauliOperator(0x0, X_op_list, Z_op_list)
    push!(f_deformators, X_loop_1 * Z_loop_1) # Up Loop

    Z_op_list = zeros(Bool, nbits)
    X_op_list = zeros(Bool, nbits)
    for cell_index = 0:L-1
        X_op_list[2*L*cell_index+1+1] = true
    end
    X_loop_2 = PauliOperator(0x0, X_op_list, Z_op_list)
    push!(f_deformators, X_loop_2 * Z_loop_2) # Right Loop

    return f_deformators

end
