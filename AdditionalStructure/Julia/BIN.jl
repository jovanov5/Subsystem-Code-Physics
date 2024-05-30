function fermion_correlator(state::MixedDestabilizer, system, f_deformator)
    """this calculates the fermionic correlator between two sites i and j of the state on a torus.
    state            : state to be examined.
    system           : system object that has knows everything!
    i, j             : indices of the sites to be correlated.
    """

    L = system.L
    k = system.k
    nbits = system.nbits
    cell_num = system.cell_num

    corr_function = zeros(Float64, L, L)

    for r_x  = 0:L-1
        for r_y = 0:L-1
            if r_x == 0 && r_y == 0
                corr_function[r_x+1, r_y+1] = 1
            else
                f_representative = get_f_reprentative((r_x, r_y), system)
                size_of_bff = 0
                dim_source, _ = zassenhausen_alg(e_representative, e_deformator, stabilizerview(state), nbits)
                dim_norm, _ = zassenhausen_alg(e_deformator, stabilizerview(state), nbits)
                corr_function[r_x+1, r_y+1] = dim_source - dim_norm
            end
        end
    end

    return corr_function

end

function e_boson_correlator(state::MixedDestabilizer, system, e_deformator)
    """this calculates the bosonic correlator between two sites i and j of the state on a torus.
    e_boson is defined as a vertex defect proliferated by Z bond-operators.
    state            : state to be examined.
    system           : system object that has knows everything!
    i, j             : indices of the sites to be correlated.
    Best is to precompute the deformator (invalid-)stabiliser state.
    """

    L = system.L
    k = system.k
    nbits = system.nbits
    cell_num = system.cell_num

    corr_function = zeros(Float64, L, L)

    for r_x  = 0:L-1
        for r_y = 0:L-1
            if r_x == 0 && r_y == 0
                corr_function[r_x+1, r_y+1] = 1
            else
                e_representative = get_e_reprentative((r_x, r_y), system)
                size_of_bff = 0
                dim_source, _ = zassenhausen_alg(e_representative, e_deformator, stabilizerview(state), nbits)
                dim_norm, _ = zassenhausen_alg(e_deformator, stabilizerview(state), nbits)
                corr_function[r_x+1, r_y+1] = dim_source - dim_norm
            end
        end
    end

    return corr_function

end