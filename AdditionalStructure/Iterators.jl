# Iterators and other

function iterate_measurements_only!(state::MixedDestabilizer, system, get_random_pauli::Function, full_steps::Integer)
    """ Evolves the state by a set number of full steps (intensive variable).
    state - the refernece to a state involved, it is mutable and will change
    system - details of the geometry
    get_random_pauli - defines the model
    full_steps - the number of steps in units of system size
    """

    nbits = system.nbits

    for i_n = 1:nbits
        for i_s = 1:full_steps
            projectrand!(state, get_random_pauli())
        end
    end
    
    return state
end