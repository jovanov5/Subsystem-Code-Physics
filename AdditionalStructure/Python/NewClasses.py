import itertools

class EdgeSquareLattice:
    """
    The class that represents a square lattice on a torus with qudits places at the edges.
    Hold geometrical data in terms of a couple of integers (ToDo: add regions ranges and projectors, for TEE and TEN).
    """

    def __init__(self, L):
        """
        The simplest constructor of the class.
        :param L: The linear size of the lattice.
        """
        self.L = L
        self.k = 2 # Assumed to always be true, don't know if to generalise at the moment, maybe in the future.
        self.cell_num = L * L
        self.nbits = self.k * self.cell_num


class ToricCode:
    """
    Toric Code model defined on a EddgeSquareLattice.

    The stabilisers are represented as one-dim vectors of length 2nbits, starting with X operators and then Z operators.
    """
    
    def __init__(self, L):
        """
        The simplest constructor of the class.
        :param L: The linear size of the lattice.
        """
        self.system = EdgeSquareLattice(L)
        self.plaquettes = [self.get_plaquette(i, j) for (i, j) in itertools.product(range(self.system.L), range(self.system.L))] # Sage will like this definition, these are the check rows 1x(2nbits)
        self.stars = [self.get_star(i, j) for (i, j) in itertools.product(range(self.system.L), range(self.system.L))] # Sage will like this definition, these are the check rows 1x(2nbits)
        self.stabilisers = self.plaquettes + self.stars # Sage will like this definition, these are the check rows 1x(2nbits)
        self.electrons = [self.get_electron(i, j, orientation) for (i, j, orientation) in itertools.product(range(self.system.L), range(self.system.L), range(2))] # Sage will like this definition, these are the check rows 1x(2nbits)
        self.fluxes = [self.get_flux(i, j, orientation) for (i, j, orientation) in itertools.product(range(self.system.L), range(self.system.L), range(2))] # Sage will like this definition, these are the check rows 1x(2nbits)
        self.fermions = [self.get_fermion(i, j, orientation) for (i, j, orientation) in itertools.product(range(self.system.L), range(self.system.L), range(2))] # Sage will like this definition, these are the check rows 1x(2nbits)
    
    def get_plaquette(self, i, j):
        """
        Get the Z-plaquette operator at the given position.
        :param i: The i-th row.
        :param j: The j-th column.
        :return: The plaquette operator as a one-dim vector of length 2nbits.
        """
        cell_lin_index = i * self.system.L + j
        cell_up_lin_index = i * self.system.L + (j + 1) % self.system.L
        cell_right_lin_index = (i + 1) % self.system.L * self.system.L + j
        z_positions = [2*cell_lin_index,
                       2*cell_lin_index + 1,
                       2*cell_right_lin_index + 1,
                       2*cell_up_lin_index] # The Z operators in the plaquette. Edges come 2 per cell.

        operator = [0] * (2 * self.system.nbits) # for both X and Z operators. Plaq uses Z.
        for pos in z_positions:
            operator[self.system.nbits + pos] = 1 # X then Z, remember.
        return operator
    
    def get_star(self, i, j):
        """
        Get the X-star operator at the given position.
        :param i: The i-th row.
        :param j: The j-th column.
        :return: The star operator as a one-dim vector of length 2nbits.
        """
        cell_lin_index = i * self.system.L + j
        cell_down_lin_index = i * self.system.L + (j - 1) % self.system.L
        cell_left_lin_index = (i - 1) % self.system.L * self.system.L + j
        x_positions = [2*cell_lin_index,
                       2*cell_lin_index + 1,
                       2*cell_left_lin_index,
                       2*cell_down_lin_index + 1]
        operator = [0] * (2 * self.system.nbits)
        for pos in x_positions:
            operator[pos] = 1 # X then Z, remember.
        return operator
    
    def get_stab(self, i, j, kind):
        """
        Get the X-star or Z-plaq operator at the given position.
        :param i: The i-th row.
        :param j: The j-th column.
        :return: The star operator as a one-dim vector of length 2nbits.
        """
        cell_lin_index = i * self.system.L + j
        cell_down_lin_index = i * self.system.L + (j - 1) % self.system.L
        cell_left_lin_index = (i - 1) % self.system.L * self.system.L + j
        cell_up_lin_index = i * self.system.L + (j + 1) % self.system.L
        cell_right_lin_index = (i + 1) % self.system.L * self.system.L + j

        if kind == 0:
            x_positions = [2*cell_lin_index,
                       2*cell_lin_index + 1,
                       2*cell_left_lin_index,
                       2*cell_down_lin_index + 1]
            z_positions = []
        else:
            x_positions = []
            z_positions = [2*cell_lin_index,
                       2*cell_lin_index + 1,
                       2*cell_right_lin_index + 1,
                       2*cell_up_lin_index] # The Z operators in the plaquette. Edges come 2 per cell.

        operator = [0] * (2 * self.system.nbits)
        for pos in x_positions:
            operator[pos] = 1 # X then Z, remember.
        for pos in z_positions:
            operator[self.system.nbits + pos] = 1
        return operator
    
    def get_electron(self, i, j, orientation):
        """
        Get the Z-electon string operator at the given position.
        :param i: The i-th row.
        :param j: The j-th column.
        :param orientation: The orientation of the electron string.
        :return: The star operator as a one-dim vector of length 2nbits.
        """
        cell_lin_index = i * self.system.L + j
        cell_up_lin_index = i * self.system.L + (j + 1) % self.system.L
        cell_right_lin_index = (i + 1) % self.system.L * self.system.L + j
        cell_down_lin_index = i * self.system.L + (j - 1) % self.system.L
        cell_left_lin_index = (i - 1) % self.system.L * self.system.L + j
        if orientation == 0:
            x_positions = []
            z_positions = [2*cell_left_lin_index]
        else:
            x_positions = []
            z_positions = [2*cell_down_lin_index + 1]
        
        operator = [0] * (2 * self.system.nbits)
        for pos in x_positions:
            operator[pos] = 1 # X then Z, remember.
        for pos in z_positions:
            operator[self.system.nbits + pos] = 1 # X then Z, remember.
        return operator
    
    def get_flux(self, i, j, orientation):
        """
        Get the X-flux string operator at the given position.
        :param i: The i-th row.
        :param j: The j-th column.
        :param orientation: The orientation of the electron string.
        :return: The star operator as a one-dim vector of length 2nbits.
        """
        cell_lin_index = i * self.system.L + j
        cell_up_lin_index = i * self.system.L + (j + 1) % self.system.L
        cell_right_lin_index = (i + 1) % self.system.L * self.system.L + j
        cell_down_lin_index = i * self.system.L + (j - 1) % self.system.L
        cell_left_lin_index = (i - 1) % self.system.L * self.system.L + j
        if orientation == 0:
            x_positions = [2*cell_lin_index + 1]
            z_positions = []
        else:
            x_positions = [2*cell_lin_index]
            z_positions = []
        
        operator = [0] * (2 * self.system.nbits)
        for pos in x_positions:
            operator[pos] = 1 # X then Z, remember.
        for pos in z_positions:
            operator[self.system.nbits + pos] = 1 # X then Z, remember.
        return operator
    
    def get_fermion(self, i, j, orientation):
        """
        Get the ZX-fermion string operator at the given position.
        :param i: The i-th row.
        :param j: The j-th column.
        :param orientation: The orientation of the electron string.
        :return: The star operator as a one-dim vector of length 2nbits.
        """
        cell_lin_index = i * self.system.L + j
        cell_up_lin_index = i * self.system.L + (j + 1) % self.system.L
        cell_right_lin_index = (i + 1) % self.system.L * self.system.L + j
        cell_down_lin_index = i * self.system.L + (j - 1) % self.system.L
        cell_left_lin_index = (i - 1) % self.system.L * self.system.L + j
        if orientation == 0:
            x_positions = [2*cell_lin_index + 1]
            z_positions = [2*cell_left_lin_index]
        else:
            x_positions = [2*cell_lin_index]
            z_positions = [2*cell_down_lin_index + 1]
        
        operator = [0] * (2 * self.system.nbits)
        for pos in x_positions:
            operator[pos] = 1 # X then Z, remember.
        for pos in z_positions:
            operator[self.system.nbits + pos] = 1 # X then Z, remember.
        return operator
