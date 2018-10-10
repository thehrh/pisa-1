# author: T. Ehrhardt
# date:   June 29, 2017
"""
OscParams: Characterize neutrino oscillation parameters
           (mixing angles, Dirac-type CP-violating phase, mass splittings)
"""

from __future__ import division

import numpy as np

from pisa import FTYPE


__all__ = ['OscParams']


class OscParams(object):
    """
    Holds neutrino oscillation parameters, i.e., mixing angles, squared-mass
    differences, and a Dirac-type CPV phase. The neutrino mixing (PMNS) matrix
    constructed from these parameters is given in the standard
    3x3 parameterization. Also holds the generalised matter potential matrix
    (divided by the matter potential a), i.e. diag(1, 0, 0) for the standard
    case.

    Parameters
    ----------
    dm21, dm31, dm41 : float
        Mass splittings (delta M^2_{21,31,41}) expected to be given in [eV^2]

    sin12, sin13, sin23 : float
        1-2, 1-3 and 2-3 mixing angles, interpreted as sin(theta_{ij})

    deltacp : float
        Value of CPV phase in [rad]


    Attributes
    ----------
    dm21, dm31, dm41 : float
        Cf. parameters

    sin12, sin13, sin23, sin14 : float
        Cf. parameters

    theta12, theta13, theta23, theta14 : float
        Mixing angles (corresponding to sinXY)

    deltacp : float
        Cf. parameters

    gen_mat_pot_matrix : 3d float array of shape (3, 3, 2)

    gen_mat_pot_matrix : 3d complex array

    mix_matrix : 3d float array of shape (3, 3, 2)
        Neutrino mixing (PMNS) matrix in standard parameterization. The third
        dimension holds the real and imaginary parts of each matrix element.

    mix_matrix_complex : 3d complex array

    mix_matrix_reparam : 3d float array of shape (3, 3, 2)
        Reparameterized neutrino mixing matrix, such that CPT invariance
        of vacuum propagation implemented by 3 simultaneous osc. param.
        transformations.

    mix_matrix_reparam_complex : 3d complex array

    dm_matrix : 2d float array of shape (3, 3)
        Antisymmetric matrix of squared-mass differences in vacuum

    """
    def __init__(self):

        self._sin12 = 0.
        self._sin13 = 0.
        self._sin23 = 0.
        self._sin14 = 0.
        self._deltacp = 0.
        self.dm21 = 0.
        self.dm31 = 0.
        self.dm41 = 0.
        self._eps_scale = 1.
        self._eps_prime = 0.
        self._sin_phi12 = 0.
        self._sin_phi13 = 0.
        self._sin_phi23 = 0.
        self._alpha1 = 0.
        self._alpha2 = 0.
        self._deltansi = 0.

    # --- theta12 ---
    @property
    def sin12(self):
        """Sine of 1-2 mixing angle"""
        return self._sin12

    @sin12.setter
    def sin12(self, value):
        assert (abs(value) <= 1)
        self._sin12 = value

    @property
    def theta12(self):
        return np.arcsin(self.sin12)

    @theta12.setter
    def theta12(self, value):
        self.sin12 = np.sin(value)

    # --- theta13 ---
    @property
    def sin13(self):
        """Sine of 1-3 mixing angle"""
        return self._sin13

    @sin13.setter
    def sin13(self, value):
        assert (abs(value) <= 1)
        self._sin13 = value

    @property
    def theta13(self):
        return np.arcsin(self.sin13)

    @theta13.setter
    def theta13(self, value):
        self.sin13 = np.sin(value)

    # --- theta23 ---
    @property
    def sin23(self):
        """Sine of 2-3 mixing angle"""
        return self._sin23

    @sin23.setter
    def sin23(self, value):
        assert (abs(value) <= 1)
        self._sin23 = value

    @property
    def theta23(self):
        return np.arcsin(self.sin23)

    @theta23.setter
    def theta23(self, value):
        self.sin23 = np.sin(value)

    # --- theta14 ---
    @property
    def sin14(self):
        """Sine of 1-4 mixing angle"""
        return self._sin14

    @sin14.setter
    def sin14(self, value):
        assert (abs(value) <= 1)
        self._sin14 = value

    @property
    def theta14(self):
        return np.arcsin(self.sin14)

    @theta14.setter
    def theta14(self, value):
        self.sin14 = np.sin(value)

    # --- deltaCP ---
    @property
    def deltacp(self):
        """CPV phase"""
        return self._deltacp

    @deltacp.setter
    def deltacp(self, value):
        assert value >= 0. and value <= 2*np.pi
        self._deltacp = value

    # --- NSI epsilons ---
    # --- overall matter potential strength ---
    @property
    def eps_scale(self):
        """Generalised matter potential strength scale"""
        return self._eps_scale

    @eps_scale.setter
    def eps_scale(self, value):
        self._eps_scale = value

    @property
    def eps_prime(self):
        """Second Hmat eigenvalue (beside eps_scale)"""
        return self._eps_prime

    @eps_prime.setter
    def eps_prime(self, value):
        self._eps_prime = value

    # --- projection phases ---
    # --- phi12 ---
    @property
    def sin_phi12(self):
        """1-2 phase"""
        return self._sin_phi12

    @sin_phi12.setter
    def sin_phi12(self, value):
        assert (abs(value) <= 1)
        self._sin_phi12 = value

    @property
    def phi12(self):
        return np.arcsin(self.sin_phi12)

    @phi12.setter
    def phi12(self, value):
        self.sin_phi12 = np.sin(value)

    # --- phi13 ---
    @property
    def sin_phi13(self):
        """1-3 phase"""
        return self._sin_phi13

    @sin_phi13.setter
    def sin_phi13(self, value):
        assert (abs(value) <= 1)
        self._sin_phi13 = value

    @property
    def phi13(self):
        return np.arcsin(self.sin_phi13)

    @phi13.setter
    def phi13(self, value):
        self.sin_phi13 = np.sin(value)

    # --- phi23 ---
    @property
    def sin_phi23(self):
        """2-3 phase"""
        return self._sin_phi23

    @sin_phi23.setter
    def sin_phi23(self, value):
        assert (abs(value) <= 1)
        self._sin_phi23 = value

    @property
    def phi23(self):
        return np.arcsin(self.sin_phi23)

    @phi23.setter
    def phi23(self, value):
        self.sin_phi23 = np.sin(value)

    # --- vacuum-matter relative phases ---
    # --- alpha1 ---
    @property
    def alpha1(self):
        """1-phase"""
        return self._alpha1

    @alpha1.setter
    def alpha1(self, value):
        assert value >= 0. and value <= 2*np.pi
        self._alpha1 = value

    # --- alpha2 ---
    @property
    def alpha2(self):
        """2-phase"""
        return self._alpha2

    @alpha2.setter
    def alpha2(self, value):
        assert value >= 0. and value <= 2*np.pi
        self._alpha2 = value

    # --- nsi phase ---
    @property
    def deltansi(self):
        """NSI phase"""
        return self._deltansi

    @deltansi.setter
    def deltansi(self, value):
        assert value >= 0. and value <= 2*np.pi
        self._deltansi = value

    # --- nsi coupling matrix elements ---
    @property
    def eps_ee(self):
        """nue-nue NSI coupling parameter"""
        return self.gen_mat_pot_matrix_complex[0, 0].real  - 1

    @property
    def eps_emu(self):
        """nue-numu NSI coupling parameter"""
        return self.gen_mat_pot_matrix_complex[0, 1]

    @property
    def eps_etau(self):
        """nue-nutau NSI coupling parameter"""
        return self.gen_mat_pot_matrix_complex[0, 2]

    @property
    def eps_mumu(self):
        """numu-numu NSI coupling parameter"""
        return self.gen_mat_pot_matrix_complex[1, 1].real

    @property
    def eps_mutau(self):
        """numu-nutau NSI coupling parameter"""
        return self.gen_mat_pot_matrix_complex[1, 2]

    @property
    def eps_tautau(self):
        """nutau-nutau NSI coupling parameter"""
        return self.gen_mat_pot_matrix_complex[2, 2].real

    @property
    def gen_mat_pot_matrix(self):
        """Matter Hamiltonian without the matter parameter a=sqrt(2)G_F N_e"""
        pot = np.zeros((3, 3, 2), dtype=FTYPE)

        sp12 = self.sin_phi12
        sp13 = self.sin_phi13
        sp23 = self.sin_phi23
        cp12 = np.sqrt(1. - sp12**2)
        cp13 = np.sqrt(1. - sp13**2)
        cp23 = np.sqrt(1. - sp23**2)

        sdnsi = np.sin(self.deltansi)
        cdnsi = np.cos(self.deltansi)

        # 1 + eps_ee - eps_mumu (real)
        pot[0, 0, 0] = (
            self.eps_scale * cp13**2 * (cp12**2 - sp12**2) +
            self.eps_prime * (
                (cp12**2 - sp12**2) * (sp13**2 * sp23**2 - cp23**2) -
                4 * cp12 * sp12 * sp13 * cp23 * sp23 * cdnsi
            )
        )
        pot[0, 0, 1] = 0.
        # eps_emu (complex)
        pot[0, 1, 0] = (
            self.eps_scale * cp12 * sp12 * cp13**2 * np.cos(self.alpha1 - self.alpha2) +
            self.eps_prime * (
                (
                    cp12 * sp12 * (sp13**2 * sp23**2 - cp23**2) +
                    sp13 * cp23 * sp23 * cdnsi * (cp12**2 - sp12**2)
                ) * np.cos(self.alpha1 - self.alpha2) -
                (
                    sp13 * cp23 * sp23 * sdnsi
                ) * np.sin(self.alpha1 - self.alpha2)
            )
        )
        pot[0, 1, 1] = (
            self.eps_scale * cp12 * sp12 * cp13**2 * np.sin(self.alpha1 - self.alpha2) +
            self.eps_prime * (
                (
                    cp12 * sp12 * (sp13**2 * sp23**2 - cp23**2) +
                    sp13 * cp23* sp23 * cdnsi * (cp12**2 - sp12**2)
                ) * np.sin(self.alpha1 - self.alpha2) +
                (
                    sp13 * cp23 * sp23 * sdnsi
                ) * np.cos(self.alpha1 - self.alpha2)
            )
        )
        # eps_etau (complex)
        pot[0, 2, 0] = (
            -self.eps_scale * cp12 * sp13 * cp13 * np.cos(2 * self.alpha1 + self.alpha2) +
            self.eps_prime * (
                (
                    cp13 * sp23 * (cp12 * sp13 * sp23 - sp12 * cp23 * cdnsi)
                ) * np.cos(2 * self.alpha1  + self.alpha2) -
                (
                    cp13 * sp12 * cp23 * sp23 * sdnsi
                ) * np.sin(2 * self.alpha1 + self.alpha2)
            )
        )
        pot[0, 2, 1] = (
            -self.eps_scale * cp12* sp13 * cp13 * np.sin(2 * self.alpha1 + self.alpha2) +
            self.eps_prime * (
                (
                    cp13 * sp23 * (cp12 * sp13 * sp23 - sp12 * cp23 * cdnsi)
                ) * np.sin(2 * self.alpha1 + self.alpha2) +
                (
                    cp13 * sp23 * sp12 * cp23 * sdnsi
                ) * np.cos(2 * self.alpha1 + self.alpha2)
            )
        )
        # eps_emu* (complex)
        pot[1, 0, 0] = pot[0, 1, 0]
        pot[1, 0, 1] = -pot[0, 1, 1]
        # eps_etau* (complex)
        pot[2, 0, 0] = pot[0, 2, 0]
        pot[2, 0, 1] = -pot[0, 2, 1]
        # eps_mumu - eps_mumu (0 by definition)
        pot[1, 1, 0] = 0.
        pot[1, 1, 1] = 0.
        # eps_mutau (complex)
        pot[1, 2, 0] = (
            -self.eps_scale * sp12 * cp13 * sp13 * np.cos(self.alpha1 + 2 * self.alpha2) +
            self.eps_prime * (
                (
                    cp13 * sp23 * (sp12 * sp13 * sp23 + cp12 * cp23 * cdnsi)
                ) * np.cos(self.alpha1 + 2 * self.alpha2) +
                (
                    cp12 * cp13 * cp23 * sp23 * sdnsi
                ) * np.sin(self.alpha1 + 2 * self.alpha2)
            )
        )
        pot[1, 2, 1] = (
            -self.eps_scale * sp12 * cp13 * sp13 * np.sin(self.alpha1 + 2 * self.alpha2) +
            self.eps_prime * (
                (
                    -cp12 * cp13 * cp23 * sp23 * sdnsi
                ) * np.cos(self.alpha1 + 2 * self.alpha2) +
                (
                    cp13 * sp23 * (sp12 * sp13 * sp23 + cp12 * cp23 * cdnsi)
                ) * np.sin(self.alpha1 + 2 * self.alpha2)
            )
        )
        # eps_mutau* (complex)
        pot[2, 1, 0] = pot[1, 2, 0]
        pot[2, 1, 1] = -pot[1, 2, 1]
        # eps_tautau - eps_mumu (real)
        pot[2, 2, 0] = (
            self.eps_scale * (sp13**2 - cp13**2 * sp12**2) +
            self.eps_prime *(
                sp23**2 * (cp13**2 - sp12**2 * sp13**2) -
                2 * cp12 * sp12 * sp13 * cp23 * sp23 * cdnsi -
                cp12**2 * cp23**2
            )
        )
        pot[2, 2, 1] = 0.

        return pot

    @property
    def gen_mat_pot_matrix_complex(self):
        """General matter potential matrix as complex 2-d array"""
        pot = self.gen_mat_pot_matrix
        pot_complex = pot[:, :, 0] + pot[:, :, 1] * 1.j

        # make sure this is a valid Hermitian potential matrix
        # before returning anything
        assert np.array_equal(pot_complex, pot_complex.conj().T)

        return pot_complex


    @property
    def mix_matrix(self):
        """Neutrino mixing matrix in its 'standard' form"""
        mix = np.zeros((3, 3, 2), dtype=FTYPE)

        sd = np.sin(self.deltacp)
        cd = np.cos(self.deltacp)

        c12 = np.sqrt(1. - self.sin12**2)
        c23 = np.sqrt(1. - self.sin23**2)
        c13 = np.sqrt(1. - self.sin13**2)

        mix[0, 0, 0] = c12 * c13
        mix[0, 0, 1] = 0.
        mix[0, 1, 0] = self.sin12 * c13
        mix[0, 1, 1] = 0.
        mix[0, 2, 0] = self.sin13 * cd
        mix[0, 2, 1] = - self.sin13 * sd
        mix[1, 0, 0] = - self.sin12 * c23 - c12 * self.sin23 * self.sin13 * cd
        mix[1, 0, 1] = - c12 * self.sin23 * self.sin13 * sd
        mix[1, 1, 0] = c12 * c23 - self.sin12 * self.sin23 * self.sin13 * cd
        mix[1, 1, 1] = - self.sin12 * self.sin23 * self.sin13 * sd
        mix[1, 2, 0] = self.sin23 * c13
        mix[1, 2, 1] = 0.
        mix[2, 0, 0] = self.sin12 * self.sin23 - c12 * c23 * self.sin13 * cd
        mix[2, 0, 1] = - c12 * c23 * self.sin13 * sd
        mix[2, 1, 0] = - c12 * self.sin23 - self.sin12 * c23 * self.sin13 * cd
        mix[2, 1, 1] = - self.sin12 * c23 * self.sin13 * sd
        mix[2, 2, 0] = c23 * c13
        mix[2, 2, 1] = 0.

        return mix

    @property
    def mix_matrix_complex(self):
        """Mixing matrix as complex 2-d array"""
        return self.mix_matrix[:, :, 0] + self.mix_matrix[:, :, 1] * 1.j

    @property
    def mix_matrix_reparam(self):
        """
        Neutrino mixing matrix reparameterised in a way
        such that the CPT trafo Hvac -> -Hvac*  is exactly implemented by
        the simultaneous transformations
            * deltamsq31 -> -deltamsq32
            * theta12 -> pi/2 - theta12
            * deltacp -> pi - deltacp

        which hence leave vacuum propagation invariant.

        This representation follows from the standard form U
        as diag(exp(i*deltacp), 0, 0) * U * diag(exp(-i*deltacp), 0, 0).

        """
        mix = np.zeros((3, 3, 2), dtype=FTYPE)

        sd = np.sin(self.deltacp)
        cd = np.cos(self.deltacp)

        c12 = np.sqrt(1. - self.sin12**2)
        c23 = np.sqrt(1. - self.sin23**2)
        c13 = np.sqrt(1. - self.sin13**2)

        mix[0, 0, 0] = c12 * c13
        mix[0, 0, 1] = 0.
        mix[0, 1, 0] = self.sin12 * c13 * cd
        mix[0, 1, 1] = self.sin12 * c13 * sd
        mix[0, 2, 0] = self.sin13
        mix[0, 2, 1] = 0.
        mix[1, 0, 0] = - self.sin12 * c23 * cd - c12 * self.sin23 * self.sin13
        mix[1, 0, 1] = self.sin12 * c23 * sd
        mix[1, 1, 0] = c12 * c23 - self.sin12 * self.sin23 * self.sin13 * cd
        mix[1, 1, 1] = - self.sin12 * self.sin23 * self.sin13 * sd
        mix[1, 2, 0] = self.sin23 * c13
        mix[1, 2, 1] = 0.
        mix[2, 0, 0] = self.sin12 * self.sin23 * cd - c12 * c23 * self.sin13
        mix[2, 0, 1] = - self.sin12 * self.sin23 * sd
        mix[2, 1, 0] = - c12 * self.sin23 - self.sin12 * c23 * self.sin13 * cd
        mix[2, 1, 1] = - self.sin12 * c23 * self.sin13 * sd
        mix[2, 2, 0] = c23 * c13
        mix[2, 2, 1] = 0.

        return mix

    @property
    def mix_matrix_reparam_complex(self):
        """Reparameterised mixing matrix as complex 2-d array"""
        return (self.mix_matrix_reparam[:, :, 0]
                + self.mix_matrix_reparam[:, :, 1] * 1.j)

    @property
    def dm_matrix(self):
        """Neutrino mass splitting matrix in vacuum"""
        dmVacVac = np.zeros((3, 3), dtype=FTYPE)
        mVac = np.zeros(3, dtype=FTYPE)
        delta = 5.e-9

        mVac[0] = 0.
        mVac[1] = self.dm21
        mVac[2] = self.dm31

        # Break any degeneracies
        if mVac[1] == 0.:
            mVac[0] -= delta
        if mVac[2] == 0.:
            mVac[2] += delta

        dmVacVac[0, 0] = 0.
        dmVacVac[1, 1] = 0.
        dmVacVac[2, 2] = 0.
        dmVacVac[0, 1] = mVac[0] - mVac[1]
        dmVacVac[1, 0] = - dmVacVac[0, 1]
        dmVacVac[0, 2] = mVac[0] - mVac[2]
        dmVacVac[2, 0] = - dmVacVac[0, 2]
        dmVacVac[1, 2] = mVac[1] - mVac[2]
        dmVacVac[2, 1] = - dmVacVac[1, 2]

        return dmVacVac

def test_nsi_parameterization():
    alpha1, alpha2, deltansi = np.random.rand(3) * 2. * np.pi
    phi12, phi13, phi23 = np.random.rand(3) * np.pi/2.
    eps_scale, eps_prime = np.random.rand(2) * 10.
    osc_params = OscParams()
    osc_params.eps_scale = eps_scale
    osc_params.eps_prime = eps_prime
    osc_params.phi12 = phi12
    osc_params.phi13 = phi13
    osc_params.phi23 = phi23
    osc_params.alpha1 = alpha1
    osc_params.alpha2 = alpha2
    osc_params.deltansi = deltansi
    # relative matter-nsi phases
    Qrel = (
        np.array([
            complex(np.cos(alpha1), np.sin(alpha1)),
            complex(np.cos(alpha2), np.sin(alpha2)),
            complex(np.cos(-(alpha1 + alpha2)), np.sin(-(alpha1 + alpha2)))
        ]) * np.eye(3, dtype=FTYPE)
    )
    # rotation matrices in right-handed convention
    R12 = np.array(
        [[np.cos(phi12), -np.sin(phi12), 0],
        [np.sin(phi12), np.cos(phi12), 0],
        [0, 0, 1]],
        dtype=FTYPE
    )
    R13 = np.array(
        [[np.cos(phi13), 0, np.sin(phi13)],
        [0, 1, 0],
        [-np.sin(phi13), 0, np.cos(phi13)]],
        dtype=FTYPE
    )
    R23_complex = np.array(
        [[1, 0, 0],
        [0, np.cos(phi23), -np.sin(phi23) * complex(np.cos(deltansi), np.sin(-deltansi))],
        [0, np.sin(phi23) * complex(np.cos(deltansi), np.sin(deltansi)), np.cos(phi23)]],
    )
    # "matter mixing matrix"
    Umat = np.matmul(R12, np.matmul(R13, R23_complex))
    # Hmat eigenvalues
    Dmat = np.array([eps_scale, eps_prime, 0], dtype=FTYPE) * np.eye(3, dtype=FTYPE)
    # matter Hamiltonian from matrix multiplication vs. analytically
    # start from the innermost product, work your way outwards
    Hmat_ref = np.matmul(
        Qrel,
        np.matmul(Umat,
            np.matmul(Dmat,
                np.matmul(Umat.conj().T, Qrel.conj().T)
            )
        )
    )
    # subtract mumu entry from diagonal entries (trace irrelevant)
    Hmat_ref = Hmat_ref - Hmat_ref[1, 1] * np.eye(3, dtype=FTYPE)
    # already subtracted for class attribute
    Hmat = osc_params.gen_mat_pot_matrix_complex
    logging.info("Matter Hamiltonian from matrix multiplication:\n%s" % Hmat_ref)
    logging.info("Analytical expansion:\n%s" % Hmat)
    if not np.all(np.isclose(Hmat, Hmat_ref)):
        raise ValueError(
            'Evaluating analytical expressions for matter Hamiltonian elements'
            ' does not give agreement with numerical calculation!'
        )

def test_sympy_mat_mult():
    """
    Sympy calculation of generalised matter Hamiltonian.
    Mainly for reference.

    """
    from sympy import (cos, sin, Matrix, eye, I, Symbol, symbols)
    from sympy.physics.quantum.dagger import Dagger
    phi12, phi13, phi23 = symbols('phi12 phi13 phi23', real=True)
    alpha1, alpha2 = symbols('alpha1 alpha2', real=True)
    eps_scale, eps_prime = symbols('eps_scale eps_prime', real=True)
    deltansi = Symbol('deltansi', real=True)

    Dmat = Matrix(
        [[eps_scale, 0, 0], [0, eps_prime, 0], [0, 0, 0]]
    )
    Qrel = Matrix(
        [[cos(alpha1) + I * sin(alpha1), 0, 0],
        [0, cos(alpha2) + I * sin(alpha2), 0],
        [0, 0, cos(-(alpha1 + alpha2)) + I * sin(-(alpha1 + alpha2))]]
    )
    R12 = Matrix(
        [[cos(phi12), -sin(phi12), 0],
        [sin(phi12), cos(phi12), 0],
        [0, 0, 1]]
    )
    R13 = Matrix(
        [[cos(phi13), 0, sin(phi13)],
        [0, 1, 0],
        [-sin(phi13), 0, cos(phi13)]]
    )
    R23_complex = Matrix(
        [[1, 0, 0],
        [0, cos(phi23), -sin(phi23) * (cos(deltansi) + I * sin(-deltansi))],
        [0, sin(phi23) * (cos(deltansi) + I * sin(deltansi)), cos(phi23)]]
    )

    Umat = R12 * R13 * R23_complex
    tmp = Dagger(Umat) * Dagger(Qrel)
    tmp2 = Dmat * tmp
    tmp3 = Umat * tmp2
    Hmat_sympy = Qrel * tmp3
    Hmat_sympy_minus_mumu = Hmat_sympy - Hmat_sympy[1, 1] * eye(3)
    return Hmat_sympy_minus_mumu


if __name__=='__main__':
    from pisa import TARGET
    from pisa.utils.log import set_verbosity, logging
    assert TARGET == 'cpu', "Cannot test functions on GPU, set PISA_TARGET to 'cpu'"
    set_verbosity(1)
    test_nsi_parameterization()
    try:
        test_sympy_mat_mult()
    except:
        pass
