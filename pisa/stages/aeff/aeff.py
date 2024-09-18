"""
PISA pi stage to apply effective area weights
"""

from __future__ import absolute_import, print_function, division

from pisa.core.stage import Stage
from pisa.utils.profiler import profile


class aeff(Stage):  # pylint: disable=invalid-name
    """
    PISA Pi stage to apply aeff weights.

    This combines the detector effective area with the flux weights calculated
    in an earlier stage to compute the weights.

    Various scalings can be applied for particular event classes. The weight is
    then multiplied by the livetime to get an event count.

    Parameters
    ----------
    params
        Expected params are .. ::

            livetime : Quantity with time units
            aeff_scale : dimensionless Quantity
            nu*_cc_norm : dimensionless Quantity
            nu*_nc_norm : dimensionelss Quantity

    """
    def __init__(
        self,
        **std_kwargs,
    ):
        expected_params = (
            'livetime',
            'aeff_scale',
            'nue_cc_norm',
            'nuebar_cc_norm',
            'numu_cc_norm',
            'numubar_cc_norm',
            'nutaubar_cc_norm',
            'nutau_cc_norm',
            'nunubar_nc_norm',
        )

        # init base class
        super().__init__(
            expected_params=expected_params,
            **std_kwargs,
        )


    @profile
    def apply_function(self):

        # read out
        aeff_scale = self.params.aeff_scale.m_as('dimensionless')
        livetime_s = self.params.livetime.m_as('sec')
        nue_cc_norm = self.params.nue_cc_norm.m_as('dimensionless')
        nuebar_cc_norm = self.params.nuebar_cc_norm.m_as('dimensionless')
        numu_cc_norm = self.params.numu_cc_norm.m_as('dimensionless')
        numubar_cc_norm = self.params.numubar_cc_norm.m_as('dimensionless')
        nutau_cc_norm = self.params.nutau_cc_norm.m_as('dimensionless')
        nutaubar_cc_norm = self.params.nutaubar_cc_norm.m_as('dimensionless')
        nunubar_nc_norm = self.params.nunubar_nc_norm.m_as('dimensionless')

        for container in self.data:
            scale = aeff_scale * livetime_s
            if container.name == 'nue_cc':
                scale *= nue_cc_norm
            elif container.name == 'nuebar_cc':
                scale *= nuebar_cc_norm
            elif container.name == 'numu_cc':
                scale *= numu_cc_norm
            elif container.name == 'numubar_cc':
                scale *= numubar_cc_norm
            elif container.name == 'nutau_cc':
                scale *= nutau_cc_norm
            elif container.name == 'nutaubar_cc':
                scale *= nutaubar_cc_norm
            elif 'nc' in container.name:
                scale *= nunubar_nc_norm

            container['weights'] *= container['weighted_aeff'] * scale
            container.mark_changed('weights')
