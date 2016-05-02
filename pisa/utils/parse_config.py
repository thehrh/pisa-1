# author: P.Eller
#         pde3+pisa@psu.edu
#
# date:   2016-04-28
"""
Parse a ConfigFile object into a dict containing a an entry for every
stages, that contains values indicated by param. as a param set, binning as
binning objects and all other values as ordinary strings
"""

from collections import OrderedDict

import numpy as np
import uncertainties
from uncertainties import unumpy as unp
from uncertainties import ufloat, ufloat_fromstr
import pint
ureg = pint.UnitRegistry()
# Config files use "uinits.xyz" to denote that "xyz" is a unit; therefore,
# ureg is also referred to as "units" in this context.
units = ureg

from pisa.core.prior import Prior
from pisa.core.param import Param, ParamSet
from pisa.utils.log import logging
from pisa.utils.fileio import from_file
from pisa.core.binning import OneDimBinning, MultiDimBinning


def parse_quantity(string):
    value = string.replace(' ', '')
    if 'units.' in value:
        value, unit = value.split('units.')
    else:
        unit = None
    value = value.rstrip('*')
    if '+/-' in value:
        value = ufloat_fromstr(value)
    else:
        value = ufloat(float(value), 0)
    value *= ureg(unit)
    return value


def parse_string_literal(string):
    return string


def list_split(string):
    list = string.split(',')
    return [x.strip() for x in list]


def parse_config(config):
    if isinstance(config, basestring):
        config = from_file(config)
    # create binning objects
    binning_dict = {}
    order = list_split(config.get('binning', 'order'))
    binnings = list_split(config.get('binning', 'binnings'))
    for binning in binnings:
        bins = []
        for bin_name in order:
            args = eval(config.get('binning', binning + '.' + bin_name))
            bins.append(OneDimBinning(bin_name, **args))
        binning_dict[binning] = MultiDimBinning(*bins)

    args_dict = OrderedDict()
    # find pipline setting
    pipeline_order = list_split(config.get('pipeline', 'order'))
    for item in pipeline_order:
        stage, service = item.split(':')
        section = 'stage:' + stage
        # get infos for stages
        args_dict[stage] = {}
        args_dict[stage]['service'] = service
        params = []
        if config.has_option(section, 'param_selector'):
            param_selector = config.get(section, 'param_selector')
        else:
            param_selector = ''

        for name, value in config.items(section):
            if name.startswith('param.'):
                # find parameter root
                if name.startswith('param.'+ param_selector + '.') and \
                        name.count('.') == 2:
                    _, _, pname = name.split('.')
                elif name.startswith('param.') and name.count('.') == 1:
                    _, pname = name.split('.')
                else:
                    continue

                # defaults
                args = {'name': pname, 'is_fixed': True, 'prior': None,
                        'range': None}
                try:
                    value = parse_quantity(value)
                    args['value'] = value.n * value.units
                except ValueError:
                    value = parse_string_literal(value)
                    args['value'] = value
                # search for explicit specifications
                if config.has_option(section, name + '.fixed'):
                    args['is_fixed'] = config.getboolean(section, name +
                                                         '.fixed')
                if config.has_option(section, name + '.prior'):
                    if config.get(section, name + '.prior') == 'uniform':
                        args['prior'] = Prior(kind='uniform')
                    elif config.get(section, name + '.prior') == 'spline':
                        priorname = pname
                        if param_selector:
                            priorname += '_' + param_selector
                        data = config.get(section, name + '.prior.data')
                        data = from_file(data)
                        data = data[priorname]
                        knots = ureg.Quantity(np.asarray(data['knots']), data['units'])
                        knots = knots.to(value.units)
                        coeffs = np.asarray(data['coeffs'])
                        deg = data['deg']
                        args['prior'] = Prior(kind='spline', knots=knots.m,
                                coeffs=coeffs,
                                deg=deg)
                    elif 'gauss' in config.get(section, name + '.prior'):
                        raise Exception('''Please use new style +/- notation for
                            gaussian priors in config''')
                    else:
                        raise Exception('Prior type unknown')
                elif hasattr(value, 's') and value.s != 0:
                    args['prior'] = Prior(kind='gaussian', fiducial=value.n,
                                          sigma = value.s)
                if config.has_option(section, name + '.range'):
                    range = config.get(section, name + '.range')
                    if 'nominal' in range:
                        nominal = value.n * value.units
                    if 'sigma' in range:
                        sigma = value.s * value.units
                    range = range.replace('[', 'np.array([')
                    range = range.replace(']', '])')
                    args['range'] = eval(range).to(value.units).m

                params.append(Param(**args))

            elif 'binning' in name:
                args_dict[stage][name] = binning_dict[value]

            elif not name == 'param_selector':
                args_dict[stage][name] = value

        if len(params) > 0:
            args_dict[stage]['params'] = ParamSet(*params)

    return args_dict
