#!/usr/bin/env python
"""Executable file for generating and solving diafiltration model."""

from diafiltration_flowsheet_model import diafiltration_model
from idaes.core.util.model_statistics import report_statistics
import pyomo.environ as pyo
import utils
import sys


def main(args):
    """Driver for creating diafiltration model."""
    # collect arguments
    args = args[1:]
    mix_style = args[0]
    num_s = int(args[1])
    num_t = int(args[2])

    # set relevant parameter values
    NS = num_s
    NT = num_t
    solutes = ['Li', 'Co']
    flux = 0.1
    sieving_coefficient = {"Li": 1.3, "Co": 0.5}
    feed = {
        'solvent': 100,    # m^3/hr of water
        'Li': 1.7*100,     # kg/hr
        'Co': 17*100,     # kg/hr
    }
    diaf = {
        'solvent': 30,     # m^3/hr of water
        'Li': 0.1*30,      # kg/hr
        'Co': 0.2*30,      # kg/hr
    }
    prec = True

    df = diafiltration_model(
        NS=NS,
        NT=NT,
        solutes=solutes,
        flux=flux,
        sieving_coefficient=sieving_coefficient,
        feed=feed,
        diaf=diaf,
        precipitate=prec,
        perc_precipitate={'permeate': {'Li': 0.81, 'Co': 0.01},
                          'retentate': {'Li': 0.20, 'Co': 0.89}}
    )

    mixing = mix_style
    m = df.build_flowsheet(mixing=mixing)
    df.initialize(m, mixing=mixing, precipitate=prec)
    df.unfix_dof(m, mixing=mixing, precipitate=prec)
    m.fs.precipitator['retentate'].V.fix(500)
    m.fs.precipitator['permeate'].V.fix(500)
    report_statistics(m)

    # solve model
    # R is used for the Li LB constraint.
    # This can be changed to any desired LB.
    m.R = 0.8
    solver = pyo.SolverFactory('ipopt')
    result = solver.solve(m, tee=True)

    # NOTE The below 2 percent recoveries are for use without precipitators
    # m.rec_perc_co.display()
    # m.rec_perc_li.display()
    # NOTE These percent recoveries are for precipitators
    m.prec_perc_co.display()
    m.prec_perc_li.display()

    # Print all relevant flow information
    utils.report_values(m)


if __name__ == '__main__':
    main(sys.argv)
