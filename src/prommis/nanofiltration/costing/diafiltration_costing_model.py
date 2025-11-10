"""Preliminary work on costing model for membrane flowsheets."""
# imports
from idaes.core.util.model_statistics import report_statistics
import pyomo.environ as pyo
import model.utils as utils

class costing_model:
    """
    Costing for the diafiltration model.

    Builds on top of the diafiltration flowsheet.
    Requires:
    working diafiltration flowsheet

    Parameters:
    precitator capital/operating cost [$/m^3]
    membrane capital/operating cost [$/m^2]
    pump capital/operating cost [$/kW]
    """

    def __init__(
            self,
            model,
            f_prec=1000,   # CEPCI 2024 Jun and 2001 NETL Process Equipment Cost Estimation
            f_prec_op={'permeate': .22, 'retentate': 0.1},  # currently not using
            f_memb=66,   # based on source from https://doi.org/10.1016/j.ijggc.2019.03.018, 798.8/607.5*50=742 using CEPCI 2024 Jun and 2019
            f_memb_op=0.2,   # based on watertap is 0.15. Just assume every 5 years is 0.2
            f_pump=700,    # $/kW, this paper was published in 2024
            f_pump_op=742.  # use 560 $/yr kw, 798.8/693.1*560=742 using CEPCI 2024 Jun and 2018
    ):
        """Set up costing model with appropriate cost factors."""
        self.m = model
        self.f_prec = f_prec
        self.f_prec_op = f_prec_op
        self.f_memb = f_memb
        self.f_memb_op = f_memb_op
        self.f_pump = f_pump
        self.f_pump_op = f_pump_op

    def build_costing(self):
        """Build the costing model."""
        m = self.m
        m.costing = pyo.Block()
        self.costing_prec(m, self.f_prec, self.f_prec_op)
        self.costing_memb(m, self.f_memb, self.f_memb_op)
        self.costing_pump(m, self.f_pump, self.f_pump_op)
        self.costing_obj(m)

        return m

    def costing_obj(self, m):
        """Set up costing objectives."""
        # make sure both Li/Co lower bounds are active.
        # precipitator=True is assumed.
        m.prec_li_lb.activate()
        m.prec_co_lb.activate()

        # remove recovery objectives
        m.prec_co_obj.deactivate()
        m.prec_li_obj.deactivate()

        # add costing objective
        # costing objective variables
        m.costing.eps = pyo.Param(initialize=0.1)

        # costing objective
        m.costing.obj = pyo.Objective(
            expr=((m.costing.RR*m.costing.C_memb
                   + m.costing.C_pump_op)
                  # + sum(m.costing.C_prec_op[i] for i in m.precips)
                  + m.costing.eps*(
                      sum(m.costing.C_prec[i] for i in m.precips)
                      + m.costing.C_memb
                      + m.costing.C_pump)
                  ),
            sense=pyo.minimize
        )

    def costing_prec(self, m, f_prec, f_prec_op):
        """Set up precipitator costing."""
        # precipitator cost varibles, parameters
        m.precips = pyo.Set(initialize=['permeate', 'retentate'])
        m.costing.C_prec = pyo.Var(m.precips)
        m.costing.f_prec = pyo.Param(initialize=f_prec)

        # Cost as linear function of volume
        m.costing.prec_cost = pyo.Constraint(
            m.precips,
            expr={
                i: m.costing.C_prec[i] == m.costing.f_prec*m.fs.precipitator[i].V + 12000
                for i in m.precips
            }
        )

    def costing_memb(self, m, f_memb, f_memb_op):
        """Set up membrane costing."""
        # membrane cost variables/parameters
        m.costing.C_memb = pyo.Var()
        m.costing.f_memb = pyo.Param(initialize=f_memb)
        m.costing.RR = pyo.Param(initialize=f_memb_op)

        # membrane cost constraints
        m.costing.memb_cost = pyo.Constraint(
            expr=m.costing.C_memb == sum(m.costing.f_memb*m.fs.stage[i].length
                                         for i in m.fs.stages)
        )

    def costing_pump(self, m, f_pump, f_pump_op):
        """Set up pump costing."""
        # pump cost variables/parameters
        m.costing.C_pump = pyo.Var()
        m.costing.f_pump = pyo.Param(initialize=f_pump)
        m.costing.P_inst = pyo.Var()
        # pump cost is a linear function of pump power specification
        m.costing.pump_cost = pyo.Constraint(
            expr=(m.costing.C_pump == m.costing.f_pump
                  * m.costing.P_inst)
        )
        # pump operating costs variables/parameters
        m.costing.C_pump_op = pyo.Var()
        m.costing.f_pump_op = pyo.Param(initialize=f_pump_op)
        m.costing.P_op = pyo.Var()

        # Assume pump operating power directly related to diafiltrate flow rate.
        m.costing.pump_P_op = pyo.Constraint(
            # [m^3/hr][1 hr / 3600 s][898415 J/m^3][1 kJ/1000J]
            # assuming change in pressure from atmospheric to 145psi
            expr=m.costing.P_op == .24956*m.fs.split_diafiltrate.inlet.flow_vol[0]
        )

        # operating cost is a linear function of operating power requirements.
        m.costing.pump_cost_op = pyo.Constraint(
            expr=(m.costing.C_pump_op == m.costing.f_pump_op
                  * m.costing.P_op)
        )

        # operating power cannot be greater than pump power specification.
        m.costing.pump_inst_op = pyo.Constraint(
            expr=m.costing.P_inst >= m.costing.P_op
        )
