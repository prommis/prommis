#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Utility functions."""

# imports
import pyomo.environ as pyo

##############################################################################
# Model Information Utility functions
##############################################################################


def report_values(m, prec=True):
    """
    Report flow values for precipitator model.

    This is generic for both stage and tube mixing.
    """
    NS = len(m.fs.stages)
    NT = len(m.fs.tubes)
    data = {}
    data["size"] = [NS, NT]
    data["solutes"] = [sol for sol in m.fs.solutes]

    print("Membrane outlet recoveries")
    m.rec_perc_co.display()
    m.rec_perc_li.display()
    data["Membrane outlet recoveries"] = [
        pyo.value(m.rec_perc_co),
        pyo.value(m.rec_perc_li),
    ]
    print("Membrane Area")
    data["Membrane Area"] = []
    for i in m.fs.stages:
        m.fs.stage[i].length.pprint()
        data["Membrane Area"].append(m.fs.stage[i].length.value)
    print("feed flows")
    data["feed flows"] = []
    for i in pyo.RangeSet(len(m.fs.stages)):
        print(f"Stage {i}")
        for j in pyo.RangeSet(len(m.fs.tubes)):
            try:
                flow = pyo.value(
                    getattr(m.fs.split_feed, f"outlet_{(i-1)*NT+j}").flow_vol[0]
                )
                print(flow)
                data["feed flows"].append(flow)
            except:
                print("no feeds here")
                data["feed flows"].append(0)
    print("\ndiaf flows")
    data["diaf flows"] = []
    for i in pyo.RangeSet(len(m.fs.stages)):
        print(f"Stage {i}")
        for j in pyo.RangeSet(len(m.fs.tubes)):
            try:
                flow = pyo.value(
                    getattr(m.fs.split_diafiltrate, f"outlet_{(i-1)*NT+j}").flow_vol[0]
                )
                print(flow)
                data["diaf flows"].append(flow)
            except:
                print("no diafiltrates here")
                data["diaf flows"].append(0)

    print("\nrecycle flows")
    data["recycle flows"] = []
    for i in pyo.RangeSet(2, len(m.fs.stages)):
        print(f"Stage {i-1}")
        for j in pyo.RangeSet(len(m.fs.tubes)):
            try:
                flow = pyo.value(
                    getattr(m.fs.recycle_splitters[i], f"outlet_{j}").flow_vol[0]
                )
                print(flow)
                data["recycle flows"].append(flow)
            except:
                print("no recycles here")
                data["recycle flows"].append(0)

    if hasattr(m.fs, "splitters"):
        print("\ninlet splitter flows")
        data["inlet splitter flows"] = []
        for i in m.fs.stages:
            print(f"Stage {i}")
            for j in m.fs.tubes:
                flow = pyo.value(getattr(m.fs.splitters[i], f"outlet_{j}").flow_vol[0])
                print(flow)
                data["inlet splitter flows"].append(flow)

    print("\nCo product flows")
    data["Co product flows"] = []
    for i in pyo.RangeSet(len(m.fs.stages)):
        try:
            flow = pyo.value(m.fs.split_retentate[i].product.flow_vol[0])
            print(flow)
            data["Co product flows"].append(flow)
            for sol in m.fs.solutes:
                print(
                    f"    {pyo.value(m.fs.split_retentate[i].product.mass_solute[0, sol])}"
                )
                data["Co product flows"].append(
                    pyo.value(m.fs.split_retentate[i].product.mass_solute[0, sol])
                )
        except:
            print("no Co product here")
            data["Co product flows"].append(0)

    print("\nLi product flows")
    data["Li product flows"] = []
    for i in pyo.RangeSet(len(m.fs.stages)):
        try:
            flow = pyo.value(m.fs.split_permeate[i].product.flow_vol[0])
            print(flow)
            data["Li product flows"].append(flow)
            for sol in m.fs.solutes:
                print(
                    f"    {pyo.value(m.fs.split_permeate[i].product.mass_solute[0, sol])}"
                )
                data["Li product flows"].append(
                    pyo.value(m.fs.split_permeate[i].product.mass_solute[0, sol])
                )
        except:
            print("no Li product here")
            data["Li product flows"].append(0)

    if prec:
        print("\ndiafiltrate inlet")
        data["diafiltrate inlet"] = []
        print(pyo.value(m.fs.split_diafiltrate.inlet.flow_vol[0]))
        data["diafiltrate inlet"].append(
            pyo.value(m.fs.split_diafiltrate.inlet.flow_vol[0])
        )
        for sol in m.fs.solutes:
            print(f"    {pyo.value(m.fs.split_diafiltrate.inlet.mass_solute[0, sol])}")
            data["diafiltrate inlet"].append(
                pyo.value(m.fs.split_diafiltrate.inlet.mass_solute[0, sol])
            )

        print("\nprecipitate products")
        data["precipitate products retentate"] = []
        print(pyo.value(m.fs.precipitator["retentate"].solid.flow_vol[0]))
        data["precipitate products retentate"].append(
            pyo.value(m.fs.precipitator["retentate"].solid.flow_vol[0])
        )
        for sol in m.fs.solutes:
            print(
                f'    {pyo.value(m.fs.precipitator["retentate"].solid.mass_solute[0, sol])}'
            )
            data["precipitate products retentate"].append(
                pyo.value(m.fs.precipitator["retentate"].solid.mass_solute[0, sol])
            )
        print("^retentate vpermeate")
        data["precipitate products permeate"] = []
        print(pyo.value(m.fs.precipitator["permeate"].solid.flow_vol[0]))
        data["precipitate products permeate"].append(
            pyo.value(m.fs.precipitator["permeate"].solid.flow_vol[0])
        )
        for sol in m.fs.solutes:
            print(
                f'    {pyo.value(m.fs.precipitator["permeate"].solid.mass_solute[0, sol])}'
            )
            data["precipitate products permeate"].append(
                pyo.value(m.fs.precipitator["permeate"].solid.mass_solute[0, sol])
            )
        print("\nprecipitate recycle")
        data["precipitate recycle retentate"] = []
        print(pyo.value(m.fs.precipitator["retentate"].recycle.flow_vol[0]))
        data["precipitate recycle retentate"].append(
            pyo.value(m.fs.precipitator["retentate"].recycle.flow_vol[0])
        )
        for sol in m.fs.solutes:
            print(
                f'    {pyo.value(m.fs.precipitator["retentate"].recycle.mass_solute[0, sol])}'
            )
            data["precipitate recycle retentate"].append(
                pyo.value(m.fs.precipitator["retentate"].recycle.mass_solute[0, sol])
            )
        print("^retentate vpermeate")
        data["precipitate recycle permeate"] = []
        print(pyo.value(m.fs.precipitator["permeate"].recycle.flow_vol[0]))
        data["precipitate recycle permeate"].append(
            pyo.value(m.fs.precipitator["permeate"].recycle.flow_vol[0])
        )
        for sol in m.fs.solutes:
            print(
                f'    {pyo.value(m.fs.precipitator["permeate"].recycle.mass_solute[0, sol])}'
            )
            data["precipitate recycle permeate"].append(
                pyo.value(m.fs.precipitator["permeate"].recycle.mass_solute[0, sol])
            )

        print("\nwaste")
        data["waste"] = []
        print(pyo.value(m.fs.split_precipitate_recycle.waste.flow_vol[0]))
        data["waste"].append(
            pyo.value(m.fs.split_precipitate_recycle.waste.flow_vol[0])
        )
        for sol in m.fs.solutes:
            print(
                f"    {pyo.value(m.fs.split_precipitate_recycle.waste.mass_solute[0, sol])}"
            )
            data["waste"].append(
                pyo.value(m.fs.split_precipitate_recycle.waste.mass_solute[0, sol])
            )

        print("\nsplit fractions for precipitator")
        data["prec sf retentate"] = []
        for sol in m.fs.solutes:
            print(m.fs.precipitator["retentate"].yields[sol, "solid"].value)
            data["prec sf retentate"].append(
                m.fs.precipitator["retentate"].yields[sol, "solid"].value
            )
        print("^retentate vpermeate")
        data["prec sf permeate"] = []
        for sol in m.fs.solutes:
            print(m.fs.precipitator["permeate"].yields[sol, "solid"].value)
            data["prec sf permeate"].append(
                m.fs.precipitator["permeate"].yields[sol, "solid"].value
            )

        print("\nactual Co/Li Recovery")
        print(pyo.value(m.prec_perc_co))
        print(pyo.value(m.prec_perc_li))
        data["actual recovery"] = [pyo.value(m.prec_perc_co), pyo.value(m.prec_perc_li)]
        m.fs.precipitator["retentate"].V.pprint()
        m.fs.precipitator["permeate"].V.pprint()
        data["prec volumes"] = [
            m.fs.precipitator["retentate"].V.value,
            m.fs.precipitator["permeate"].V.value,
        ]
    return data
