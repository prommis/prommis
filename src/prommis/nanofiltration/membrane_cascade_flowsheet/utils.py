#####################################################################################################
# “PrOMMiS” was produced under the DOE Process Optimization and Modeling for Minerals Sustainability
# (“PrOMMiS”) initiative, and is copyright (c) 2023-2025 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory, et al. All rights reserved.
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license information.
#####################################################################################################
"""Utility functions."""

# imports
import pyomo.environ as pyo

import matplotlib.patches as patches
from matplotlib import pyplot as plt
from matplotlib.colors import to_rgba


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
        data["Membrane Area"].append(pyo.value(m.fs.stage[i].length))

    print("feed flows")
    data["feed flows"] = []
    for i in m.fs.inlet_mixers:
        if not isinstance(i, tuple):
            print(f"Stage {i}")
            flow = pyo.value(m.fs.inlet_mixers[i].feed_state[0].flow_vol)
            print(flow)
            data["feed flows"].append(flow)
        else:
            if i[1] == 1:
                print(f"Stage {i[0]}")
            flow = pyo.value(m.fs.inlet_mixers[i].feed_state[0].flow_vol)
            print(flow)
            data["feed flows"].append(flow)

    print("\ndiaf flows")
    data["diaf flows"] = []
    for i in m.fs.inlet_mixers:
        if not isinstance(i, tuple):
            print(f"Stage {i}")
            flow = pyo.value(m.fs.inlet_mixers[i].diafiltrate_state[0].flow_vol)
            print(flow)
            data["diaf flows"].append(flow)
        else:
            if i[1] == 1:
                print(f"Stage {i[0]}")
            flow = pyo.value(m.fs.inlet_mixers[i].diafiltrate_state[0].flow_vol)
            print(flow)
            data["diaf flows"].append(flow)
    print("\nrecycle flows")
    data["recycle flows"] = []
    for i in m.fs.inlet_mixers:
        if not isinstance(i, tuple):
            print(f"Stage {i}")
            flow = pyo.value(m.fs.inlet_mixers[i].recycle_state[0].flow_vol)
            print(flow)
            data["recycle flows"].append(flow)
        else:
            if i[1] == 1:
                print(f"Stage {i[0]}")
            flow = pyo.value(m.fs.inlet_mixers[i].recycle_state[0].flow_vol)
            print(flow)
            data["recycle flows"].append(flow)

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
                    f"    {pyo.value(m.fs.split_retentate[i].product.flow_mass_solute[0, sol])}"
                )
                data["Co product flows"].append(
                    pyo.value(m.fs.split_retentate[i].product.flow_mass_solute[0, sol])
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
                    f"    {pyo.value(m.fs.split_permeate[i].product.flow_mass_solute[0, sol])}"
                )
                data["Li product flows"].append(
                    pyo.value(m.fs.split_permeate[i].product.flow_mass_solute[0, sol])
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
            print(
                f"    {pyo.value(m.fs.split_diafiltrate.inlet.flow_mass_solute[0, sol])}"
            )
            data["diafiltrate inlet"].append(
                pyo.value(m.fs.split_diafiltrate.inlet.flow_mass_solute[0, sol])
            )

        print("\nprecipitate products")
        data["precipitate products retentate"] = []
        print(pyo.value(m.fs.precipitator["retentate"].solid.flow_vol[0]))
        data["precipitate products retentate"].append(
            pyo.value(m.fs.precipitator["retentate"].solid.flow_vol[0])
        )
        for sol in m.fs.solutes:
            print(
                f'    {pyo.value(m.fs.precipitator["retentate"].solid.flow_mass_solute[0, sol])}'
            )
            data["precipitate products retentate"].append(
                pyo.value(m.fs.precipitator["retentate"].solid.flow_mass_solute[0, sol])
            )
        print("^retentate vpermeate")
        data["precipitate products permeate"] = []
        print(pyo.value(m.fs.precipitator["permeate"].solid.flow_vol[0]))
        data["precipitate products permeate"].append(
            pyo.value(m.fs.precipitator["permeate"].solid.flow_vol[0])
        )
        for sol in m.fs.solutes:
            print(
                f'    {pyo.value(m.fs.precipitator["permeate"].solid.flow_mass_solute[0, sol])}'
            )
            data["precipitate products permeate"].append(
                pyo.value(m.fs.precipitator["permeate"].solid.flow_mass_solute[0, sol])
            )
        print("\nprecipitate recycle")
        data["precipitate recycle retentate"] = []
        print(pyo.value(m.fs.precipitator["retentate"].recycle.flow_vol[0]))
        data["precipitate recycle retentate"].append(
            pyo.value(m.fs.precipitator["retentate"].recycle.flow_vol[0])
        )
        for sol in m.fs.solutes:
            print(
                f'    {pyo.value(m.fs.precipitator["retentate"].recycle.flow_mass_solute[0, sol])}'
            )
            data["precipitate recycle retentate"].append(
                pyo.value(
                    m.fs.precipitator["retentate"].recycle.flow_mass_solute[0, sol]
                )
            )
        print("^retentate vpermeate")
        data["precipitate recycle permeate"] = []
        print(pyo.value(m.fs.precipitator["permeate"].recycle.flow_vol[0]))
        data["precipitate recycle permeate"].append(
            pyo.value(m.fs.precipitator["permeate"].recycle.flow_vol[0])
        )
        for sol in m.fs.solutes:
            print(
                f'    {pyo.value(m.fs.precipitator["permeate"].recycle.flow_mass_solute[0, sol])}'
            )
            data["precipitate recycle permeate"].append(
                pyo.value(
                    m.fs.precipitator["permeate"].recycle.flow_mass_solute[0, sol]
                )
            )

        print("\nwaste")
        data["waste"] = []
        print(pyo.value(m.fs.split_precipitate_recycle.waste.flow_vol[0]))
        data["waste"].append(
            pyo.value(m.fs.split_precipitate_recycle.waste.flow_vol[0])
        )
        for sol in m.fs.solutes:
            print(
                f"    {pyo.value(m.fs.split_precipitate_recycle.waste.flow_mass_solute[0, sol])}"
            )
            data["waste"].append(
                pyo.value(m.fs.split_precipitate_recycle.waste.flow_mass_solute[0, sol])
            )

        print("\nsplit fractions for precipitator")
        data["prec sf retentate"] = []
        for sol in m.fs.solutes:
            print(pyo.value(m.fs.precipitator["retentate"].yields[sol, "solid"]))
            data["prec sf retentate"].append(
                pyo.value(m.fs.precipitator["retentate"].yields[sol, "solid"])
            )
        print("^retentate vpermeate")
        data["prec sf permeate"] = []
        for sol in m.fs.solutes:
            print(pyo.value(m.fs.precipitator["permeate"].yields[sol, "solid"]))
            data["prec sf permeate"].append(
                pyo.value(m.fs.precipitator["permeate"].yields[sol, "solid"])
            )

        print("\nactual Co/Li Recovery")
        print(pyo.value(m.prec_perc_co))
        print(pyo.value(m.prec_perc_li))
        data["actual recovery"] = [pyo.value(m.prec_perc_co), pyo.value(m.prec_perc_li)]
        m.fs.precipitator["retentate"].volume.pprint()
        m.fs.precipitator["permeate"].volume.pprint()
        data["prec volumes"] = [
            pyo.value(m.fs.precipitator["retentate"].volume),
            pyo.value(m.fs.precipitator["permeate"].volume),
        ]
    return data


def visualize_flows(num_boxes, num_sub_boxes, conf="stage", model=None):
    """Plot flowsheet of diafiltration system."""
    # Parameters for the boxes
    box_size = 1  # Size of each box (width and height)
    spacing = 0.3  # Space between the boxes
    thick_line_width = 2  # Thickness of the outer box lines
    tol = 0.001
    legflagret = True
    legflagperm = True
    legflagrec = True
    legflagfeed = True
    legflagdiaf = True
    if model is None:
        model = {
            "size": [1] * num_boxes * num_sub_boxes,
            "solutes": [1] * num_boxes * num_sub_boxes,
            "Membrane outlet recoveries": [1] * num_boxes * num_sub_boxes,
            "Membrane Area": [1] * num_boxes * num_sub_boxes,
            "feed flows": [1] * num_boxes * num_sub_boxes,
            "diaf flows": [1] * num_boxes * num_sub_boxes,
            "recycle flows": [1] * num_boxes * num_sub_boxes,
            "inlet splitter flows": [1] * num_boxes * num_sub_boxes,
            "Co product flows": [1] * num_boxes * 3,
            "Li product flows": [1] * num_boxes * 3,
            "diafiltrate inlet": [1] * num_boxes * num_sub_boxes,
            "precipitate products retentate": [1] * num_boxes * num_sub_boxes,
            "precipitate products permeate": [1] * num_boxes * num_sub_boxes,
            "precipitate recycle retentate": [1] * num_boxes * num_sub_boxes,
            "precipitate recycle permeate": [1] * num_boxes * num_sub_boxes,
            "waste": [1] * num_boxes * num_sub_boxes,
            "prec sf retentate": [1] * num_boxes * num_sub_boxes,
            "prec sf permeate": [1] * num_boxes * num_sub_boxes,
            "actual recovery": [1] * num_boxes * num_sub_boxes,
            "prec volumes": [1] * num_boxes * num_sub_boxes,
        }

    # Create a new figure and axis
    fig, ax = plt.subplots(
        figsize=(2 * (num_boxes * (box_size + spacing)), 2 * (box_size + 1))
    )

    # Draw each box with sub-boxes
    arrow_positions = []  # Collect arrow positions to draw later
    for i in range(num_boxes):
        # Calculate the position of the current box
        x_start = i * (box_size + spacing)
        y_start = 0

        # Draw the outer box
        outer_box = patches.Rectangle(
            (x_start, y_start),
            box_size,
            box_size,
            linewidth=thick_line_width,
            edgecolor="black",
            facecolor="none",
        )
        ax.add_patch(outer_box)

        if i != num_boxes - 1:
            # Draw arrows from bottom of box to top of next box
            ax.arrow(
                x_start + box_size + 0.03,  # Start X position (shifted further right)
                y_start + 0.05,  # Y position (above bottom border)
                (spacing - 0.12) / 2,  # X direction (shorter arrow)
                0,  # Y direction
                head_width=0.02,
                head_length=0.03,
                fc=to_rgba("black"),
                ec=to_rgba("black"),
            )
            ax.arrow(
                x_start
                + box_size
                + (spacing) / 2,  # Start X position (shifted further right)
                y_start + 0.05,  # Y position (above bottom border)
                0,  # X direction (shorter arrow)
                box_size - 0.13,  # Y direction
                head_width=0.02,
                head_length=0.03,
                fc=to_rgba("black"),
                ec=to_rgba("black"),
            )
            ax.arrow(
                x_start
                + box_size
                + (spacing) / 2,  # Start X position (shifted further right)
                y_start + box_size - 0.05,  # Y position (above bottom border)
                (spacing - 0.07) / 2,  # X direction (shorter arrow)
                0,  # Y direction
                head_width=0.02,
                head_length=0.03,
                fc=to_rgba("black"),
                ec=to_rgba("black"),
            )

        # Draw arrows for product streams
        if i == num_boxes - 1 or conf != "classic":
            if model["Li product flows"][i * 3] > tol:
                ax.arrow(
                    x_start
                    + box_size
                    + 0.03,  # Start X position (shifted further right)
                    y_start + 0.05,  # Y position (above bottom border)
                    0,  # X direction (shorter arrow)
                    -0.25,  # Y direction
                    head_width=0.02,
                    head_length=0.03,
                    fc=to_rgba("darkorange"),
                    ec=to_rgba("darkorange"),
                    label="Permeate Product" if legflagperm else "",
                )
                legflagperm = False
        if i == 0:
            prodlength = 0.3 + (num_boxes + 1) * 0.1
            prodcol = to_rgba("blue")
        else:
            prodlength = 0.3 + i * 0.1
            prodcol = to_rgba("black")
        ax.arrow(
            x_start + box_size + 0.03,  # Start X position (shifted further right)
            y_start + box_size - 0.05,  # Y position (above bottom border)
            0,  # X direction (shorter arrow)
            prodlength,  # Y direction
            head_width=0.02,
            head_length=0.03,
            fc=prodcol,
            ec=prodcol,
            label="Retentate Product" if legflagret else "",
        )
        legflagret = False
        if conf != "classic":
            if i != 0 and model["Co product flows"][i * 3] > tol:
                ax.arrow(
                    x_start
                    + box_size
                    + 0.03,  # Start X position (shifted further right)
                    y_start
                    + box_size
                    - 0.02
                    + prodlength,  # Y position (above bottom border)
                    0,  # X direction (shorter arrow)
                    0.27 + (num_boxes + 1) * 0.1 - prodlength,  # Y direction
                    head_width=0.02,
                    head_length=0.03,
                    fc=to_rgba("blue"),
                    ec=to_rgba("blue"),
                    label="Retentate Product" if legflagret else "",
                )
                legflagret = False

        # Draw sub-boxes and middle rectangles inside the current box
        sub_box_width = box_size / num_sub_boxes
        for j in range(num_sub_boxes):
            # Draw sub-box
            sub_box = patches.Rectangle(
                (x_start + j * sub_box_width, y_start),
                sub_box_width,
                box_size,
                linewidth=1,
                edgecolor="darkblue",
                facecolor=to_rgba("lightblue"),
            )
            ax.add_patch(sub_box)

            # Draw the middle rectangle inside each sub-box
            middle_rect = patches.Rectangle(
                (x_start + j * sub_box_width, y_start + box_size / 2 - 0.0333),
                sub_box_width,
                0.0667,
                linewidth=1,
                edgecolor="darkblue",
                facecolor=to_rgba("bisque"),
            )
            ax.add_patch(middle_rect)

            # Collect positions for arrows
            arrow_positions.append((x_start, y_start, sub_box_width, box_size, j, i))

    # Draw arrows after all boxes are plotted
    for x_start, y_start, sub_box_width, box_size, j, i in arrow_positions:
        if conf == "tube":
            # Add a vertical arrow at the top of the box
            ax.arrow(
                x_start + j * sub_box_width + sub_box_width / 2,  # X position
                y_start
                + box_size
                + 0.15,  # Start Y position (slightly above middle rect)
                0,  # X direction
                -0.12,  # Y direction (longer arrow)
                head_width=0.02,
                head_length=0.03,
                fc=to_rgba("black"),
                ec=to_rgba("black"),
            )
            # add inlet mixers
            xloc = x_start + j * sub_box_width + sub_box_width / 2
            yloc = y_start + box_size + 0.15
            ax.add_patch(
                patches.Polygon(
                    [
                        [xloc, yloc],
                        [xloc + sub_box_width / 2 - sub_box_width / 8, yloc + 0.05],
                        [xloc - sub_box_width / 2 + sub_box_width / 8, yloc + 0.05],
                    ],
                    edgecolor="black",
                    facecolor="none",
                )
            )

            # Add recycle streams
            if i != num_boxes - 1:
                if model["recycle flows"][(num_sub_boxes) * i + j] > tol:
                    ax.arrow(
                        xloc - sub_box_width / 2 + sub_box_width / 8,  # X position
                        yloc
                        + 0.23
                        + i * 0.1,  # Start Y position (slightly above middle rect)
                        0,  # X direction
                        -0.14 - 0.1 * i,  # Y direction (longer arrow)
                        head_width=0.02,
                        head_length=0.03,
                        fc=to_rgba("darkcyan"),
                        ec=to_rgba("darkcyan"),
                        label="Recycle" if legflagrec else "",
                    )
                    legflagrec = False
                if model["recycle flows"][(num_sub_boxes) * i + j] > tol:
                    ax.plot(
                        [
                            xloc - sub_box_width / 2 + sub_box_width / 8,
                            (i + 1) * (box_size + spacing) + box_size + 0.03,
                        ],
                        [yloc + 0.23 + i * 0.1, yloc + 0.23 + i * 0.1],
                        c=to_rgba("darkcyan"),
                    )

            # Add feed
            if model["feed flows"][(num_sub_boxes) * i + j] > tol:
                ax.arrow(
                    x_start + j * sub_box_width + sub_box_width / 2,  # X position
                    yloc
                    + (num_boxes + 1)
                    * 0.1,  # Start Y position (slightly above middle rect)
                    0,  # X direction
                    0.09 - (num_boxes + 1) * 0.1,  # Y direction (longer arrow)
                    head_width=0.02,
                    head_length=0.03,
                    fc=to_rgba("red"),
                    ec=to_rgba("red"),
                    label="Feed" if legflagfeed else "",
                )
                legflagfeed = False
            # connect tops of the arrows
            if j == 0:
                if i == 0 and model["feed flows"][(num_sub_boxes) * i + j] > tol:
                    ax.plot(
                        [xloc, 0],
                        [yloc + (num_boxes + 1) * 0.1, yloc + (num_boxes + 1) * 0.1],
                        c=to_rgba("red"),
                    )
                elif model["feed flows"][(num_sub_boxes) * i + j] > tol:
                    ax.plot(
                        [xloc, 0],
                        [yloc + (num_boxes + 1) * 0.1, yloc + (num_boxes + 1) * 0.1],
                        c=to_rgba("red"),
                    )
            else:
                if model["feed flows"][(num_sub_boxes) * i + j] > tol:
                    ax.plot(
                        [xloc, 0],
                        [yloc + (num_boxes + 1) * 0.1, yloc + (num_boxes + 1) * 0.1],
                        c=to_rgba("red"),
                    )

            # Add diafiltrate
            if model["diaf flows"][(num_sub_boxes) * i + j] > tol:
                ax.arrow(
                    xloc + sub_box_width / 2 - sub_box_width / 8,  # X position
                    yloc
                    + 0.06
                    + (num_boxes + 1)
                    * 0.1,  # Start Y position (slightly above middle rect)
                    0,  # X direction
                    0.03 - (num_boxes + 1) * 0.1,  # Y direction (longer arrow)
                    head_width=0.02,
                    head_length=0.03,
                    fc=to_rgba("darkviolet"),
                    ec=to_rgba("darkviolet"),
                    label="Diafiltrate" if legflagdiaf else "",
                )
                legflagdiaf = False
            # connect tops of the arrows
            if j == 0:
                if i == 0 and model["diaf flows"][(num_sub_boxes) * i + j] > tol:
                    ax.plot(
                        [xloc + sub_box_width / 2 - sub_box_width / 8, 0],
                        [
                            yloc + 0.06 + (num_boxes + 1) * 0.1,
                            yloc + 0.06 + (num_boxes + 1) * 0.1,
                        ],
                        c=to_rgba("darkviolet"),
                    )
                elif model["diaf flows"][(num_sub_boxes) * i + j] > tol:
                    ax.plot(
                        [xloc + sub_box_width / 2 - sub_box_width / 8, 0],
                        [
                            yloc + 0.06 + (num_boxes + 1) * 0.1,
                            yloc + 0.06 + (num_boxes + 1) * 0.1,
                        ],
                        c=to_rgba("darkviolet"),
                    )
            else:
                if model["diaf flows"][(num_sub_boxes) * i + j] > tol:
                    ax.plot(
                        [xloc + sub_box_width / 2 - sub_box_width / 8, 0],
                        [
                            yloc + 0.06 + (num_boxes + 1) * 0.1,
                            yloc + 0.06 + (num_boxes + 1) * 0.1,
                        ],
                        c=to_rgba("darkviolet"),
                    )

        if conf == "stage":
            # Add a vertical arrow at the top of the box
            if model["inlet splitter flows"][(num_sub_boxes) * i + j] > tol:
                ax.arrow(
                    x_start + j * sub_box_width + sub_box_width / 2,  # X position
                    y_start
                    + box_size
                    + 0.15,  # Start Y position (slightly above middle rect)
                    0,  # X direction
                    -0.12,  # Y direction (longer arrow)
                    head_width=0.02,
                    head_length=0.03,
                    fc=to_rgba("black"),
                    ec=to_rgba("black"),
                )
            # add inlet mixers
            xloc = x_start + j * sub_box_width + sub_box_width / 2
            yloc = y_start + box_size + 0.15
            if j == 0:
                ax.add_patch(
                    patches.Polygon(
                        [
                            [xloc, yloc],
                            [xloc + sub_box_width / 2 - sub_box_width / 8, yloc + 0.05],
                            [xloc - sub_box_width / 2 + sub_box_width / 8, yloc + 0.05],
                        ],
                        edgecolor="black",
                        facecolor="none",
                    )
                )
            # connect tops of the arrows
            if j != 0 and model["inlet splitter flows"][(num_sub_boxes) * i + j] > tol:
                ax.plot([xloc, x_start + sub_box_width / 2], [yloc, yloc], c="k")

            # Add recycle streams
            if i != num_boxes - 1:
                if j == 0 and model["recycle flows"][i] > tol:
                    ax.arrow(
                        xloc - sub_box_width / 2 + sub_box_width / 8,  # X position
                        yloc
                        + 0.23
                        + i * 0.1,  # Start Y position (slightly above middle rect)
                        0,  # X direction
                        -0.14 - 0.1 * i,  # Y direction (longer arrow)
                        head_width=0.02,
                        head_length=0.03,
                        fc=to_rgba("darkcyan"),
                        ec=to_rgba("darkcyan"),
                        label="Recycle" if legflagrec else "",
                    )
                    legflagrec = False
                if model["recycle flows"][i] > tol:
                    ax.plot(
                        [
                            xloc - sub_box_width / 2 + sub_box_width / 8,
                            (i + 1) * (box_size + spacing) + box_size + 0.03,
                        ],
                        [yloc + 0.23 + i * 0.1, yloc + 0.23 + i * 0.1],
                        c=to_rgba("darkcyan"),
                    )

            # Add feed
            if j == 0 and model["feed flows"][i] > tol:
                ax.arrow(
                    x_start + j * sub_box_width + sub_box_width / 2,  # X position
                    yloc
                    + (num_boxes + 1)
                    * 0.1,  # Start Y position (slightly above middle rect)
                    0,  # X direction
                    0.09 - (num_boxes + 1) * 0.1,  # Y direction (longer arrow)
                    head_width=0.02,
                    head_length=0.03,
                    fc=to_rgba("red"),
                    ec=to_rgba("red"),
                    label="Feed" if legflagfeed else "",
                )
                legflagfeed = False
            # connect tops of the arrows
            if j == 0:
                if i == 0 and model["feed flows"][i] > tol:
                    ax.plot(
                        [xloc, 0],
                        [yloc + (num_boxes + 1) * 0.1, yloc + (num_boxes + 1) * 0.1],
                        c=to_rgba("red"),
                    )
                elif model["feed flows"][i] > tol:
                    ax.plot(
                        [xloc, 0],
                        [yloc + (num_boxes + 1) * 0.1, yloc + (num_boxes + 1) * 0.1],
                        c=to_rgba("red"),
                    )

            # Add diafiltrate
            if j == 0 and model["diaf flows"][i] > tol:
                ax.arrow(
                    xloc + sub_box_width / 2 - sub_box_width / 8,  # X position
                    yloc
                    + 0.06
                    + (num_boxes + 1)
                    * 0.1,  # Start Y position (slightly above middle rect)
                    0,  # X direction
                    0.03 - (num_boxes + 1) * 0.1,  # Y direction (longer arrow)
                    head_width=0.02,
                    head_length=0.03,
                    fc=to_rgba("darkviolet"),
                    ec=to_rgba("darkviolet"),
                    label="Diafiltrate" if legflagdiaf else "",
                )
                legflagdiaf = False
            # connect tops of the arrows
            if j == 0:
                if i == 0 and model["diaf flows"][i] > tol:
                    ax.plot(
                        [xloc + sub_box_width / 2 - sub_box_width / 8, 0],
                        [
                            yloc + 0.06 + (num_boxes + 1) * 0.1,
                            yloc + 0.06 + (num_boxes + 1) * 0.1,
                        ],
                        c=to_rgba("darkviolet"),
                    )
                elif model["diaf flows"][i] > tol:
                    ax.plot(
                        [xloc + sub_box_width / 2 - sub_box_width / 8, 0],
                        [
                            yloc + 0.06 + (num_boxes + 1) * 0.1,
                            yloc + 0.06 + (num_boxes + 1) * 0.1,
                        ],
                        c=to_rgba("darkviolet"),
                    )
        if conf == "classic":
            # Add a vertical arrow at the top of the box
            if j == 0 and model["inlet splitter flows"][(num_sub_boxes) * i + j] > tol:
                ax.arrow(
                    x_start + j * sub_box_width + sub_box_width / 2,  # X position
                    y_start
                    + box_size
                    + 0.15,  # Start Y position (slightly above middle rect)
                    0,  # X direction
                    -0.12,  # Y direction (longer arrow)
                    head_width=0.02,
                    head_length=0.03,
                    fc=to_rgba("black"),
                    ec=to_rgba("black"),
                )
            # add inlet mixers
            xloc = x_start + j * sub_box_width + sub_box_width / 2
            yloc = y_start + box_size + 0.15
            if j == 0:
                ax.add_patch(
                    patches.Polygon(
                        [
                            [xloc, yloc],
                            [xloc + sub_box_width / 2 - sub_box_width / 8, yloc + 0.05],
                            [xloc - sub_box_width / 2 + sub_box_width / 8, yloc + 0.05],
                        ],
                        edgecolor="black",
                        facecolor="none",
                    )
                )
            # Add recycle streams
            if i != num_boxes - 1:
                if j == 0 and model["recycle flows"][i] > tol:
                    ax.arrow(
                        xloc - sub_box_width / 2 + sub_box_width / 8,  # X position
                        yloc
                        + 0.23
                        + i * 0.1,  # Start Y position (slightly above middle rect)
                        0,  # X direction
                        -0.14 - 0.1 * i,  # Y direction (longer arrow)
                        head_width=0.02,
                        head_length=0.03,
                        fc=to_rgba("darkcyan"),
                        ec=to_rgba("darkcyan"),
                        label="Recycle" if legflagrec else "",
                    )
                    legflagrec = False
                if model["recycle flows"][i] > tol:
                    ax.plot(
                        [
                            xloc - sub_box_width / 2 + sub_box_width / 8,
                            (i + 1) * (box_size + spacing) + box_size + 0.03,
                        ],
                        [yloc + 0.23 + i * 0.1, yloc + 0.23 + i * 0.1],
                        c=to_rgba("darkcyan"),
                    )

            # Add feed
            if j == 0 and model["feed flows"][i] > tol:
                ax.arrow(
                    x_start + j * sub_box_width + sub_box_width / 2,  # X position
                    yloc
                    + (num_boxes + 1)
                    * 0.1,  # Start Y position (slightly above middle rect)
                    0,  # X direction
                    0.09 - (num_boxes + 1) * 0.1,  # Y direction (longer arrow)
                    head_width=0.02,
                    head_length=0.03,
                    fc=to_rgba("red"),
                    ec=to_rgba("red"),
                    label="Feed" if legflagfeed else "",
                )
                legflagfeed = False
            # connect tops of the arrows
            if j == 0:
                if i == 0 and model["feed flows"][i] > tol:
                    ax.plot(
                        [xloc, 0],
                        [yloc + (num_boxes + 1) * 0.1, yloc + (num_boxes + 1) * 0.1],
                        c=to_rgba("red"),
                    )
                elif model["feed flows"][i] > tol:
                    ax.plot(
                        [xloc, 0],
                        [yloc + (num_boxes + 1) * 0.1, yloc + (num_boxes + 1) * 0.1],
                        c=to_rgba("red"),
                    )

            # Add diafiltrate
            if j == 0 and model["diaf flows"][i] > tol:
                ax.arrow(
                    xloc + sub_box_width / 2 - sub_box_width / 8,  # X position
                    yloc
                    + 0.06
                    + (num_boxes + 1)
                    * 0.1,  # Start Y position (slightly above middle rect)
                    0,  # X direction
                    0.03 - (num_boxes + 1) * 0.1,  # Y direction (longer arrow)
                    head_width=0.02,
                    head_length=0.03,
                    fc=to_rgba("darkviolet"),
                    ec=to_rgba("darkviolet"),
                    label="Diafiltrate" if legflagdiaf else "",
                )
                legflagdiaf = False
            # connect tops of the arrows
            if j == 0:
                if i == 0 and model["diaf flows"][i] > tol:
                    ax.plot(
                        [xloc + sub_box_width / 2 - sub_box_width / 8, 0],
                        [
                            yloc + 0.06 + (num_boxes + 1) * 0.1,
                            yloc + 0.06 + (num_boxes + 1) * 0.1,
                        ],
                        c=to_rgba("darkviolet"),
                    )
                elif model["diaf flows"][i] > tol:
                    ax.plot(
                        [xloc + sub_box_width / 2 - sub_box_width / 8, 0],
                        [
                            yloc + 0.06 + (num_boxes + 1) * 0.1,
                            yloc + 0.06 + (num_boxes + 1) * 0.1,
                        ],
                        c=to_rgba("darkviolet"),
                    )

        # Add a vertical arrow through the middle rectangle
        ax.arrow(
            x_start + j * sub_box_width + sub_box_width / 2,  # X position
            y_start
            + box_size / 2
            + 0.07,  # Start Y position (slightly above middle rect)
            0,  # X direction
            -0.15,  # Y direction (longer arrow)
            head_width=0.02,
            head_length=0.03,
            fc=to_rgba("sienna"),
            ec=to_rgba("sienna"),
        )

        # Add a horizontal arrow slightly below the top border of each sub-box
        ax.arrow(
            x_start
            + j * sub_box_width
            + 0.5 * sub_box_width,  # Start X position (shifted further right)
            y_start + box_size - 0.05,  # Y position (below top border)
            0.45 * sub_box_width,  # X direction (shorter arrow)
            0,  # Y direction
            head_width=0.02,
            head_length=0.03,
            fc=to_rgba("blue"),
            ec=to_rgba("blue"),
        )

        # Add a horizontal arrow slightly above the bottom border of each sub-box
        ax.arrow(
            x_start
            + j * sub_box_width
            + 0.5 * sub_box_width,  # Start X position (shifted further right)
            y_start + 0.05,  # Y position (above bottom border)
            0.45 * sub_box_width,  # X direction (shorter arrow)
            0,  # Y direction
            head_width=0.02,
            head_length=0.03,
            fc=to_rgba("darkorange"),
            ec=to_rgba("darkorange"),
        )

    # Adjust axis limits and hide the axes
    ax.set_xlim(-spacing, num_boxes * (box_size + spacing))
    ax.set_ylim(-spacing, box_size + spacing + num_boxes * 0.2)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0), ncol=3)

    # Show the plot
    plt.show()
