# WARNING : Call this script using pvbatch or pvpython not python3 !
import argparse
import logging
import json
import os
import numbers

import numpy as np
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

logging.basicConfig(
    format="[%(asctime)s] - %(levelname)s - %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)


def delete_constant_column(data, column_id_to_check):
    to_delete = []

    for col_id in column_id_to_check:
        if (data[0][col_id] == data[:, col_id]).all():
            to_delete.append(col_id)
    
    d = np.delete(data, to_delete, 1)
    return d


def pvbatch_subprocess_caller(output_folder):
    from subprocess import run, PIPE

    p = run(
        [
            "pvbatch",
            os.path.abspath(__file__),
            "--d",
            os.path.abspath(output_folder),
        ],
        stdout=PIPE,
        encoding="ascii",
    )
    # @TODO : Check return code !
    p.returncode
    print(p.stdout)


def postprocess_csv_line(csv_filename):
    data = np.genfromtxt(csv_filename, delimiter=",", skip_header=1)
    d = delete_constant_column(data, [0, 1, 2])

    l = d.shape[0]

    # Create alias
    x = d[:, 0]
    w = d[:, 1]
    Ix = d[:, 2]
    Iy = d[:, 3]
    Iz = d[:, 4]
    In = np.sqrt(Ix * Ix + Iy * Iy + Iz * Iz)

    f = open(csv_filename, "w")
    for i in range(l):
        f.write(
            "{xi} {wi} {Ixi} {Iyi} {Izi} {Ini}\n".format(
                xi=x[i], wi=w[i], Ixi=Ix[i], Iyi=Iy[i], Izi=Iz[i], Ini=In[i]
            )
        )
    f.close()


def find_clostest_t(time_step_list, time):
    # Find closest time value
    closest_t = min(time_step_list, key=lambda x: abs(x - time))

    if abs(closest_t - time) > 1e-6:
        logger.warning(
            "Time value: {:>6.8f} is not found !\n \t\t\t\t   - Time value: {:>6.8f} is used !".format(
                time,
                closest_t,
            )
        )
    return closest_t


def paraview_dump_line(
    csv_filename, pv_xdmf3, pv_animation_scene, pv_time_keeper, line_dic, time
):
    pv_line = PlotOverLine(Input=pv_xdmf3, Source="High Resolution Line Source")

    pv_line.Source.Point1 = line_dic["point_1"]
    pv_line.Source.Point2 = line_dic["point_2"]

    
    # Setup time step
    pv_animation_scene.AnimationTime = time
    pv_time_keeper.Time = time
    
    # Create and configure a new 'SpreadSheet View'
    spreadsheet_view = CreateView("SpreadSheetView")
    spreadsheet_view.ColumnToSort = ""
    spreadsheet_view.BlockSize = 1024

    spreadsheet_view.HiddenColumnLabels = [
        "Point ID",
        "Points_Magnitude",
        "arc_length",
        "vtkValidPointMask",
    ]
    # Generate the view
    display = Show(pv_line, spreadsheet_view)

    # Write csv file
    ExportView(csv_filename, view=spreadsheet_view)


def postprocess_csv_slice(csv_filename):
    data = np.genfromtxt(csv_filename, delimiter=",", skip_header=1)
    d = delete_constant_column(data, [0, 1, 2])

    # Sort data in lexicographic order x,y
    d = d[np.lexsort((d[:, 1], d[:, 0]))]

    # Create alias
    x = d[:, 0]
    y = d[:, 1]
    w = d[:, 2]
    Ix = d[:, 3]
    Iy = d[:, 4]
    Iz = d[:, 5]

    # Compute norm of the intensity
    In = np.sqrt(Ix * Ix + Iy * Iy + Iz * Iz)

    # Find bounding box
    xmin = np.min(x)
    xmax = np.max(x)
    ymin = np.min(y)
    ymax = np.max(y)

    # logger.info(
    #     "{:<30}:\n\t- xmin: {:>6.8f}\n\t- xmax: {:>6.8f}\n\t- ymin: {:>6.8f}\n\t- ymax: {:>6.8f}".format(
    #         "Found bounding box", xmin, xmax, ymin, ymax
    #     )
    # )

    # Find dx value (data are sorted on x)
    dx = 0
    for val in x:
        if abs(val - x[0]) > 1e-12:
            dx = abs(x[0] - val)
            break

    # Find dy value
    dy = np.abs(y[1] - y[0])

    # logger.info("{:<30}:\n\t- dx: {:>6.8f}\n\t- dy: {:>6.8f}".format("Found spacial steps", dx, dy))

    # Transfert node data extracted with paraview (on the computationnal domain)
    # on a equally spaced grid that is generated on the bounding box.
    # Trick : When a point of the bounding box grid is located outside the
    # computationnal domain a NaN value is written.
    # This value will be skipped later using pgfplot and displayed in white.

    Nx = np.int32(abs(xmax - xmin) / dx)
    Ny = np.int32(abs(ymax - ymin) / dy)
    z = np.full((Nx + 1, Ny + 1, 5), np.nan, dtype=np.float64)

    idx_geo = 0
    for i in range(0, Nx + 1):
        grid_x = xmin + i * dx
        if abs(grid_x - x[idx_geo]) < 1e-6:
            for k in range(Ny + 1):
                grid_y = ymin + k * dy
                if abs(grid_y - y[idx_geo]) < 1e-6:
                    z[i][k][:] = [
                        w[idx_geo],
                        Ix[idx_geo],
                        Iy[idx_geo],
                        Iz[idx_geo],
                        In[idx_geo],
                    ]
                    idx_geo += 1

                if idx_geo == len(x):
                    break

    f = open(csv_filename, "w")
    for i in range(Nx + 1):
        for k in range(Ny + 1):
            f.write(
                "{x} {y} {w} {Ix} {Iy} {Iz} {In}\n".format(
                    x=xmin + i * dx,
                    y=ymin + k * dy,
                    w=z[i][k][0],
                    Ix=z[i][k][1],
                    Iy=z[i][k][2],
                    Iz=z[i][k][3],
                    In=z[i][k][4],
                )
            )

    f.close()


def paraview_dump_slice(
    csv_filename, pv_xdmf3, pv_animation_scene, pv_time_keeper, slice_dic, time
):
    # Create a new 'Slice'
    pv_slice = Slice(Input=pv_xdmf3)
    pv_slice.SliceType = "Plane"
    pv_slice.SliceOffsetValues = [0.0]

    # Init the 'Plane' selected for 'SliceType'
    pv_slice.SliceType.Origin = slice_dic["point"]

    # Properties modified on pv_slice.SliceType
    pv_slice.SliceType.Normal = slice_dic["normal"]

    # Smooth data using CellDatatoPointData filter
    pv_data = CellDatatoPointData(Input=pv_slice)
    pv_data.CellDataArraytoprocess = ["w0", "w1", "w2", "w3"]

    # Setup time step
    pv_animation_scene.AnimationTime = time
    pv_time_keeper.Time = time

    # Create and configure a new 'SpreadSheet View'
    spreadsheet_view = CreateView("SpreadSheetView")
    spreadsheet_view.ColumnToSort = ""
    spreadsheet_view.BlockSize = 1024
    spreadsheet_view.HiddenColumnLabels = ["_Magnitude", "Point ID"]
    spreadsheet_view.HiddenColumnLabels = ["Points_Magnitude", "Point ID"]

    # Generate the view
    display = Show(pv_data, spreadsheet_view)

    # Write csv file
    ExportView(csv_filename, view=spreadsheet_view)


def paraview_extract_data_to_csv(output_folder):
    xdmf_filename = os.path.join(output_folder, "results.xmf")
    json_filename = os.path.join(output_folder, "postprocess_parameters.json")

    with open(json_filename, "r") as json_file:
        postprocess_parameters = json.load(json_file)

    pv_xdmf3 = Xdmf3ReaderT(FileName=[xdmf_filename])

    if postprocess_parameters["dim"] == 3:
        pv_xdmf3.CellArrays = ["w0", "w1", "w2", "w3"]
    else:
        pv_xdmf3.CellArrays = ["w0", "w1", "w2"]

    pv_animation_scene = GetAnimationScene()
    pv_time_keeper = GetTimeKeeper()

    # Select time step
    time_step_list = pv_animation_scene.TimeKeeper.TimestepValues

    if "time" in postprocess_parameters:
        assert all(
            isinstance(t, numbers.Number) for t in postprocess_parameters["time"]
        )

        for t in postprocess_parameters["time"]:
            close_t = find_clostest_t(time_step_list, t)

            if "line" in postprocess_parameters:
                
                # Smooth data using CellDatatoPointData filter
                pv_data = CellDatatoPointData(Input=pv_xdmf3)
                
                if postprocess_parameters["dim"] == 3:
                    pv_data.CellDataArraytoprocess = ["w0", "w1", "w2", "w3"]
                else :
                    pv_xdmf3.CellArrays = ["w0", "w1", "w2"]
      
                
                for line in postprocess_parameters["line"]:
                    base_folder = os.path.basename(output_folder)

                    csv_filename = "{}_{}_pt1{}_pt2{}_t{}.csv".format(
                        base_folder,
                        "line",
                        "_".join(map(str, line["point_1"])).replace(".", ""),
                        "_".join(map(str, line["point_2"])).replace(".", ""),
                        str(t).replace(".", ""),
                    )

                    csv_filename = os.path.join(output_folder, csv_filename)
                    paraview_dump_line(
                        csv_filename,
                        pv_data,
                        pv_animation_scene,
                        pv_time_keeper,
                        line,
                        t,
                    )

                    postprocess_csv_line(csv_filename)

                    logger.info(
                        "{:<30}: {:<30}".format(
                            "Line saved in", os.path.basename(csv_filename)
                        )
                    )

            if "slice" in postprocess_parameters and postprocess_parameters["dim"] == 3:
                for slice in postprocess_parameters["slice"]:
                    # Check vector value conformity
                    base_folder = os.path.basename(output_folder)

                    csv_filename = "{}_{}_pt{}_vn{}_t{}.csv".format(
                        base_folder,
                        "slice",
                        "_".join(map(str, slice["point"])).replace(".", ""),
                        "_".join(map(str, slice["normal"])).replace(".", ""),
                        str(t).replace(".", ""),
                    )

                    csv_filename = os.path.join(output_folder, csv_filename)
                    paraview_dump_slice(
                        csv_filename,
                        pv_xdmf3,
                        pv_animation_scene,
                        pv_time_keeper,
                        slice,
                        t,
                    )
                    postprocess_csv_slice(csv_filename)

                    logger.info(
                        "{:<30}: {:<30}".format(
                            "Slice saved in", os.path.basename(csv_filename)
                        )
                    )


if __name__ == "__main__":
    from paraview.simple import *

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--d", help="Path to the directory containing xdmf files", type=str
    )
    args = parser.parse_args()

    paraview_extract_data_to_csv(args.d)
