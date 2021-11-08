import paraview
from paraview.simple import *
from paraview import vtk
import numpy as np
import matplotlib.pyplot as plt

import subprocess


def post_process_3d_data(csv_filename):
    data = np.genfromtxt(csv_filename, delimiter=",", skip_header=1)
    d = np.delete(data, 2, axis=1)



    d = d[np.lexsort((d[:,1],d[:,0]))]
    print(d)
    
    x = d[:, 0]
    y = d[:, 1]

    w = d[:, 2]
    Ix = d[:, 3]
    Iy = d[:, 4]
    Iz = d[:, 5]
    In = np.sqrt(Ix * Ix + Iy * Iy + Iz * Iz)


    xmin = np.min(x)
    xmax = np.max(x)

    ymin = np.min(y)
    ymax = np.max(y)

    print("BoundingBox", xmin, xmax, ymin, ymax )


    dx = 0
    #Search for dx, data are sorted on x
    for val in x:
        if abs(val - x[0]) > 1e-12:
            dx = abs(x[0] - val)
            break

    #Search for dy
    dy = np.abs(y[1] - y[0])

    print("dx dy", dx, dy)
    
    Nx = np.int32(abs(xmax - xmin) / dx)
    Ny = np.int32(abs(ymax - ymin) / dy)
    

    z = np.full((Nx + 1, Ny + 1), np.nan, dtype=np.float64)
    
    idx_geo = 0
    for i in range(0, Nx + 1):
        grid_x = xmin + i * dx
        
        if abs(grid_x - x[idx_geo]) < 1e-6:
            for k in range(Ny + 1):
                grid_y = ymin + k * dy
                if abs(grid_y - y[idx_geo]) < 1e-6:
                    # print(grid_x, grid_y, w[idx_geo])
                    z[i][k] = w[idx_geo]
                    idx_geo += 1
        
                if idx_geo == len(x):
                    break

    print(Nx)
    print(Ny)


    f = open(csv_filename, "w")
    for i in range(Nx + 1):
         for k in range(Ny + 1):
            f.write(f"{xmin + i * dx} {ymin + k * dy} {z[i][k]}\n")
    f.close()

# Points_0,Points_1,Points_2,w0,w1,w2,w3
# 0,0.75,0.5,0.0001,0,0,0
# 0,0.8125,0.5,0.0001,0,0,0
# 0,0.875,0.5,0.0001,0,0,0
# 0,0.9375,0.5,0.0001,0,0,0
# 0,1,0.5,0.0001,0,0,0
# 0,0.1875,0.5,0.0001,0,0,0
# 0,0.25,0.5,0.0001,0,0,0
# 0,0.125,0.5,0.0001,0,0,0
# 0,0.0625,0.5,0.0001,0,0,0
# 0,0,0.5,0.0001,0,0,0
# 0,0.5,0.5,0.0001,0,0,0
# 0,0.5625,0.5,0.0001,0,0,0
# 0,0.625,0.5,0.0001,0,0,0
# 0,0.6875,0.5,0.0001,0,0,0
# 0,0.3125,0.5,0.0001,0,0,0
# 0,0.375,0.5,0.0001,0,0,0
# 0,0.4375,0.5,0.0001,0,0,0



def paraview_extract_3d_slice_to_csv(
    xdmf_filename, csv_filename, time, slice_origin_pt, slice_normal_vec
):
    xdmf3 = Xdmf3ReaderT(FileName=[xdmf_filename])

    xdmf3.CellArrays = ["w0", "w1", "w2", "w3"]
    # Create a new 'Slice'
    slice = Slice(Input=xdmf3)
    slice.SliceType = "Plane"
    slice.SliceOffsetValues = [0.0]

    # Init the 'Plane' selected for 'SliceType'
    slice.SliceType.Origin = slice_origin_pt

    # Properties modified on slice.SliceType
    slice.SliceType.Normal = slice_normal_vec

    # Smooth data using mean over the cell
    data = CellDatatoPointData(Input=slice)
    data.CellDataArraytoprocess = ["w0", "w1", "w2", "w3"]

    animation_scene = GetAnimationScene()
    time_keeper = GetTimeKeeper()

    # Select time step
    time_step_list = animation_scene.TimeKeeper.TimestepValues

    # Find closest time value
    closest_t = min(time_step_list, key=lambda x: abs(x - time))
    if abs(closest_t - time) > 1e-6:
        print("Warning")

    # Setup time step
    animation_scene.AnimationTime = closest_t
    time_keeper.Time = closest_t

    # Create and configure a new 'SpreadSheet View'
    spreadsheet_view = CreateView("SpreadSheetView")
    spreadsheet_view.ColumnToSort = ""
    spreadsheet_view.BlockSize = 1024
    spreadsheet_view.HiddenColumnLabels = ["_Magnitude", "Point ID"]

    # Ascending sort by x value
    spreadsheet_view.HiddenColumnLabels = ["Points_Magnitude", "Point ID"]
    # spreadsheet_view.ColumnToSort = "Points_1"
    spreadsheet_view.ColumnToSort = "Points_0"
    spreadsheet_view.InvertOrder = 0

    # Generate the view
    display = Show(data, spreadsheet_view)

    # Write csv file
    ExportView(csv_filename, view=spreadsheet_view)


def paraview_extract_2d_to_csv(xdmf_filename, csv_filename, time):
    xdmf3 = Xdmf3ReaderT(FileName=[xdmf_filename])

    xdmf3.CellArrays = ["w0", "w1", "w2"]
    # Smooth data using mean over the cell
    data = CellDatatoPointData(Input=xdmf3)
    data.CellDataArraytoprocess = ["w0", "w1", "w2"]

    animation_scene = GetAnimationScene()
    time_keeper = GetTimeKeeper()

    # Select time step
    time_step_list = animation_scene.TimeKeeper.TimestepValues

    # Find closest time value
    closest_t = min(time_step_list, key=lambda x: abs(x - time))
    if abs(closest_t - time) > 1e-6:
        print("Warning")

    # Setup time step
    animation_scene.AnimationTime = closest_t
    time_keeper.Time = closest_t

    # Create and configure a new 'SpreadSheet View'
    spreadsheet_view = CreateView("SpreadSheetView")
    spreadsheet_view.ColumnToSort = ""
    spreadsheet_view.BlockSize = 1024
    spreadsheet_view.HiddenColumnLabels = ["_Magnitude", "Point ID"]

    # Ascending sort by x value : '_0'
    spreadsheet_view.HiddenColumnLabels = ["_Magnitude", "Point ID"]
    spreadsheet_view.ColumnToSort = "_0"

    spreadsheet_view.InvertOrder = 0
    spreadsheet_view.InvertOrder = 1
    spreadsheet_view.InvertOrder = 0
    # Generate the view
    display = Show(data, spreadsheet_view)

    # Write csv file
    ExportView(csv_filename, view=spreadsheet_view)


# xdmf_filename = "/home/gerhard/Desktop/github/reduced-kinetic-models-solver/examples/spherical_harmonics_2d_order1/results.xmf"
# xdmf_filename = "/home/gerhard/Desktop/github/reduced-kinetic-models-solver/examples/spherical_harmonics_3d_order1/results.xmf"
xdmf_filename = "/home/gerhard/Desktop/github/reduced-kinetic-models-solver/examples/ordinates_3d_lebedev_order6/results.xmf"
xdmf_filename = "/home/gerhard/Desktop/github/reduced-kinetic-models-solver/examples/moments_3d_s2_order1/results.xmf"
csv_filename = "foo.csv"
time = 39


# p = subprocess.Popen(['pvbatch', 'test',
#                               self.mod_name, self.func_name],
#                               stdin=subprocess.PIPE,
#                               stdout=subprocess.PIPE)

slice_origin_pt = [5, 1.5, 0.5]
slice_normal_vec = [0, 0, 1]

paraview_extract_3d_slice_to_csv(
    xdmf_filename, csv_filename, time, slice_origin_pt, slice_normal_vec
)
post_process_3d_data(csv_filename)
