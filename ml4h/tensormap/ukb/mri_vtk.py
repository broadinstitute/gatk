# VTK-specific tensormaps
import os
import h5py
import numpy as np
import vtk
import vtk.util.numpy_support
from ml4h.TensorMap import TensorMap, Interpretation
from ml4h.defines import MRI_SEGMENTED, MRI_LAX_SEGMENTED, MRI_FRAMES, MRI_PIXEL_WIDTH, MRI_PIXEL_HEIGHT, MRI_SLICE_THICKNESS, MRI_PATIENT_ORIENTATION, MRI_PATIENT_POSITION
from ml4h.tensormap.ukb.mri import mri_tensor_2d


def _mri_tensor_4d(hd5, name):
    """
    Returns MRI image tensors from HD5 as 4-D numpy arrays. Useful for raw SAX and LAX images and segmentations.
    """
    if isinstance(hd5[name], h5py.Group):
        nslices = len(hd5[name]) // MRI_FRAMES
        for img in hd5[name]:
            img_shape = hd5[name][img].shape
            break
        shape = (img_shape[0], img_shape[1], nslices, MRI_FRAMES)
        arr = np.zeros(shape)
        t = 0
        s = 0
        for k in sorted(hd5[name], key=int):
            arr[:, :, s, t] = np.array(hd5[name][k]).T
            t += 1
            if t == MRI_FRAMES:
                s += 1
                t = 0
    elif isinstance(hd5[name], h5py.Dataset):
        nslices = 1
        shape = (hd5[name].shape[0], hd5[name].shape[1], nslices, MRI_FRAMES)
        arr = np.zeros(shape)
        for t in range(MRI_FRAMES):
            arr[:, :, 0, t] = np.array(hd5[name][:, :, t]).T
    else:
        raise ValueError(f'{name} is neither a HD5 Group nor a HD5 dataset')
    return arr


def _mri_hd5_to_structured_grids(hd5, name, save_path=None, order='F'):
    """
    Returns MRI tensors as list of VTK structured grids aligned to the reference system of the patient
    """
    arr = _mri_tensor_4d(hd5, name)
    width = hd5['_'.join([MRI_PIXEL_WIDTH, name])]
    height = hd5['_'.join([MRI_PIXEL_HEIGHT, name])]
    positions = mri_tensor_2d(hd5, '_'.join([MRI_PATIENT_POSITION, name]))
    orientations = mri_tensor_2d(
        hd5,
        '_'.join([MRI_PATIENT_ORIENTATION, name]),
    )
    thickness = hd5['_'.join([MRI_SLICE_THICKNESS, name])]
    _, dataset_indices, dataset_counts = np.unique(
        orientations,
        axis=1,
        return_index=True,
        return_counts=True,
    )
    grids = []
    for d_idx, d_cnt in zip(dataset_indices, dataset_counts):
        grids.append(vtk.vtkStructuredGrid())
        nslices = d_cnt
        # If multislice, override thickness as distance between voxel centers. Note: removes eventual gaps between slices
        if nslices > 1:
            thickness = np.linalg.norm(
                positions[:, d_idx] -
                positions[:, d_idx + 1],
            )
        transform = vtk.vtkTransform()
        n_orientation = np.cross(
            orientations[3:, d_idx], orientations[
                :3,
                d_idx
            ],
        )
        # 4x4 transform matrix to align to the patient reference system
        transform.SetMatrix([
            orientations[3, d_idx] * height,
            orientations[0, d_idx] * width,
            n_orientation[0] * thickness,
            positions[0, d_idx],
            orientations[4, d_idx] * height,
            orientations[1, d_idx] * width,
            n_orientation[1] * thickness,
            positions[1, d_idx],
            orientations[5, d_idx] * height,
            orientations[2, d_idx] * width,
            n_orientation[2] * thickness,
            positions[2, d_idx],
            0,
            0,
            0,
            1,
        ])
        x_coors = np.arange(0, arr.shape[0] + 1) - 0.5
        y_coors = np.arange(0, arr.shape[1] + 1) - 0.5
        z_coors = np.arange(0, d_cnt + 1) - 0.5
        xyz_meshgrid = np.meshgrid(x_coors, y_coors, z_coors)
        xyz_pts = np.zeros(
            ((arr.shape[0] + 1) * (arr.shape[1] + 1) * (d_cnt + 1), 3),
        )
        for dim in range(3):
            xyz_pts[:, dim] = xyz_meshgrid[dim].ravel(order=order)
        vtk_pts = vtk.vtkPoints()
        vtk_pts.SetData(vtk.util.numpy_support.numpy_to_vtk(xyz_pts))
        grids[-1].SetPoints(vtk_pts)
        grids[-1].SetDimensions(len(x_coors), len(y_coors), len(z_coors))
        grids[-1].SetExtent(
            0,
            len(x_coors) - 1, 0,
            len(y_coors) - 1, 0,
            len(z_coors) - 1,
        )
        for t in range(MRI_FRAMES):
            arr_vtk = vtk.util.numpy_support.numpy_to_vtk(
                arr[:, :, d_idx:d_idx + d_cnt, t].ravel(order=order),
                deep=True,
            )
            arr_vtk.SetName(f'{name}_{t}')
            grids[-1].GetCellData().AddArray(arr_vtk)
        transform_filter = vtk.vtkTransformFilter()
        transform_filter.SetInputData(grids[-1])
        transform_filter.SetTransform(transform)
        transform_filter.Update()
        grids[-1].DeepCopy(transform_filter.GetOutput())
        if save_path:
            writer = vtk.vtkXMLStructuredGridWriter()
            writer.SetFileName(
                os.path.join(save_path, f'grid_{name}_{d_idx}.vts'),
            )
            writer.SetInputData(grids[-1])
            writer.Update()
    return grids


def _cut_through_plane(dataset, plane_center, plane_orientation):
    plane = vtk.vtkPlane()
    plane.SetOrigin(plane_center)
    plane.SetNormal(plane_orientation)
    cutter = vtk.vtkCutter()
    cutter.SetInputData(dataset)
    cutter.SetCutFunction(plane)
    poly = vtk.vtkDataSetSurfaceFilter()
    poly.SetInputConnection(cutter.GetOutputPort())
    poly.Update()
    return poly.GetOutput()


def _map_points_to_cells(pts, dataset, tol=1e-3):
    locator = vtk.vtkCellLocator()
    locator.SetDataSet(dataset)
    locator.BuildLocator()
    closest_pt = np.zeros(3)
    generic_cell = vtk.vtkGenericCell()
    cell_id, sub_id, dist2, inside = vtk.mutable(0), vtk.mutable(
        0,
    ), vtk.mutable(0.0), vtk.mutable(0)
    map_to_cells = np.zeros(len(pts), dtype=np.int64)
    for pt_id, pt in enumerate(pts):
        if locator.FindClosestPointWithinRadius(
            pt, tol, closest_pt,
            generic_cell, cell_id, sub_id,
            dist2, inside,
        ):
            map_to_cells[pt_id] = cell_id.get()
    return map_to_cells


def _make_mri_projected_segmentation_from_file(
    to_segment_name,
    segmented_name,
    save_path=None,
):
    def mri_projected_segmentation(tm, hd5):
        if segmented_name not in [MRI_SEGMENTED, MRI_LAX_SEGMENTED]:
            raise ValueError(
                f'{segmented_name} is recognized neither as SAX nor LAX segmentation',
            )
        cine_segmented_grids = _mri_hd5_to_structured_grids(
            hd5, segmented_name,
        )
        cine_to_segment_grids = _mri_hd5_to_structured_grids(
            hd5, to_segment_name,
        )
        tensor = np.zeros(tm.shape, dtype=np.float32)
        # Loop through segmentations and datasets
        for ds_i, ds_segmented in enumerate(cine_segmented_grids):
            for ds_j, ds_to_segment in enumerate(cine_to_segment_grids):
                dims = ds_to_segment.GetDimensions()
                pts = vtk.util.numpy_support.vtk_to_numpy(
                    ds_to_segment.GetPoints().GetData(),
                )
                npts_per_slice = dims[0] * dims[1]
                ncells_per_slice = (dims[0] - 1) * (dims[1] - 1)
                n_orientation = (pts[npts_per_slice] - pts[0])
                n_orientation /= np.linalg.norm(n_orientation)
                cell_centers = vtk.vtkCellCenters()
                cell_centers.SetInputData(ds_to_segment)
                cell_centers.Update()
                cell_pts = vtk.util.numpy_support.vtk_to_numpy(
                    cell_centers.GetOutput().GetPoints().GetData(),
                )
                # Loop through dataset slices
                for s in range(dims[2] - 1):
                    slice_center = np.mean(
                        pts[
                            s * npts_per_slice:(s + 2) *
                            npts_per_slice
                        ],
                        axis=0,
                    )
                    slice_cell_pts = cell_pts[
                        s * ncells_per_slice:(s + 1) *
                        ncells_per_slice
                    ]
                    slice_segmented = _cut_through_plane(
                        ds_segmented, slice_center, n_orientation,
                    )
                    map_to_segmented = _map_points_to_cells(
                        slice_cell_pts, slice_segmented,
                    )
                    # Loop through time
                    for t in range(MRI_FRAMES):
                        arr_name = f'{segmented_name}_{t}'
                        segmented_arr = vtk.util.numpy_support.vtk_to_numpy(
                            slice_segmented.GetCellData().GetArray(arr_name),
                        )
                        projected_arr = segmented_arr[map_to_segmented]
                        if len(tm.shape) == 3:
                            tensor[:, :, t] = np.maximum(
                                tensor[:, :, t],
                                projected_arr.reshape(
                                    tm.shape[0],
                                    tm.shape[1],
                                ),
                            )
                        elif len(tm.shape) == 4:
                            tensor[:, :, s, t] = np.maximum(
                                tensor[:, :, s, t],
                                projected_arr.reshape(
                                    tm.shape[0],
                                    tm.shape[1],
                                ),
                            )
                    if save_path:
                        writer_segmented = vtk.vtkXMLPolyDataWriter()
                        writer_segmented.SetInputData(slice_segmented)
                        writer_segmented.SetFileName(
                            os.path.join(
                                save_path,
                                f'{tm.name}_segmented_{ds_i}_{ds_j}_{s}.vtp',
                            ),
                        )
                        writer_segmented.Update()
        return tensor

    return mri_projected_segmentation


cine_segmented_lax_2ch_proj_from_sax = TensorMap(
    'cine_segmented_lax_2ch_proj_from_sax',
    Interpretation.CONTINUOUS,
    shape=(256, 256, 50),
    loss='logcosh',
    tensor_from_file=_make_mri_projected_segmentation_from_file(
        'cine_segmented_lax_2ch', MRI_SEGMENTED,
    ),
)
cine_segmented_lax_3ch_proj_from_sax = TensorMap(
    'cine_segmented_lax_3ch_proj_from_sax',
    Interpretation.CONTINUOUS,
    shape=(256, 256, 50),
    loss='logcosh',
    tensor_from_file=_make_mri_projected_segmentation_from_file(
        'cine_segmented_lax_3ch', MRI_SEGMENTED,
    ),
)
cine_segmented_lax_4ch_proj_from_sax = TensorMap(
    'cine_segmented_lax_4ch_proj_from_sax',
    Interpretation.CONTINUOUS,
    shape=(256, 256, 50),
    loss='logcosh',
    tensor_from_file=_make_mri_projected_segmentation_from_file(
        'cine_segmented_lax_4ch', MRI_SEGMENTED,
    ),
)
cine_segmented_lax_2ch_proj_from_lax = TensorMap(
    'cine_segmented_lax_2ch_proj_from_lax',
    Interpretation.CONTINUOUS,
    shape=(256, 256, 50),
    loss='logcosh',
    tensor_from_file=_make_mri_projected_segmentation_from_file(
        'cine_segmented_lax_2ch', MRI_LAX_SEGMENTED,
    ),
)
cine_segmented_lax_3ch_proj_from_lax = TensorMap(
    'cine_segmented_lax_3ch_proj_from_lax',
    Interpretation.CONTINUOUS,
    shape=(256, 256, 50),
    loss='logcosh',
    tensor_from_file=_make_mri_projected_segmentation_from_file(
        'cine_segmented_lax_3ch', MRI_LAX_SEGMENTED,
    ),
)
cine_segmented_lax_4ch_proj_from_lax = TensorMap(
    'cine_segmented_lax_4ch_proj_from_lax',
    Interpretation.CONTINUOUS,
    shape=(256, 256, 50),
    loss='logcosh',
    tensor_from_file=_make_mri_projected_segmentation_from_file(
        'cine_segmented_lax_4ch', MRI_LAX_SEGMENTED,
    ),
)
