
CNV PoN fields (v5)
-------------------

Groups:

- ``fnt_control_matrix`` -- normal samples proportional coverage divided by the target factors.  These have been filtered by the targetFactorPercentileThreshold, so there will be fewer targets than what was specified in the initial call to CreatePanelOfNormals.  Size:  \[targets x samples\]  
- ``target_factors`` --  the median coverage per target.  These have been filtered by the targetFactorPercentileThreshold, so there will be fewer targets than what was specified in the initial call to CreatePanelOfNormals.  Size:  \[targets x 1\]
- ``log_normals`` -- normals that have been log'd (base 2) and centered by the median (per column).  These have been filtered through several steps and will have fewer targets and samples than in ``target_factors`` and ``fnt_control_matrix``.  Size: \[targets x samples\] (filtered)
- ``log_normals_pinv`` -- pseudoinverse of ``log_normals``.  Size is transpose of ``log_normals``.
- ``log_normals_targets`` -- the targets used in ``log_normals`` as a string array, in row major form.  See ``num_target_cols``
- ``num_target_cols`` -- targets are encoded as an array, even though a table makes more sense.  This is the number of columns in any target field.  Stored as a single number (3.0) Size: \[1x1\]
- ``raw_targets`` -- Unfiltered target list that was used as input to the CreatePanelOfNormals call.  These are stored as string array, in row major form.  See ``num_target_cols``.  Size \[initial_targets x 3 (num_target_cols)\]
- ``version`` -- PoN version.  Stored as a single double (5.0) Size:  \[1 x 1\]
- ``reduced_pon`` -- right eigenvectors of the SVD of ``log_normals`` with high variance components (columns) preserved.  The number of preserved columns (eigensamples) depends on the mean of the eigenvalues.  Size:  \[targets x eigensamples\]  The corresponding targets are in ``log_normals_targets``.
- ``reduced_pon_inverse`` -- pseudoinverse of the ``reduced_pon``.  Size is transpose of ``reduced_pon``.
- ``target_variances`` -- variance for each target after tangent normalizing each normal (column in ``fnt_control_matrix``) into the ``reduced_pon``.  Same size as ``target_factors``.

Data:

- ``index`` -- names for each row in the ``block0_values`` or ``values``
- ``block0_values`` or ``values`` -- the actual data

