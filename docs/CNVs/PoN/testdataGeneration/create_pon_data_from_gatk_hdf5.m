function [reduced_panel, reduced_pinv, fnt_control_matrix, log_normals] = create_pon_data_from_gatk_hdf5(pon_filename)
%
% [reduced_panel, reduced_pinv, fnt_control_matrix] = create_pon_data_from_recapseg_hdf5(pon_filename)
%
% pon_filename -- v5 PoN created by GATK CNV
% 
% reduced_panel -- the panel of eigensamples.  pon.getReducedPanelCounts()
% reduced_pinv -- the pseudoinverse of the reduced panel.  pon.getReducedPanelPInverseCounts()
% fnt_control_matrix -- {targets x samples} Each value is normalized by the target factor for the given row.
%       pon.getNormalizedCounts();
%

% Read the fnt_control_matrix
fnt_field_name = '/fnt_control_matrix/block0_values';
fnt_control_matrix = h5read(pon_filename, fnt_field_name)';

% Read the reduced pinv
reduced_pinv_field_name = '/reduced_pon_pinv/block0_values';
reduced_pinv = h5read(pon_filename, reduced_pinv_field_name)';

% Read the reduced panel
reduced_panel_field_name = '/reduced_pon/block0_values';
reduced_panel = h5read(pon_filename, reduced_panel_field_name)';

% Read the reduced panel
log_normals_field_name = '/log_normals/block0_values';
log_normals = h5read(pon_filename, log_normals_field_name)';
