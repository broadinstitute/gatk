% THIS FILE IS FOR RECORD KEEPING ONLY
%  THERE IS NO SUPPORT FOR THIS FILE AND IS MEANT FOR GATK DEV ONLY
%  For usage of this data, see TangentNormalizerUnitTest.java
% LTL January 6, 2016
% Using Matlab 2015a

filenames = [{'create-pon-some-targets.pon'} {'create-pon-all-targets.pon'}];

% tangent-normalize the normals (normalized by target factors) into the reduced panel
for j = 1:length(filenames)
    pon_file = filenames{j};
    [reduced_panel, reduced_pinv, fnt_control_matrix, log_normals] = create_pon_data_from_gatk_hdf5(pon_file);

    % Tangent normalize the fnt_control_matrix
    tangent_norm_results = zeros(size(reduced_panel,1), size(fnt_control_matrix, 2));
    for i = 1:size(fnt_control_matrix, 2)
        sample = fnt_control_matrix(:,i);
        sample = log2(sample);
        sample = sample - median(sample);
        
        beta_hat = reduced_pinv * sample;
        orig_plane_projection = reduced_panel * beta_hat;
        tangent_norm_results(:,i) = sample - orig_plane_projection;
    end

    % If this is being done on windows a dos2unix command is still required
    % for each output
    disp('Output files will still need a dos2unix')
    save([pon_file '.normal_projection'], 'tangent_norm_results', '-ascii', '-tabs', '-double');
end