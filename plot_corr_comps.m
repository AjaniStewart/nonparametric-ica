load("simulated_voxel_data.mat")

indices = [0,1,2,3,4,5,7,8,10,13,20];
noisy_voxels = struct();
K = 12;
N_RANDOM_INITS = 10;
PLOT_FIGURES = 0;
RAND_SEED = 4;

for i = 1:length(indices)
    varname = sprintf('noisy_voxels_sig_%d', indices(i));
    noisy_voxels.(sprintf('sig_%d', indices(i))) = eval(varname);
end

% Now iterate safely:
fields = fieldnames(noisy_voxels);
for i = 1:length(fields)
    current_voxels = noisy_voxels.(fields{i});
    % voxel decomp
    [R_inferred, W_inferred] = nonparametric_ica(current_voxels, K, N_RANDOM_INITS, PLOT_FIGURES, RAND_SEED);
    % align w/ real components
    [perm, ~, ~] = align_components(real_components, R_inferred);
    R_inferred = R_inferred(:,perm);
    % plot correlation
    figure;
    imagesc(corr(R_inferred, real_components));
    ylabel('Inferred Response Profile');
    xlabel('True Response Profile');
    colorbar;
    axis square;
    title(sprintf("Correlation between inferred and ground truth components, noise level: %d", indices(i)));



end