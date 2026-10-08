function export_temperature_field(modelFile, outputDir, x_nm, y_nm, time_s, z_nm)
% Read a solved model only; coordinate units follow its nm geometry.
validateattributes(x_nm, {'numeric'}, {'vector','finite'});
validateattributes(y_nm, {'numeric'}, {'vector','finite'});
assert(all(diff(x_nm)>0) && all(diff(y_nm)>0), 'Axes must increase strictly.');
validateattributes(time_s, {'numeric'}, {'scalar','positive'});
validateattributes(z_nm, {'numeric'}, {'scalar','finite'});
assert(isfile(modelFile), 'Solved model not found.');
assert(~isfolder(outputDir), 'Choose a new output directory.');
addpath(comsol_mli_path());
mphstart(comsol_server_port());
tag = ['FieldExport_' char(datetime('now','Format','HHmmssSSS'))];
model = mphload(modelFile, tag);
cleanup = onCleanup(@() com.comsol.model.util.ModelUtil.remove(tag));
[X, Y] = meshgrid(x_nm, y_nm);
coord = [X(:)'; Y(:)'; z_nm * ones(1,numel(X))];
T = mphinterp(model, 'T', 'coord', coord, 't', time_s);
assert(all(isfinite(T(:))), 'Export includes points outside the model or nonfinite values.');
T = reshape(T, size(X));
xCenter = (x_nm(1)+x_nm(end))/2;
yCenter = (y_nm(1)+y_nm(end))/2;
Tx = mphinterp(model, 'T', 'coord', [x_nm(:)'; yCenter*ones(1,numel(x_nm)); z_nm*ones(1,numel(x_nm))], 't', time_s);
Ty = mphinterp(model, 'T', 'coord', [xCenter*ones(1,numel(y_nm)); y_nm(:)'; z_nm*ones(1,numel(y_nm))], 't', time_s);
mkdir(outputDir);
writematrix(T, fullfile(outputDir,'temperature_map_K.csv'));
writematrix(x_nm(:), fullfile(outputDir,'x_axis_nm.csv'));
writematrix(y_nm(:), fullfile(outputDir,'y_axis_nm.csv'));
writetable(table(x_nm(:),Tx(:),'VariableNames',{'position_nm','temperature_K'}), fullfile(outputDir,'cutline_along_x.csv'));
writetable(table(y_nm(:),Ty(:),'VariableNames',{'position_nm','temperature_K'}), fullfile(outputDir,'cutline_along_y.csv'));
clear cleanup
end
