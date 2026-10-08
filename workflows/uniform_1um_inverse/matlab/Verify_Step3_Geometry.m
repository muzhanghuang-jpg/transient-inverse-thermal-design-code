function Verify_Step3_Geometry()
%VERIFY_STEP3_GEOMETRY Build Step3 geometry only and report boundary centroids.
%
% Replicates the Step3 wp1/extrude geometry without mesh / study / solve.
% Prints the centroid of every boundary, plus a direct readout of the four
% hard-coded IDs (gnd1=21, pot1=35, temp1=3, hf1=13) so we can verify whether
% the diode-as-highdope reassignment and the X-shrink to a square heater
% leave the boundary numbering intact.

import com.comsol.model.*
import com.comsol.model.util.*

mliPath = comsol_mli_path();
if exist(mliPath, 'dir') == 7
    addpath(mliPath);
end
try
    mphstart(comsol_server_port());
    fprintf('Connected to COMSOL Server on port 2036.\n');
catch
    fprintf('mphstart already connected or running inside MATLAB-LiveLink session.\n');
end

% Array sizing must match Step1a/Step3.
Nx_core = 21; Ny_core = 21;
W_core = 2000; H_core = 2000;
W_edge = 4000; H_edge = 4000;
Nx = Nx_core + 2; Ny = Ny_core + 2;

W_arr = ones(1, Nx) * W_core; W_arr(1) = W_edge; W_arr(end) = W_edge;
H_arr = ones(1, Ny) * H_core; H_arr(1) = H_edge; H_arr(end) = H_edge;

X_bl = zeros(1, Nx); X_bl(1) = -W_arr(1);
for i = 2:Nx, X_bl(i) = X_bl(i-1) + W_arr(i-1); end
Y_bl = zeros(1, Ny); Y_bl(1) = -H_arr(1);
for j = 2:Ny, Y_bl(j) = Y_bl(j-1) + H_arr(j-1); end
Total_W = sum(W_arr); Total_H = sum(H_arr);

% Same parameter block as Step3.
lead_x_overhang = 10000;
lead_height = 10000;
pad_height_top = 40000;
pad_height_bot_main = 35000;
pad_height_bot_transition = 5000;
pad_width = Total_W + 2 * lead_x_overhang;
array_x_left = X_bl(1);
array_x_right = X_bl(1) + Total_W;
array_y_bot = Y_bl(1);
array_y_top = Y_bl(1) + Total_H;
sq1_pos_y = array_y_top + lead_height + pad_height_top / 2;
sq2_pos_y = array_y_bot - lead_height - pad_height_bot_transition - pad_height_bot_main / 2;
sq4_pos_y = array_y_bot - lead_height - pad_height_bot_transition / 2;
pol2_table = [array_x_left,                    array_y_top;
              array_x_right,                   array_y_top;
              array_x_right + lead_x_overhang, array_y_top + lead_height;
              array_x_left  - lead_x_overhang, array_y_top + lead_height];
pol1_table = [array_x_left,                    array_y_bot;
              array_x_left  - lead_x_overhang, array_y_bot - lead_height;
              array_x_right + lead_x_overhang, array_y_bot - lead_height;
              array_x_right,                   array_y_bot];
centx_val = X_bl(1) + Total_W / 2;
centy_val = Y_bl(1) + Total_H / 2;

% Dummy uniform wire width well inside the side-gap limit.
width_uniform = 1000;

try ModelUtil.remove('Verify_Geom'); catch, end
model = ModelUtil.create('Verify_Geom');
model.modelPath(pwd);
model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 3);
geom = model.component('comp1').geom('geom1');
geom.lengthUnit('nm');

model.param.set('L_m', '300000[nm]');
model.param.set('centx', sprintf('%d[nm]', centx_val));
model.param.set('centy', sprintf('%d[nm]', centy_val));

geom.create('wp1', 'WorkPlane');
geom.feature('wp1').set('unite', true);
wp1 = geom.feature('wp1').geom;

wp1.create('sq1', 'Rectangle');
wp1.feature('sq1').set('size', [pad_width pad_height_top]);
wp1.feature('sq1').set('base', 'center');
wp1.feature('sq1').set('pos', {'centx' num2str(sq1_pos_y)});

wp1.create('pol2', 'Polygon');
wp1.feature('pol2').set('source', 'table');
wp1.feature('pol2').set('table', pol2_table);

wp1.create('sq2', 'Rectangle');
wp1.feature('sq2').set('size', [pad_width pad_height_bot_main]);
wp1.feature('sq2').set('base', 'center');
wp1.feature('sq2').set('pos', {'centx' num2str(sq2_pos_y)});

wp1.create('sq4', 'Rectangle');
wp1.feature('sq4').label('Rectangle 2.1');
wp1.feature('sq4').set('size', [pad_width pad_height_bot_transition]);
wp1.feature('sq4').set('base', 'center');
wp1.feature('sq4').set('pos', {'centx' num2str(sq4_pos_y)});

wp1.create('pol1', 'Polygon');
wp1.feature('pol1').set('source', 'table');
wp1.feature('pol1').set('table', pol1_table);

wp1.create('sq3', 'Square');
wp1.feature('sq3').set('size', 'L_m');
wp1.feature('sq3').set('base', 'center');
wp1.feature('sq3').set('pos', {'centx' 'centy'});

for i = 1:Nx
    for j = 1:Ny
        xpos = X_bl(i) + W_arr(i)/2;
        ypos = Y_bl(j) + H_arr(j)/2;
        bg_tag = sprintf('bg_%d_%d', i, j);
        wp1.create(bg_tag, 'Rectangle');
        wp1.feature(bg_tag).set('size', [W_arr(i) H_arr(j)]);
        wp1.feature(bg_tag).set('base', 'center');
        wp1.feature(bg_tag).set('pos', [xpos ypos]);
        r_tag = sprintf('r_%d_%d', i, j);
        wp1.create(r_tag, 'Rectangle');
        wp1.feature(r_tag).set('size', [width_uniform H_arr(j)]);
        wp1.feature(r_tag).set('base', 'center');
        wp1.feature(r_tag).set('pos', [xpos ypos]);
    end
end
geom.run('wp1');

geom.feature.create('ext1', 'Extrude'); geom.feature('ext1').set('workplane', 'wp1'); geom.feature('ext1').selection('input').set({'wp1'}); geom.feature('ext1').setIndex('distance', 220, 0); geom.run('ext1');
geom.create('wp2', 'WorkPlane'); geom.feature('wp2').set('unite', true); wp2 = geom.feature('wp2').geom; wp2.create('sq1', 'Square'); wp2.feature('sq1').set('size', 'L_m'); wp2.feature('sq1').set('base', 'center'); wp2.feature('sq1').set('pos', {'centx' 'centy'}); geom.run('wp2');
geom.feature.create('ext2', 'Extrude'); geom.feature('ext2').set('workplane', 'wp2'); geom.feature('ext2').selection('input').set({'wp2'}); geom.feature('ext2').setIndex('distance', -2000, 0); geom.run('ext2');
geom.create('wp3', 'WorkPlane'); geom.feature('wp3').set('unite', true); geom.feature('wp3').set('quickz', 220); wp3 = geom.feature('wp3').geom; wp3.create('sq1', 'Square'); wp3.feature('sq1').set('size', 'L_m'); wp3.feature('sq1').set('base', 'center'); wp3.feature('sq1').set('pos', {'centx' 'centy'}); geom.run('wp3');
geom.feature.create('ext3', 'Extrude'); geom.feature('ext3').set('workplane', 'wp3'); geom.feature('ext3').selection('input').set({'wp3'}); geom.feature('ext3').setIndex('distance', 10, 0); geom.run('ext3');
geom.create('wp4', 'WorkPlane'); geom.feature('wp4').set('quickz', -2000); wp4 = geom.feature('wp4').geom; wp4.create('sq1', 'Square'); wp4.feature('sq1').set('size', 'L_m'); wp4.feature('sq1').set('base', 'center'); wp4.feature('sq1').set('pos', {'centx' 'centy'}); geom.run('wp4');
geom.feature.create('ext4', 'Extrude'); geom.feature('ext4').set('workplane', 'wp4'); geom.feature('ext4').selection('input').set({'wp4'}); geom.feature('ext4').setIndex('distance', -20000, 0); geom.run('ext4');
geom.run;

nDomains = geom.getNDomains();
nBoundaries = geom.getNBoundaries();
fprintf('\n========== Geometry stats ==========\n');
fprintf('Total domains   : %d\n', nDomains);
fprintf('Total boundaries: %d\n', nBoundaries);
fprintf('Expected pad y-ranges (nm):\n');
fprintf('  top pad sq1 (z=220 face): y in [%.0f, %.0f]\n', sq1_pos_y - pad_height_top/2, sq1_pos_y + pad_height_top/2);
fprintf('  bot pad sq2 (z=220 face): y in [%.0f, %.0f]\n', sq2_pos_y - pad_height_bot_main/2, sq2_pos_y + pad_height_bot_main/2);
fprintf('  bot trans sq4 (z=220 face): y in [%.0f, %.0f]\n', sq4_pos_y - pad_height_bot_transition/2, sq4_pos_y + pad_height_bot_transition/2);
fprintf('  array y-range: [%.0f, %.0f], x-range: [%.0f, %.0f]\n', array_y_bot, array_y_top, array_x_left, array_x_right);

ids_to_check = [21, 35, 3, 13];
labels = {'gnd1 (Ground, expect z=220 over bottom pad)', ...
          'pot1 (Potential, expect z=220 over top pad)', ...
          'temp1 (293.15K, expect z=-22000 substrate bottom)', ...
          'hf1 (convection, expect z=230 silica thin top)'};
fprintf('\n========== Hard-coded boundary IDs ==========\n');
for k = 1:numel(ids_to_check)
    id = ids_to_check(k);
    if id > nBoundaries
        fprintf('  ID %2d (%s): OUT OF RANGE (max %d)\n', id, labels{k}, nBoundaries);
        continue;
    end
    try
        c = mphgetcoords(model, 'geom1', 'boundary', id);
        cmean = mean(c, 2);
        fprintf('  ID %2d (%s):\n', id, labels{k});
        fprintf('       centroid (x, y, z) = (%.0f, %.0f, %.0f) nm\n', cmean(1), cmean(2), cmean(3));
    catch ME
        fprintf('  ID %2d (%s): mphgetcoords error %s\n', id, labels{k}, ME.message);
    end
end

fprintf('\n========== All boundaries on key z-planes ==========\n');
fprintf('(Looking for candidates of each expected face)\n');
for id = 1:nBoundaries
    try
        c = mphgetcoords(model, 'geom1', 'boundary', id);
        cmean = mean(c, 2);
        z = cmean(3); x = cmean(1); y = cmean(2);
        tag = '';
        if abs(z - (-22000)) < 5
            tag = 'CANDIDATE temp1 (substrate bottom z=-22000)';
        elseif abs(z - 230) < 5
            tag = 'CANDIDATE hf1 (silica thin top z=230)';
        elseif abs(z - 220) < 5 && y > array_y_top + lead_height - 1
            tag = 'CANDIDATE pot1 (top pad/lead face z=220)';
        elseif abs(z - 220) < 5 && y < array_y_bot - lead_height + 1
            tag = 'CANDIDATE gnd1 (bottom pad/trans/lead face z=220)';
        end
        if ~isempty(tag)
            fprintf('  boundary %3d  (x,y,z)=(%6.0f,%6.0f,%6.0f) nm  %s\n', id, x, y, z, tag);
        end
    catch
    end
end

mphsave(model, 'Verify_Step3_Geometry.mph');
fprintf('\nSaved geometry to Verify_Step3_Geometry.mph\n');
fprintf('Done.\n');
end
