function Build_Plain_Heater_Geom()
%BUILD_PLAIN_HEATER_GEOM Geometry-only build + domain/boundary ID prospecting.
%
% Plain (un-patterned) heater benchmark: the same pad / lead / substrate / silica
% stack as the inverse-design Step3 model, but the patterned pixel array is
% replaced by ONE solid 50x50 um low-dope heater block. NO materials / physics /
% mesh are set here -- this only builds the geometry and prints the centroid +
% bounding box of every DOMAIN and the centroid of every BOUNDARY, so the IDs can
% be mapped by hand before writing the real model (IDs differ completely from the
% patterned Step3 and must NOT be reused).

import com.comsol.model.*
import com.comsol.model.util.*

mliPath = comsol_mli_path();
if exist(mliPath, 'dir') == 7, addpath(mliPath); end
try, mphstart(comsol_server_port()); fprintf('Connected to COMSOL Server on port 2036.\n'); catch, end

% --- Geometry params: replicate the inverse-design (44x44) Step3 stack exactly,
%     so the thermal environment matches; only the heater interior is solid. ---
Total_W = 50000; Total_H = 50000;          % 50 um x 50 um heater footprint
array_x_left = -4000; array_x_right = array_x_left + Total_W;   % -4000 .. 46000
array_y_bot  = -4000; array_y_top  = array_y_bot  + Total_H;    % -4000 .. 46000
centx_val = array_x_left + Total_W/2;       % 21000
centy_val = array_y_bot  + Total_H/2;       % 21000

lead_x_overhang = 10000; lead_height = 10000;
pad_height_top = 40000; pad_height_bot_main = 35000; pad_height_bot_transition = 5000;
pad_width = Total_W + 2*lead_x_overhang;    % 70000
sq1_pos_y = array_y_top + lead_height + pad_height_top/2;
sq2_pos_y = array_y_bot - lead_height - pad_height_bot_transition - pad_height_bot_main/2;
sq4_pos_y = array_y_bot - lead_height - pad_height_bot_transition/2;
pol2_table = [array_x_left,                    array_y_top;
              array_x_right,                   array_y_top;
              array_x_right + lead_x_overhang, array_y_top + lead_height;
              array_x_left  - lead_x_overhang, array_y_top + lead_height];
pol1_table = [array_x_left,                    array_y_bot;
              array_x_left  - lead_x_overhang, array_y_bot - lead_height;
              array_x_right + lead_x_overhang, array_y_bot - lead_height;
              array_x_right,                   array_y_bot];

try, ModelUtil.remove('Plain_Geom'); catch, end
model = ModelUtil.create('Plain_Geom');
model.modelPath(pwd);
model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 3);
geom = model.component('comp1').geom('geom1');
geom.lengthUnit('nm');
model.param.set('L_m', '300000[nm]');
model.param.set('centx', sprintf('%d[nm]', centx_val));
model.param.set('centy', sprintf('%d[nm]', centy_val));

geom.create('wp1', 'WorkPlane'); geom.feature('wp1').set('unite', true);
wp1 = geom.feature('wp1').geom;

% Pads / leads (same order as Step3 so the base stack matches)
wp1.create('sq1', 'Rectangle'); wp1.feature('sq1').set('size', [pad_width pad_height_top]); wp1.feature('sq1').set('base', 'center'); wp1.feature('sq1').set('pos', {'centx' num2str(sq1_pos_y)});
wp1.create('pol2', 'Polygon'); wp1.feature('pol2').set('source', 'table'); wp1.feature('pol2').set('table', pol2_table);
wp1.create('sq2', 'Rectangle'); wp1.feature('sq2').set('size', [pad_width pad_height_bot_main]); wp1.feature('sq2').set('base', 'center'); wp1.feature('sq2').set('pos', {'centx' num2str(sq2_pos_y)});
wp1.create('sq4', 'Rectangle'); wp1.feature('sq4').set('size', [pad_width pad_height_bot_transition]); wp1.feature('sq4').set('base', 'center'); wp1.feature('sq4').set('pos', {'centx' num2str(sq4_pos_y)});
wp1.create('pol1', 'Polygon'); wp1.feature('pol1').set('source', 'table'); wp1.feature('pol1').set('table', pol1_table);
wp1.create('sq3', 'Square'); wp1.feature('sq3').set('size', 'L_m'); wp1.feature('sq3').set('base', 'center'); wp1.feature('sq3').set('pos', {'centx' 'centy'});

% THE heater: one solid 50x50 um block (replaces the patterned pixel array)
wp1.create('heater', 'Rectangle'); wp1.feature('heater').set('size', [Total_W Total_H]); wp1.feature('heater').set('base', 'center'); wp1.feature('heater').set('pos', {'centx' 'centy'});
geom.run('wp1');

% Same z-stack extrudes as Step3
geom.feature.create('ext1', 'Extrude'); geom.feature('ext1').set('workplane', 'wp1'); geom.feature('ext1').selection('input').set({'wp1'}); geom.feature('ext1').setIndex('distance', 220, 0); geom.run('ext1');
geom.create('wp2', 'WorkPlane'); geom.feature('wp2').set('unite', true); wp2 = geom.feature('wp2').geom; wp2.create('sq1', 'Square'); wp2.feature('sq1').set('size', 'L_m'); wp2.feature('sq1').set('base', 'center'); wp2.feature('sq1').set('pos', {'centx' 'centy'}); geom.run('wp2');
geom.feature.create('ext2', 'Extrude'); geom.feature('ext2').set('workplane', 'wp2'); geom.feature('ext2').selection('input').set({'wp2'}); geom.feature('ext2').setIndex('distance', -2000, 0); geom.run('ext2');
geom.create('wp3', 'WorkPlane'); geom.feature('wp3').set('unite', true); geom.feature('wp3').set('quickz', 220); wp3 = geom.feature('wp3').geom; wp3.create('sq1', 'Square'); wp3.feature('sq1').set('size', 'L_m'); wp3.feature('sq1').set('base', 'center'); wp3.feature('sq1').set('pos', {'centx' 'centy'}); geom.run('wp3');
geom.feature.create('ext3', 'Extrude'); geom.feature('ext3').set('workplane', 'wp3'); geom.feature('ext3').selection('input').set({'wp3'}); geom.feature('ext3').setIndex('distance', 10, 0); geom.run('ext3');
geom.create('wp4', 'WorkPlane'); geom.feature('wp4').set('quickz', -2000); wp4 = geom.feature('wp4').geom; wp4.create('sq1', 'Square'); wp4.feature('sq1').set('size', 'L_m'); wp4.feature('sq1').set('base', 'center'); wp4.feature('sq1').set('pos', {'centx' 'centy'}); geom.run('wp4');
geom.feature.create('ext4', 'Extrude'); geom.feature('ext4').set('workplane', 'wp4'); geom.feature('ext4').selection('input').set({'wp4'}); geom.feature('ext4').setIndex('distance', -20000, 0); geom.run('ext4');
geom.run;

nD = geom.getNDomains(); nB = geom.getNBoundaries();
fprintf('\n===== Geometry stats: %d domains, %d boundaries =====\n', nD, nB);
fprintf('Expected: heater block z[0,220] xy=[%.0f,%.0f]x[%.0f,%.0f]; top pad y~%.0f; bot pad y~%.0f\n', ...
    array_x_left, array_x_right, array_y_bot, array_y_top, sq1_pos_y, sq2_pos_y);
fprintf('Stack z: substrate[-22000,-2000], silica[-2000,0], doped[0,220], thin-silica[220,230]\n');

fprintf('\n----- DOMAINS (centroid + bbox) -----\n');
for id = 1:nD
    try
        c = mphgetcoords(model, 'geom1', 'domain', id);
        cm = mean(c, 2);
        fprintf(' dom %2d: centroid(%7.0f,%7.0f,%7.0f)  x[%7.0f,%7.0f] y[%7.0f,%7.0f] z[%7.0f,%7.0f]\n', ...
            id, cm(1), cm(2), cm(3), min(c(1,:)), max(c(1,:)), min(c(2,:)), max(c(2,:)), min(c(3,:)), max(c(3,:)));
    catch ME
        fprintf(' dom %2d: error %s\n', id, ME.message);
    end
end

fprintf('\n----- BOUNDARIES of interest (z-plane candidates) -----\n');
for id = 1:nB
    try
        c = mphgetcoords(model, 'geom1', 'boundary', id);
        cm = mean(c, 2); z = cm(3); y = cm(2);
        tag = '';
        if abs(z-(-22000))<5, tag='CAND temp1 (substrate bottom)';
        elseif abs(z-230)<5, tag='CAND hf1 (thin-silica top)';
        elseif abs(z-220)<5 && y> array_y_top+lead_height-1, tag='CAND pot1 (top pad face)';
        elseif abs(z-220)<5 && y< array_y_bot-lead_height+1, tag='CAND gnd1 (bottom pad face)'; end
        if ~isempty(tag)
            fprintf(' bnd %3d: (%7.0f,%7.0f,%7.0f)  %s\n', id, cm(1), cm(2), cm(3), tag);
        end
    catch
    end
end

mphsave(model, 'Plain_Heater_Geom.mph');
fprintf('\nSaved Plain_Heater_Geom.mph. Done.\n');
end
