function Step1a_Generate_ICM_Model_Hetero(k_SOI_val, k_Hdoped_val, k_Ldoped_val)
% =========================================================================
% Step 1a: build the thermal response model on a nonuniform rectangular grid.
% Six virtual groups: center, top, bottom, side edges, top corners, bottom corners.
% Core-cell and peripheral-cell dimensions are specified separately.
% =========================================================================

import com.comsol.model.*
import com.comsol.model.util.*

% --- 0. Grid dimensions ---
Nx_core = 31;    % Number of core columns
Ny_core = 21;    % Number of core rows

W_core = 2000;   % Core-cell width (nm)
H_core = 2000;   % Core-cell height (nm)

W_edge = 4000;   % Peripheral-cell width on the left and right (nm)
H_edge = 4000;   % Peripheral-cell height on the top and bottom (nm)

% Add one peripheral cell at each side of the core grid.
Nx_new = Nx_core + 2; 
Ny_new = Ny_core + 2; 

ModelName = sprintf('ICM_Hetero_6G_%dx%d', Nx_new, Ny_new);
try model = ModelUtil.create(ModelName); catch, model = ModelUtil.model(ModelName); end
model.modelPath(pwd); model.hist.disable;

% --- 1. Model parameters and material functions ---
fprintf('Configure parameters and analytical functions.\n');
centx_val = Nx_core * W_core / 2;
centy_val = Ny_core * H_core / 2;

model.param.set('L_m', '300000[nm]'); 
model.param.set('centx', sprintf('%d[nm]', centx_val)); 
model.param.set('centy', sprintf('%d[nm]', centy_val)); 

model.param.set('P_probe', '0.001000[W]'); 
model.param.set('case_id', '1');
model.param.set('Target_temperature', '900 [K]');

% Use the material parameters supplied by the calling workflow.
model.param.set('k_SOI', sprintf('%f [W/(m*K)]', k_SOI_val));
model.param.set('k_Hdoped', num2str(k_Hdoped_val));
model.param.set('k_Ldoped', num2str(k_Ldoped_val));

% Define the step function, silicon heat capacity, and thermal conductivity.
model.func.create('step1', 'Step'); model.func('step1').set('smooth', 10);

model.func.create('an1', 'Analytic'); model.func('an1').set('funcname', 'Cp_Si');
model.func('an1').set('expr', '1000/28.09*(22.82+3.9*(x/1000)-0.0829*(x/1000)^2+0.0421*(x/1000)^3-0.354*(x/1000)^(-2))');
model.func('an1').set('argunit', {'K'});
model.func('an1').set('plotargs', {'x' '200' '1200'});

model.func.create('an2', 'Analytic'); model.func('an2').set('funcname', 'k_Si');
model.func('an2').set('expr', 'k_SOI*293.15/x');
model.func('an2').set('argunit', {'K'});
model.func('an2').set('plotargs', {'x' '200' '1200'});

% --- 2. Nonuniform grid geometry ---
fprintf('Build the %dx%d nonuniform grid.\n', Nx_new, Ny_new);
model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 3);
geom = model.component('comp1').geom('geom1');
geom.lengthUnit('nm');
geom.create('wp1', 'WorkPlane'); geom.feature('wp1').set('unite', true); 
wp1 = geom.feature('wp1').geom;

% Construct the cell-dimension arrays.
W_arr = ones(1, Nx_new) * W_core; 
W_arr(1) = W_edge; W_arr(end) = W_edge; 

H_arr = ones(1, Ny_new) * H_core;
H_arr(1) = H_edge; H_arr(end) = H_edge;

X_bl = zeros(1, Nx_new); X_bl(1) = -W_arr(1); 
for i = 2:Nx_new, X_bl(i) = X_bl(i-1) + W_arr(i-1); end
Y_bl = zeros(1, Ny_new); Y_bl(1) = -H_arr(1); 
for j = 2:Ny_new, Y_bl(j) = Y_bl(j-1) + H_arr(j-1); end

% Position the pads relative to the full grid extent.
pad_width = sum(W_arr) + 20000; 
pad_height = 50000; 
pad_offset = pad_height / 2;
pos_sq1_y = Ny_core * H_core + H_edge + pad_offset; 
pos_sq2_y = -H_edge - pad_offset;

wp1.create('sq1', 'Rectangle'); wp1.feature('sq1').set('size', [pad_width pad_height]); wp1.feature('sq1').set('base', 'center'); wp1.feature('sq1').set('pos', {'centx' num2str(pos_sq1_y)});
wp1.create('sq2', 'Rectangle'); wp1.feature('sq2').set('size', [pad_width pad_height]); wp1.feature('sq2').set('base', 'center'); wp1.feature('sq2').set('pos', {'centx' num2str(pos_sq2_y)});
wp1.create('sq3', 'Square'); wp1.feature('sq3').set('size', 'L_m'); wp1.feature('sq3').set('base', 'center'); wp1.feature('sq3').set('pos', {'centx' 'centy'});

for i = 1:Nx_new
    for j = 1:Ny_new
        xpos = X_bl(i) + W_arr(i)/2;
        ypos = Y_bl(j) + H_arr(j)/2;
        
        bg_tag = sprintf('bg_%d_%d', i, j);
        wp1.create(bg_tag, 'Rectangle'); wp1.feature(bg_tag).set('size', [W_arr(i) H_arr(j)]); wp1.feature(bg_tag).set('base', 'center'); wp1.feature(bg_tag).set('pos', [xpos ypos]);

        current_wire_width = W_arr(i) / 2;
        r_tag = sprintf('r_%d_%d', i, j);
        wp1.create(r_tag, 'Rectangle'); wp1.feature(r_tag).set('size', [current_wire_width H_arr(j)]); wp1.feature(r_tag).set('base', 'center'); wp1.feature(r_tag).set('pos', [xpos ypos]);
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

% --- 3. Material domains and properties ---
fprintf('Assign thermal-conductivity properties.\n');
high_dope = []; no_dope = []; count = 7; 
domain_id_map = zeros(Nx_new, Ny_new);
for i = 1:Nx_new   
    temp1 = count : (count + Ny_new - 1);
    temp2 = (count + Ny_new) : (count + 2*Ny_new - 1);
    temp3 = (count + 2*Ny_new) : (count + 3*Ny_new - 1);
    count = count + 3*Ny_new;
    high_dope = [high_dope temp2]; no_dope = [no_dope temp1 temp3]; 
    domain_id_map(i, :) = temp2;
end

% mat1: highdope_Si
model.component('comp1').material.create('mat1', 'Common'); model.component('comp1').material('mat1').label('highdope_Si');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'k_Hdoped*k_Si(293.15[K])' '0' '0' '0' 'k_Hdoped*k_Si(293.15[K])' '0' '0' '0' 'k_Hdoped*k_Si(293.15[K])'});
model.component('comp1').material('mat1').propertyGroup('def').set('density', '2329');
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', '700');
model.component('comp1').material('mat1').selection.set([5 6]);

% mat4: lowdope_Si
model.component('comp1').material.create('mat4', 'Common'); model.component('comp1').material('mat4').label('lowdope_Si');
model.component('comp1').material('mat4').propertyGroup('def').set('thermalconductivity', {'k_Ldoped*k_Si(Target_temperature)' '0' '0' '0' 'k_Ldoped*k_Si(Target_temperature)' '0' '0' '0' 'k_Ldoped*k_Si(Target_temperature)'});
model.component('comp1').material('mat4').propertyGroup('def').set('density', '2329');
model.component('comp1').material('mat4').propertyGroup('def').set('heatcapacity', '700');
model.component('comp1').material('mat4').selection.set(high_dope);

% mat5: Silicon_substrate
model.component('comp1').material.create('mat5', 'Common'); model.component('comp1').material('mat5').label('Silicon_substrate');
model.component('comp1').material('mat5').propertyGroup('def').set('thermalconductivity', {'130' '0' '0' '0' '130' '0' '0' '0' '130'});
model.component('comp1').material('mat5').propertyGroup('def').set('density', '2329');
model.component('comp1').material('mat5').propertyGroup('def').set('heatcapacity', '700');
model.component('comp1').material('mat5').selection.set([1]);

% mat2: undoped_Si
model.component('comp1').material.create('mat2', 'Common'); model.component('comp1').material('mat2').label('undoped_Si');
model.component('comp1').material('mat2').propertyGroup('def').set('thermalconductivity', {'k_Si(Target_temperature)' '0' '0' '0' 'k_Si(Target_temperature)' '0' '0' '0' 'k_Si(Target_temperature)'});
model.component('comp1').material('mat2').propertyGroup('def').set('density', '2329');
model.component('comp1').material('mat2').propertyGroup('def').set('heatcapacity', '700');
model.component('comp1').material('mat2').selection.set(no_dope);

% mat6: SOI_Si
model.component('comp1').material.create('mat6', 'Common'); model.component('comp1').material('mat6').label('SOI_Si');
model.component('comp1').material('mat6').propertyGroup('def').set('thermalconductivity', {'k_Si(293.15[K])' '0' '0' '0' 'k_Si(293.15[K])' '0' '0' '0' 'k_Si(293.15[K])'});
model.component('comp1').material('mat6').propertyGroup('def').set('density', '2329');
model.component('comp1').material('mat6').propertyGroup('def').set('heatcapacity', '700');
model.component('comp1').material('mat6').selection.set([3]);

% mat3: Silica glass
model.component('comp1').material.create('mat3', 'Common'); model.component('comp1').material('mat3').label('Silica glass');
model.component('comp1').material('mat3').propertyGroup('def').set('thermalconductivity', {'1.38' '0' '0' '0' '1.38' '0' '0' '0' '1.38'});
model.component('comp1').material('mat3').propertyGroup('def').set('density', '2203');
model.component('comp1').material('mat3').propertyGroup('def').set('heatcapacity', '703');
model.component('comp1').material('mat3').selection.set([2 4]);

% --- 4. Heat-transfer physics and virtual sources ---
model.component('comp1').physics.create('ht', 'HeatTransfer', 'geom1');
masks = generate_group_masks(Nx_new, Ny_new);
Area_matrix = (W_arr' * H_arr) / (1000^2);

% Create one source for each of the six groups.
for k = 1:6
    sel_tag = sprintf('sel_group_%d', k);
    model.component('comp1').selection.create(sel_tag, 'Explicit');
    group_mask = masks{k};
    model.component('comp1').selection(sel_tag).set(domain_id_map(group_mask));
    
    hs_tag = sprintf('hs_%d', k);
    model.component('comp1').physics('ht').create(hs_tag, 'HeatSource', 3);
    model.component('comp1').physics('ht').feature(hs_tag).selection.named(sel_tag);
    
    equivalent_pixels = sum(Area_matrix(group_mask));
    model.component('comp1').physics('ht').feature(hs_tag).set('heatSourceType', 'TotalPower');
    logic_expr = sprintf('P_probe * %f * (case_id == %d) * step1(t/(1[ns]))', equivalent_pixels, k);
    model.component('comp1').physics('ht').feature(hs_tag).set('Ptot', logic_expr);
end

model.component('comp1').physics('ht').create('temp1', 'TemperatureBoundary', 2);
model.component('comp1').physics('ht').feature('temp1').selection.set([3]); 
model.component('comp1').physics('ht').feature('temp1').set('T0', '293.15[K]');
model.component('comp1').physics('ht').create('hf1', 'HeatFluxBoundary', 2);
model.component('comp1').physics('ht').feature('hf1').selection.set([13]); 
model.component('comp1').physics('ht').feature('hf1').set('HeatFluxType', 'ConvectiveHeatFlux');
model.component('comp1').physics('ht').feature('hf1').set('h', 10);

% --- 5. Mesh ---
fprintf('Generate the mesh.\n');
model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').create('ftet1', 'FreeTet');
model.component('comp1').mesh('mesh1').feature('ftet1').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('ftet1').selection.set([high_dope no_dope]);
model.component('comp1').mesh('mesh1').feature('ftet1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('custom', true);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmax', H_edge / 5); 
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmin', H_edge / 200);   

model.component('comp1').mesh('mesh1').create('ftet2', 'FreeTet');
model.component('comp1').mesh('mesh1').feature('ftet2').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('ftet2').selection.set([5 6]); 
model.component('comp1').mesh('mesh1').feature('ftet2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hauto', 6); 

model.component('comp1').mesh('mesh1').create('ftet3', 'FreeTet');
model.component('comp1').mesh('mesh1').feature('ftet3').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('ftet3').selection.set([3]); 
model.component('comp1').mesh('mesh1').feature('ftet3').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftet3').feature('size1').set('hauto', 8); 

model.component('comp1').mesh('mesh1').create('swe1', 'Sweep');
model.component('comp1').mesh('mesh1').run;

% --- 6. Results and model output ---
y_mid = -H_edge / 2; 
model.result.dataset.create('cln1', 'CutLine3D');
model.result.dataset('cln1').set('genpoints', [-50000, y_mid, 110; pad_width+50000, y_mid, 110]);
model.result.create('pg1', 'PlotGroup3D'); model.result('pg1').label('Temperature (ht)');
model.result('pg1').create('vol1', 'Volume'); model.result('pg1').feature('vol1').set('expr', 'T');
model.result('pg1').feature('vol1').set('colortable', 'HeatCameraLight');
model.result.create('pg2', 'PlotGroup1D'); model.result('pg2').label('1D Temperature Profile');
model.result('pg2').set('data', 'cln1'); model.result('pg2').create('lngr1', 'LineGraph');
model.result('pg2').feature('lngr1').set('data', 'cln1'); model.result('pg2').feature('lngr1').set('expr', 'T');

model.study.create('std1');
model.study('std1').create('time', 'Transient'); 
model.study('std1').feature('time').set('tlist', 'range(0, 1.0e-6, 1.3e-5)');
model.study('std1').create('param', 'Parametric');
model.study('std1').feature('param').set('pname', {'case_id'});
% Sweep the six virtual heat-source groups.
model.study('std1').feature('param').set('plistarr', {'1 2 3 4 5 6'});

model.sol.create('sol1');
model.sol('sol1').attach('std1');
model.sol('sol1').createAutoSequence('std1');
model.sol('sol1').feature('t1').set('rtol', 0.001);
model.sol('sol1').feature('t1').set('tstepsbdf', 'strict');
model.sol('sol1').feature('t1').set('maxstepconstraintbdf', 'const');
model.sol('sol1').feature('t1').set('maxstepbdf', '0.5E-6');

% Solve the thermal response model.
fprintf('Solve the thermal response model.\n');
model.study('std1').run;

mphsave(model, ModelName);
fprintf('Thermal response model solved and saved: %s.mph\n', ModelName);
end

% --- Six-group masks ---
function masks = generate_group_masks(Nx, Ny)
    masks = cell(1, 6);
    
    % 1. Center
    mask_c = false(Nx, Ny); 
    mask_c(2:Nx-1, 2:Ny-1) = true; 
    masks{1} = mask_c;
    
    % 2. Top edge
    mask_top = false(Nx, Ny); 
    mask_top(2:Nx-1, Ny) = true; 
    masks{2} = mask_top;
    
    % 3. Bottom edge
    mask_bottom = false(Nx, Ny); 
    mask_bottom(2:Nx-1, 1) = true; 
    masks{3} = mask_bottom;
    
    % 4. Left and right edges
    mask_lr = false(Nx, Ny); 
    mask_lr([1, Nx], 2:Ny-1) = true; 
    masks{4} = mask_lr;
    
    % 5. Top corners
    mask_top_corner = false(Nx, Ny); 
    mask_top_corner([1, Nx], Ny) = true; 
    masks{5} = mask_top_corner;
    
    % 6. Bottom corners
    mask_bot_corner = false(Nx, Ny); 
    mask_bot_corner([1, Nx], 1) = true; 
    masks{6} = mask_bot_corner;
end
