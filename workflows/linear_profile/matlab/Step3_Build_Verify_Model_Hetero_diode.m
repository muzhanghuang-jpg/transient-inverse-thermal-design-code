function Step3_Build_Verify_Model_Hetero_diode(iter, k_SOI_val, k_Hdoped_val, k_Ldoped_val)
% =========================================================================
% Step 3: build the transient COMSOL verification model for the heterogeneous array.
% 
% Notes:
% 1. Load heterogeneous array dimensions and optimized wire widths.
% 2. Generate geometry with the same bg_i_j and r_i_j loop structure as Step1a.
% 3. Preserve domain numbering through the count = count + 3*Ny progression.
% =========================================================================

import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.showProgress(true);

% --- 0. Load heterogeneous grid and optimization data ---
matFiles_G = dir('G_Matrix_Hetero_*.mat');
if isempty(matFiles_G)
    error('Base array file not found. Confirm that G_Matrix_Hetero_*.mat exists.');
end
[~, latest_idx] = max([matFiles_G.datenum]);
G_Data_File = matFiles_G(latest_idx).name;
load(G_Data_File, 'W_arr', 'H_arr'); 
Nx = length(W_arr);
Ny = length(H_arr);

% Choose the initial optimization result or the previous correction by iter.
if iter == 0
    DataFile = sprintf('Initial_Optimization_Solution_Hetero_%dx%d.mat', Nx, Ny);
else
    DataFile = sprintf('Correction_Step_Latest_Hetero_%dx%d_diode_%d.mat', Nx, Ny, iter);
end

OutFile = sprintf('Verify_Model_Hetero_%dx%d_diode_%d.mph', Nx, Ny, iter);

tempData = load(DataFile);

if isfield(tempData, 'width_map_New')
    width_map = tempData.width_map_New;
elseif isfield(tempData, 'width_map_final')
    width_map = tempData.width_map_final;
else
    error('Wire-width data was not found in %s.', DataFile);
end

Voltage_Use = tempData.Voltage_Use;
Target_Time = tempData.Target_Time;
targetSettings = get_thermal_target_settings();
Time_Step = targetSettings.Time_Step_s;
geometrySettings = get_geometry_constraint_settings();
Min_Wire_Width_nm = geometrySettings.Min_Wire_Width_nm;
Min_Side_Gap_nm = geometrySettings.Min_Side_Gap_nm;


fprintf('Loaded heterogeneous array data: Nx = %d, Ny = %d\n', Nx, Ny);
fprintf('Operating voltage: %.1f V, target time: %.2f us\n', Voltage_Use, Target_Time*1e6);

% --- 0.5 Coordinate precomputation ---
X_bl = zeros(1, Nx); X_bl(1) = -W_arr(1); 
for i = 2:Nx, X_bl(i) = X_bl(i-1) + W_arr(i-1); end

Y_bl = zeros(1, Ny); Y_bl(1) = -H_arr(1); 
for j = 2:Ny, Y_bl(j) = Y_bl(j-1) + H_arr(j-1); end

Total_W = sum(W_arr);
Total_H = sum(H_arr);
Time_List = sprintf('range(0, %.16g, %.16g)', Time_Step, Target_Time);

% Pad / lead geometry parameters, derived from array Total_W / Total_H so the
% wp1 shapes scale with the heater. The sq4 rectangle (formerly Diode) is kept
% geometrically and reassigned to highdope_Si in the materials block.
lead_x_overhang           = 10000;
lead_height               = 10000;
pad_height_top            = 40000;
pad_height_bot_main       = 35000;
pad_height_bot_transition = 5000;
pad_width                 = Total_W + 2 * lead_x_overhang;
array_x_left              = X_bl(1);
array_x_right             = X_bl(1) + Total_W;
array_y_bot               = Y_bl(1);
array_y_top               = Y_bl(1) + Total_H;
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

% --- 1. Initialization and setup ---
try
    ModelUtil.remove('Model_Verify');
catch
end
model = ModelUtil.create('Model_Verify');
model.modelPath(pwd);
model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 3);
model.component('comp1').mesh.create('mesh1');


% --- 1.5 Nonlinear functions ---
model.func.create('step1', 'Step');
model.func.create('int3', 'Interpolation');
model.func.create('an1', 'Analytic');
model.func.create('an2', 'Analytic');

model.func('step1').set('smooth', 10);
model.func('step1').set('smoothactive', false);
model.func('step1').set('locationdef', 'beginning');

model.func('int3').set('funcname', 'cond');
model.func('int3').set('table', {'345.2295' '126262.6263'; '355.3665' '126869.6581'; '379.9968' '125588.1464'; '414.7557' '122822.7791'; '452.0016' '118411.6809'; '508.4174' '113110.2694'; '573.504' '106721.9817'; '698.3717' '98288.57485'; '711.813' '92898.35679'; '788.35' '86846.02435'; '836.63' '80496.10322'; '905.1918' '75271.95027'});
model.func('int3').set('fununit', {'S/m'});
model.func('int3').set('argunit', {'K'});

model.func('an1').set('funcname', 'Cp_Si');
model.func('an1').set('expr', '1000/28.09*(22.82+3.9*(x/1000)-0.0829*(x/1000)^2+0.0421*(x/1000)^3-0.354*(x/1000)^(-2))');
model.func('an1').set('argunit', {'K'});
model.func('an1').set('plotfixedvalue', {'200'});

model.func('an2').set('funcname', 'k_Si');
model.func('an2').set('expr', 'k_SOI*293.15/T');
model.func('an2').set('args', {'T'});
model.func('an2').set('argunit', {'K'});

% --- 2. Physics setup ---
model.component('comp1').physics.create('ec', 'ConductiveMedia', 'geom1');
model.component('comp1').physics.create('ht', 'HeatTransfer', 'geom1');
model.component('comp1').multiphysics.create('emh1', 'ElectromagneticHeating', 3);
model.component('comp1').multiphysics('emh1').set('EMHeat_physics', 'ec');
model.component('comp1').multiphysics('emh1').set('Heat_physics', 'ht');
model.component('comp1').multiphysics('emh1').selection.all;

% --- 3. Geometry construction aligned with Step1a order ---
fprintf('Building heterogeneous geometry...\n');
geom = model.component('comp1').geom('geom1');
geom.lengthUnit('nm');
geom.create('wp1', 'WorkPlane'); 
geom.feature('wp1').set('unite', true); 
wp1 = geom.feature('wp1').geom;

% Use conductivity parameters supplied by the workflow controller.
model.param.set('k_SOI', sprintf('%f [W/(m*K)]', k_SOI_val));
model.param.set('k_Hdoped', num2str(k_Hdoped_val));
model.param.set('k_Ldoped', num2str(k_Ldoped_val));

model.param.set('L_m', '300000[nm]'); 
centx_val = X_bl(1) + Total_W / 2;
centy_val = Y_bl(1) + Total_H / 2;
model.param.set('centx', sprintf('%d[nm]', centx_val)); 
model.param.set('centy', sprintf('%d[nm]', centy_val));          

% 3.1 Base geometry (WP1). Creation order matches Step1a / original Step3 so
% domain numbering (count = 10) stays unchanged. The sq4 rectangle (formerly
% Diode) is kept geometrically; mat1 now claims it as highdope_Si.
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

num_width_clamped = 0;
for i = 1:Nx
    for j = 1:Ny
        xpos = X_bl(i) + W_arr(i)/2;
        ypos = Y_bl(j) + H_arr(j)/2;
        
        bg_tag = sprintf('bg_%d_%d', i, j);
        wp1.create(bg_tag, 'Rectangle'); 
        wp1.feature(bg_tag).set('size', [W_arr(i) H_arr(j)]); 
        wp1.feature(bg_tag).set('base', 'center'); 
        wp1.feature(bg_tag).set('pos', [xpos ypos]);

        width_rec_raw = width_map(i,j);
        width_rec_max = max(Min_Wire_Width_nm, W_arr(i) - 2 * Min_Side_Gap_nm);
        width_rec = min(max(width_rec_raw, Min_Wire_Width_nm), width_rec_max);
        if abs(width_rec - width_rec_raw) > 1e-9
            num_width_clamped = num_width_clamped + 1;
        end
        r_tag = sprintf('r_%d_%d', i, j);
        wp1.create(r_tag, 'Rectangle'); 
        wp1.feature(r_tag).set('size', [width_rec H_arr(j)]); 
        wp1.feature(r_tag).set('base', 'center'); 
        wp1.feature(r_tag).set('pos', [xpos ypos]);
    end
end
if num_width_clamped > 0
    fprintf('Adjusted %d wire widths to preserve finite no-dope side domains.\n', num_width_clamped);
end
geom.run('wp1');

% 3.2 3D extrusion. Keep Step1a order to preserve domain numbering.
geom.feature.create('ext1', 'Extrude'); geom.feature('ext1').set('workplane', 'wp1'); geom.feature('ext1').selection('input').set({'wp1'}); geom.feature('ext1').setIndex('distance', 220, 0); geom.run('ext1');
geom.create('wp2', 'WorkPlane'); geom.feature('wp2').set('unite', true); wp2 = geom.feature('wp2').geom; wp2.create('sq1', 'Square'); wp2.feature('sq1').set('size', 'L_m'); wp2.feature('sq1').set('base', 'center'); wp2.feature('sq1').set('pos', {'centx' 'centy'}); geom.run('wp2');
geom.feature.create('ext2', 'Extrude'); geom.feature('ext2').set('workplane', 'wp2'); geom.feature('ext2').selection('input').set({'wp2'}); geom.feature('ext2').setIndex('distance', -2000, 0); geom.run('ext2');
geom.create('wp3', 'WorkPlane'); geom.feature('wp3').set('unite', true); geom.feature('wp3').set('quickz', 220); wp3 = geom.feature('wp3').geom; wp3.create('sq1', 'Square'); wp3.feature('sq1').set('size', 'L_m'); wp3.feature('sq1').set('base', 'center'); wp3.feature('sq1').set('pos', {'centx' 'centy'}); geom.run('wp3');
geom.feature.create('ext3', 'Extrude'); geom.feature('ext3').set('workplane', 'wp3'); geom.feature('ext3').selection('input').set({'wp3'}); geom.feature('ext3').setIndex('distance', 10, 0); geom.run('ext3');
geom.create('wp4', 'WorkPlane'); geom.feature('wp4').set('quickz', -2000); wp4 = geom.feature('wp4').geom; wp4.create('sq1', 'Square'); wp4.feature('sq1').set('size', 'L_m'); wp4.feature('sq1').set('base', 'center'); wp4.feature('sq1').set('pos', {'centx' 'centy'}); geom.run('wp4');
geom.feature.create('ext4', 'Extrude'); geom.feature('ext4').set('workplane', 'wp4'); geom.feature('ext4').selection('input').set({'wp4'}); geom.feature('ext4').setIndex('distance', -20000, 0); geom.run('ext4');
geom.run; 

% --- 4. Material domain assignment ---
fprintf('Assigning materials by topology sequence...\n');
high_dope = []; 
no_dope = []; 
count = 10; 
for i = 1:Nx   
    temp1 = count : (count + Ny - 1);
    temp2 = (count + Ny) : (count + 2*Ny - 1);
    temp3 = (count + 2*Ny) : (count + 3*Ny - 1);
    count = count + 3*Ny;
    high_dope = [high_dope temp2]; 
    no_dope = [no_dope temp1 temp3]; 
end

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').label('highdope_Si');
model.component('comp1').material('mat1').propertyGroup('def').set('density', '2329');
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', 'Cp_Si(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'k_Hdoped*k_Si(T)' '0' '0' '0' 'k_Hdoped*k_Si(T)' '0' '0' '0' 'k_Hdoped*k_Si(T)'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'13' '0' '0' '0' '13' '0' '0' '0' '13'});
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'cond(T)' '0' '0' '0' 'cond(T)' '0' '0' '0' 'cond(T)'});
% Domain 6 (former Diode rectangle sq4) is folded into highdope_Si.
model.component('comp1').material('mat1').selection.set([5 6 7 8 9]);

model.component('comp1').material.create('mat4', 'Common');
model.component('comp1').material('mat4').label('lowdope_Si');
model.component('comp1').material('mat4').propertyGroup('def').set('density', '2329');
model.component('comp1').material('mat4').propertyGroup('def').set('heatcapacity', 'Cp_Si(T)');
model.component('comp1').material('mat4').propertyGroup('def').set('thermalconductivity', {'k_Ldoped*k_Si(T)' '0' '0' '0' 'k_Ldoped*k_Si(T)' '0' '0' '0' 'k_Ldoped*k_Si(T)'});
model.component('comp1').material('mat4').propertyGroup('def').set('relpermittivity', {'12' '0' '0' '0' '12' '0' '0' '0' '12'});
model.component('comp1').material('mat4').propertyGroup('def').set('electricconductivity', {'cond(T)/6' '0' '0' '0' 'cond(T)/6' '0' '0' '0' 'cond(T)/6'}); 
model.component('comp1').material('mat4').selection.set(high_dope);       

model.component('comp1').material.create('mat2', 'Common');
model.component('comp1').material('mat2').label('Silicon');
model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', {'1e-12[S/m]' '0' '0' '0' '1e-12[S/m]' '0' '0' '0' '1e-12[S/m]'});
model.component('comp1').material('mat2').propertyGroup('def').set('thermalconductivity', {'k_Si(T)' '0' '0' '0' 'k_Si(T)' '0' '0' '0' 'k_Si(T)'});
model.component('comp1').material('mat2').propertyGroup('def').set('heatcapacity', 'Cp_Si(T)');
model.component('comp1').material('mat2').propertyGroup('def').set('density', '2329[kg/m^3]');
model.component('comp1').material('mat2').propertyGroup('def').set('relpermittivity', {'11.7' '0' '0' '0' '11.7' '0' '0' '0' '11.7'});
model.component('comp1').material('mat2').selection.set([3 no_dope]);   

model.component('comp1').material.create('mat3', 'Common');
model.component('comp1').material('mat3').label('Silica glass');
model.component('comp1').material('mat3').propertyGroup('def').set('electricconductivity', {'1e-14[S/m]' '0' '0' '0' '1e-14[S/m]' '0' '0' '0' '1e-14[S/m]'});
model.component('comp1').material('mat3').propertyGroup('def').set('heatcapacity', '703[J/(kg*K)]');
model.component('comp1').material('mat3').propertyGroup('def').set('relpermittivity', {'3.75' '0' '0' '0' '3.75' '0' '0' '0' '3.75'});
model.component('comp1').material('mat3').propertyGroup('def').set('density', '2203[kg/m^3]');
model.component('comp1').material('mat3').propertyGroup('def').set('thermalconductivity', {'1.38[W/(m*K)]' '0' '0' '0' '1.38[W/(m*K)]' '0' '0' '0' '1.38[W/(m*K)]'});
model.component('comp1').material('mat3').selection.set([2 4]);   

model.component('comp1').material.create('mat5', 'Common');
model.component('comp1').material('mat5').label('Silicon substrate');
model.component('comp1').material('mat5').propertyGroup('def').set('electricconductivity', {'1e-12[S/m]' '0' '0' '0' '1e-12[S/m]' '0' '0' '0' '1e-12[S/m]'});
model.component('comp1').material('mat5').propertyGroup('def').set('thermalconductivity', {'130[W/(m*K)]*293.15[K]/T' '0' '0' '0' '130[W/(m*K)]*293.15[K]/T' '0' '0' '0' '130[W/(m*K)]*293.15[K]/T'});
model.component('comp1').material('mat5').propertyGroup('def').set('heatcapacity', 'Cp_Si(T)');
model.component('comp1').material('mat5').propertyGroup('def').set('density', '2329[kg/m^3]');
model.component('comp1').material('mat5').propertyGroup('def').set('relpermittivity', {'11.7' '0' '0' '0' '11.7' '0' '0' '0' '11.7'});
model.component('comp1').material('mat5').selection.set([1]);

% --- 5. Physics boundary conditions ---
model.component('comp1').physics('ec').create('gnd1', 'Ground', 2);
model.component('comp1').physics('ec').feature('gnd1').selection.set([21]);

model.component('comp1').physics('ec').create('pot1', 'ElectricPotential', 2);
model.component('comp1').physics('ec').feature('pot1').selection.set([35]);
model.component('comp1').physics('ec').feature('pot1').set('V0', sprintf('%f*step1(t/(1[ns]))', Voltage_Use)); 

model.component('comp1').physics('ht').create('temp1', 'TemperatureBoundary', 2);
model.component('comp1').physics('ht').feature('temp1').selection.set([3]);

model.component('comp1').physics('ht').create('hf1', 'HeatFluxBoundary', 2);
model.component('comp1').physics('ht').feature('hf1').selection.set([13]);
model.component('comp1').physics('ht').feature('hf1').set('HeatFluxType', 'ConvectiveHeatFlux');
model.component('comp1').physics('ht').feature('hf1').set('h', 10);

% --- 6. Mesh ---
model.component('comp1').mesh('mesh1').create('ftet1', 'FreeTet');
model.component('comp1').mesh('mesh1').feature('ftet1').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('ftet1').selection.set([high_dope no_dope]);
model.component('comp1').mesh('mesh1').feature('ftet1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('custom', true);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmax', 400); 
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('ftet1').feature('size1').set('hmin', 10);  

model.component('comp1').mesh('mesh1').create('ftet2', 'FreeTet');
model.component('comp1').mesh('mesh1').feature('ftet2').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('ftet2').selection.set([5 6 7 8 9]);
model.component('comp1').mesh('mesh1').feature('ftet2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftet2').feature('size1').set('hauto', 6);

model.component('comp1').mesh('mesh1').create('ftet3', 'FreeTet');
model.component('comp1').mesh('mesh1').feature('ftet3').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('ftet3').selection.set([3]);
model.component('comp1').mesh('mesh1').feature('ftet3').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftet3').feature('size1').set('hauto', 8);

model.component('comp1').mesh('mesh1').create('swe1', 'Sweep');
model.component('comp1').mesh('mesh1').run;

% --- 7. Study setup ---
model.study.create('std1');
model.study('std1').create('time', 'Transient');
model.study('std1').feature('time').set('tlist', Time_List);
model.study('std1').feature('time').setSolveFor('/physics/ec', true);
model.study('std1').feature('time').setSolveFor('/physics/ht', true);
model.study('std1').feature('time').setSolveFor('/multiphysics/emh1', true);

model.sol.create('sol1');
model.sol('sol1').attach('std1');
model.sol('sol1').createAutoSequence('std1');

model.sol('sol1').feature('t1').set('tstepsbdf', 'strict');

% --- 8. Automated result plots and data tables ---
try
    model.result.table.create('evl2', 'Table');
    model.result.table('evl2').label('Evaluation 2D');
    
    model.result.dataset.create('cpl1', 'CutPlane');
    model.result.dataset('cpl1').set('quickplane', 'xy');
    model.result.dataset('cpl1').set('quickz', 110);
    
    model.result.dataset.create('cln1', 'CutLine3D');
    model.result.dataset('cln1').set('genpoints', [centx_val, -35000, 110; centx_val, 75000, 110]);

    model.result.dataset.create('cln2', 'CutLine3D');
    model.result.dataset('cln2').set('genpoints', [-35000, centy_val, 110; 95000, centy_val, 110]);
    
    model.result.create('pg1', 'PlotGroup3D');
    model.result('pg1').label('Electric Potential (ec)');
    model.result('pg1').create('vol1', 'Volume');
    model.result('pg1').feature('vol1').set('colortable', 'Dipole');
    
    model.result.create('pg2', 'PlotGroup3D');
    model.result('pg2').label('Electric Field (ec)');
    model.result('pg2').create('mslc1', 'Multislice');
    model.result('pg2').feature('mslc1').set('expr', 'ec.normE');
    model.result('pg2').feature('mslc1').set('colortable', 'Prism');
    model.result('pg2').create('strmsl1', 'StreamlineMultislice');
    model.result('pg2').feature('strmsl1').set('expr', {'ec.Ex' 'ec.Ey' 'ec.Ez'});
    model.result('pg2').feature('strmsl1').set('color', 'black');
    
    model.result.create('pg3', 'PlotGroup3D');
    model.result('pg3').label('Temperature (ht)');
    model.result('pg3').create('vol1', 'Volume');
    model.result('pg3').feature('vol1').set('expr', 'T');
    model.result('pg3').feature('vol1').set('colortable', 'HeatCameraLight');
    
    model.result.create('pg4', 'PlotGroup1D');
    model.result('pg4').set('xlabel', 'Y-coordinate (nm)');
    model.result('pg4').set('ylabel', 'Temperature (K)');
    model.result('pg4').create('lngr1', 'LineGraph');
    model.result('pg4').feature('lngr1').set('data', 'cln1');
    model.result('pg4').feature('lngr1').set('expr', 'T');
    
    model.result.create('pg5', 'PlotGroup1D');
    model.result('pg5').set('xlabel', 'X-coordinate (nm)');
    model.result('pg5').set('ylabel', 'Temperature (K)');
    model.result('pg5').create('lngr1', 'LineGraph');
    model.result('pg5').feature('lngr1').set('data', 'cln2');
    model.result('pg5').feature('lngr1').set('expr', 'T');
    
    model.result.create('pg6', 'PlotGroup2D');
    model.result('pg6').create('surf1', 'Surface');
    model.result('pg6').feature('surf1').set('data', 'cpl1');
    model.result('pg6').feature('surf1').set('expr', 'T');
    
    fprintf('Automated result plots, tables, and probes configured successfully.\n');
catch
    fprintf('Skipped post-processing plot configuration. The main physics model is unaffected.\n');
end

% Run the COMSOL solve for the verification model.
fprintf('Submitting verification model solve to COMSOL (iter = %d)...\n', iter);
model.study('std1').run;


% --- 9. Save model ---
mphsave(model, OutFile);
fprintf('Heterogeneous model build complete. Saved to: %s\n', OutFile);

end
