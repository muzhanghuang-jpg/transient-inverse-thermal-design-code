function Build_Plain_Heater(Voltage_Use, mesh_only)
if nargin < 2 || isempty(mesh_only), mesh_only = false; end
%BUILD_PLAIN_HEATER Plain uniformly-doped 50x50 um heater benchmark.
%
% Same pad / lead / substrate / silica stack and thermal/electrical BCs as the
% inverse-design Step3, but the patterned pixel array is a single solid low-dope
% heater block (domain 10). Domain/boundary IDs were verified by
% Build_Plain_Heater_Geom (DO NOT reuse the patterned Step3 IDs blindly):
%   materials  : substrate=1, silica=[2 4], field(undoped)=3, pads=[5 6 7 8 9], HEATER=10
%   boundaries : gnd1=21, pot1=35, temp1=3, hf1=13
%
% Applies Voltage_Use (default 21 V, sized for total power = inverse-design 1.224 W),
% solves the 13 us transient, then prints total dissipated power and the heater
% temperature uniformity so the voltage can be fine-tuned to match the power.

if nargin < 1 || isempty(Voltage_Use), Voltage_Use = 21; end

mliPath = comsol_mli_path();
if exist(mliPath, 'dir') == 7, addpath(mliPath); end
try, mphstart(comsol_server_port()); catch, end
import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.showProgress(false);

k_SOI = 130; k_Hdoped = 1; k_Ldoped = 1;
targetSettings = get_thermal_target_settings();
Target_Temp = targetSettings.Base_Temp_K;
Target_Time = targetSettings.Target_Time_s;
Time_Step   = targetSettings.Time_Step_s;

% --- geometry params (match inverse-design stack) ---
Total_W = 50000; Total_H = 50000;
array_x_left = -4000; array_x_right = array_x_left + Total_W;
array_y_bot  = -4000; array_y_top  = array_y_bot  + Total_H;
centx_val = array_x_left + Total_W/2; centy_val = array_y_bot + Total_H/2;
lead_x_overhang = 10000; lead_height = 10000;
pad_height_top = 40000; pad_height_bot_main = 35000; pad_height_bot_transition = 5000;
pad_width = Total_W + 2*lead_x_overhang;
sq1_pos_y = array_y_top + lead_height + pad_height_top/2;
sq2_pos_y = array_y_bot - lead_height - pad_height_bot_transition - pad_height_bot_main/2;
sq4_pos_y = array_y_bot - lead_height - pad_height_bot_transition/2;
pol2_table = [array_x_left, array_y_top; array_x_right, array_y_top; array_x_right+lead_x_overhang, array_y_top+lead_height; array_x_left-lead_x_overhang, array_y_top+lead_height];
pol1_table = [array_x_left, array_y_bot; array_x_left-lead_x_overhang, array_y_bot-lead_height; array_x_right+lead_x_overhang, array_y_bot-lead_height; array_x_right, array_y_bot];

try, ModelUtil.remove('Plain_Heater'); catch, end
model = ModelUtil.create('Plain_Heater');
model.modelPath(pwd);
model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 3);
geom = model.component('comp1').geom('geom1');
geom.lengthUnit('nm');

model.param.set('k_SOI', sprintf('%f [W/(m*K)]', k_SOI));
model.param.set('k_Hdoped', num2str(k_Hdoped));
model.param.set('k_Ldoped', num2str(k_Ldoped));
model.param.set('L_m', '300000[nm]');
model.param.set('centx', sprintf('%d[nm]', centx_val));
model.param.set('centy', sprintf('%d[nm]', centy_val));

% --- nonlinear material functions (same as Step3) ---
model.func.create('step1', 'Step'); model.func('step1').set('smooth', 10); model.func('step1').set('smoothactive', false); model.func('step1').set('locationdef', 'beginning');
model.func.create('int3', 'Interpolation'); model.func('int3').set('funcname', 'cond');
model.func('int3').set('table', {'345.2295' '126262.6263'; '355.3665' '126869.6581'; '379.9968' '125588.1464'; '414.7557' '122822.7791'; '452.0016' '118411.6809'; '508.4174' '113110.2694'; '573.504' '106721.9817'; '698.3717' '98288.57485'; '711.813' '92898.35679'; '788.35' '86846.02435'; '836.63' '80496.10322'; '905.1918' '75271.95027'});
model.func('int3').set('fununit', {'S/m'}); model.func('int3').set('argunit', {'K'});
model.func.create('an1', 'Analytic'); model.func('an1').set('funcname', 'Cp_Si');
model.func('an1').set('expr', '1000/28.09*(22.82+3.9*(x/1000)-0.0829*(x/1000)^2+0.0421*(x/1000)^3-0.354*(x/1000)^(-2))');
model.func('an1').set('argunit', {'K'}); model.func('an1').set('plotfixedvalue', {'200'});
model.func.create('an2', 'Analytic'); model.func('an2').set('funcname', 'k_Si'); model.func('an2').set('expr', 'k_SOI*293.15/T'); model.func('an2').set('args', {'T'}); model.func('an2').set('argunit', {'K'});

% --- physics (create before geometry like Step3) ---
model.component('comp1').physics.create('ec', 'ConductiveMedia', 'geom1');
model.component('comp1').physics.create('ht', 'HeatTransfer', 'geom1');
model.component('comp1').multiphysics.create('emh1', 'ElectromagneticHeating', 3);
model.component('comp1').multiphysics('emh1').set('EMHeat_physics', 'ec');
model.component('comp1').multiphysics('emh1').set('Heat_physics', 'ht');
model.component('comp1').multiphysics('emh1').selection.all;

% --- geometry ---
fprintf('Building plain-heater geometry...\n');
geom.create('wp1', 'WorkPlane'); geom.feature('wp1').set('unite', true); wp1 = geom.feature('wp1').geom;
wp1.create('sq1', 'Rectangle'); wp1.feature('sq1').set('size', [pad_width pad_height_top]); wp1.feature('sq1').set('base', 'center'); wp1.feature('sq1').set('pos', {'centx' num2str(sq1_pos_y)});
wp1.create('pol2', 'Polygon'); wp1.feature('pol2').set('source', 'table'); wp1.feature('pol2').set('table', pol2_table);
wp1.create('sq2', 'Rectangle'); wp1.feature('sq2').set('size', [pad_width pad_height_bot_main]); wp1.feature('sq2').set('base', 'center'); wp1.feature('sq2').set('pos', {'centx' num2str(sq2_pos_y)});
wp1.create('sq4', 'Rectangle'); wp1.feature('sq4').set('size', [pad_width pad_height_bot_transition]); wp1.feature('sq4').set('base', 'center'); wp1.feature('sq4').set('pos', {'centx' num2str(sq4_pos_y)});
wp1.create('pol1', 'Polygon'); wp1.feature('pol1').set('source', 'table'); wp1.feature('pol1').set('table', pol1_table);
wp1.create('sq3', 'Square'); wp1.feature('sq3').set('size', 'L_m'); wp1.feature('sq3').set('base', 'center'); wp1.feature('sq3').set('pos', {'centx' 'centy'});
wp1.create('heater', 'Rectangle'); wp1.feature('heater').set('size', [Total_W Total_H]); wp1.feature('heater').set('base', 'center'); wp1.feature('heater').set('pos', {'centx' 'centy'});
geom.run('wp1');
geom.feature.create('ext1', 'Extrude'); geom.feature('ext1').set('workplane', 'wp1'); geom.feature('ext1').selection('input').set({'wp1'}); geom.feature('ext1').setIndex('distance', 220, 0); geom.run('ext1');
geom.create('wp2', 'WorkPlane'); geom.feature('wp2').set('unite', true); wp2 = geom.feature('wp2').geom; wp2.create('sq1', 'Square'); wp2.feature('sq1').set('size', 'L_m'); wp2.feature('sq1').set('base', 'center'); wp2.feature('sq1').set('pos', {'centx' 'centy'}); geom.run('wp2');
geom.feature.create('ext2', 'Extrude'); geom.feature('ext2').set('workplane', 'wp2'); geom.feature('ext2').selection('input').set({'wp2'}); geom.feature('ext2').setIndex('distance', -2000, 0); geom.run('ext2');
geom.create('wp3', 'WorkPlane'); geom.feature('wp3').set('unite', true); geom.feature('wp3').set('quickz', 220); wp3 = geom.feature('wp3').geom; wp3.create('sq1', 'Square'); wp3.feature('sq1').set('size', 'L_m'); wp3.feature('sq1').set('base', 'center'); wp3.feature('sq1').set('pos', {'centx' 'centy'}); geom.run('wp3');
geom.feature.create('ext3', 'Extrude'); geom.feature('ext3').set('workplane', 'wp3'); geom.feature('ext3').selection('input').set({'wp3'}); geom.feature('ext3').setIndex('distance', 10, 0); geom.run('ext3');
geom.create('wp4', 'WorkPlane'); geom.feature('wp4').set('quickz', -2000); wp4 = geom.feature('wp4').geom; wp4.create('sq1', 'Square'); wp4.feature('sq1').set('size', 'L_m'); wp4.feature('sq1').set('base', 'center'); wp4.feature('sq1').set('pos', {'centx' 'centy'}); geom.run('wp4');
geom.feature.create('ext4', 'Extrude'); geom.feature('ext4').set('workplane', 'wp4'); geom.feature('ext4').selection('input').set({'wp4'}); geom.feature('ext4').setIndex('distance', -20000, 0); geom.run('ext4');
geom.run;

% --- materials (VERIFIED IDs) ---
fprintf('Assigning materials (heater=dom10, pads=[5 6 7 8 9], field=3, silica=[2 4], substrate=1)...\n');
model.component('comp1').material.create('mat1', 'Common'); model.component('comp1').material('mat1').label('highdope_Si');
model.component('comp1').material('mat1').propertyGroup('def').set('density', '2329');
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', 'Cp_Si(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'k_Hdoped*k_Si(T)' '0' '0' '0' 'k_Hdoped*k_Si(T)' '0' '0' '0' 'k_Hdoped*k_Si(T)'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'13' '0' '0' '0' '13' '0' '0' '0' '13'});
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'cond(T)' '0' '0' '0' 'cond(T)' '0' '0' '0' 'cond(T)'});
model.component('comp1').material('mat1').selection.set([5 6 7 8 9]);

model.component('comp1').material.create('mat4', 'Common'); model.component('comp1').material('mat4').label('lowdope_Si');
model.component('comp1').material('mat4').propertyGroup('def').set('density', '2329');
model.component('comp1').material('mat4').propertyGroup('def').set('heatcapacity', 'Cp_Si(T)');
model.component('comp1').material('mat4').propertyGroup('def').set('thermalconductivity', {'k_Ldoped*k_Si(T)' '0' '0' '0' 'k_Ldoped*k_Si(T)' '0' '0' '0' 'k_Ldoped*k_Si(T)'});
model.component('comp1').material('mat4').propertyGroup('def').set('relpermittivity', {'12' '0' '0' '0' '12' '0' '0' '0' '12'});
model.component('comp1').material('mat4').propertyGroup('def').set('electricconductivity', {'cond(T)/6' '0' '0' '0' 'cond(T)/6' '0' '0' '0' 'cond(T)/6'});
model.component('comp1').material('mat4').selection.set([10]);   % THE HEATER

model.component('comp1').material.create('mat2', 'Common'); model.component('comp1').material('mat2').label('Silicon');
model.component('comp1').material('mat2').propertyGroup('def').set('electricconductivity', {'1e-12[S/m]' '0' '0' '0' '1e-12[S/m]' '0' '0' '0' '1e-12[S/m]'});
model.component('comp1').material('mat2').propertyGroup('def').set('thermalconductivity', {'k_Si(T)' '0' '0' '0' 'k_Si(T)' '0' '0' '0' 'k_Si(T)'});
model.component('comp1').material('mat2').propertyGroup('def').set('heatcapacity', 'Cp_Si(T)');
model.component('comp1').material('mat2').propertyGroup('def').set('density', '2329[kg/m^3]');
model.component('comp1').material('mat2').propertyGroup('def').set('relpermittivity', {'11.7' '0' '0' '0' '11.7' '0' '0' '0' '11.7'});
model.component('comp1').material('mat2').selection.set([3]);

model.component('comp1').material.create('mat3', 'Common'); model.component('comp1').material('mat3').label('Silica glass');
model.component('comp1').material('mat3').propertyGroup('def').set('electricconductivity', {'1e-14[S/m]' '0' '0' '0' '1e-14[S/m]' '0' '0' '0' '1e-14[S/m]'});
model.component('comp1').material('mat3').propertyGroup('def').set('heatcapacity', '703[J/(kg*K)]');
model.component('comp1').material('mat3').propertyGroup('def').set('relpermittivity', {'3.75' '0' '0' '0' '3.75' '0' '0' '0' '3.75'});
model.component('comp1').material('mat3').propertyGroup('def').set('density', '2203[kg/m^3]');
model.component('comp1').material('mat3').propertyGroup('def').set('thermalconductivity', {'1.38[W/(m*K)]' '0' '0' '0' '1.38[W/(m*K)]' '0' '0' '0' '1.38[W/(m*K)]'});
model.component('comp1').material('mat3').selection.set([2 4]);

model.component('comp1').material.create('mat5', 'Common'); model.component('comp1').material('mat5').label('Silicon substrate');
model.component('comp1').material('mat5').propertyGroup('def').set('electricconductivity', {'1e-12[S/m]' '0' '0' '0' '1e-12[S/m]' '0' '0' '0' '1e-12[S/m]'});
model.component('comp1').material('mat5').propertyGroup('def').set('thermalconductivity', {'130[W/(m*K)]*293.15[K]/T' '0' '0' '0' '130[W/(m*K)]*293.15[K]/T' '0' '0' '0' '130[W/(m*K)]*293.15[K]/T'});
model.component('comp1').material('mat5').propertyGroup('def').set('heatcapacity', 'Cp_Si(T)');
model.component('comp1').material('mat5').propertyGroup('def').set('density', '2329[kg/m^3]');
model.component('comp1').material('mat5').propertyGroup('def').set('relpermittivity', {'11.7' '0' '0' '0' '11.7' '0' '0' '0' '11.7'});
model.component('comp1').material('mat5').selection.set([1]);

% --- boundary conditions (VERIFIED IDs) ---
model.component('comp1').physics('ec').create('gnd1', 'Ground', 2); model.component('comp1').physics('ec').feature('gnd1').selection.set([21]);
model.component('comp1').physics('ec').create('pot1', 'ElectricPotential', 2); model.component('comp1').physics('ec').feature('pot1').selection.set([35]);
model.component('comp1').physics('ec').feature('pot1').set('V0', sprintf('%f*step1(t/(1[ns]))', Voltage_Use));
model.component('comp1').physics('ht').create('temp1', 'TemperatureBoundary', 2); model.component('comp1').physics('ht').feature('temp1').selection.set([3]); model.component('comp1').physics('ht').feature('temp1').set('T0', '293.15[K]');
model.component('comp1').physics('ht').create('hf1', 'HeatFluxBoundary', 2); model.component('comp1').physics('ht').feature('hf1').selection.set([13]); model.component('comp1').physics('ht').feature('hf1').set('HeatFluxType', 'ConvectiveHeatFlux'); model.component('comp1').physics('ht').feature('hf1').set('h', 10);

% --- mesh: global SWEPT mesh (clean z-layered stack; free-tet would explode on
%     the 300 um thin layers). One element-size on the auto source face. ---
fprintf('Meshing (swept)...\n');
model.component('comp1').mesh.create('mesh1');
% Graded mesh: fine on heater (2 um), gentle growth (hgrad 1.3) out to a 15 um
% far-field cap -- smooth transition for publication-quality field plots.
model.component('comp1').mesh('mesh1').feature('size').set('custom', true);
model.component('comp1').mesh('mesh1').feature('size').set('hmax', 15000);
model.component('comp1').mesh('mesh1').feature('size').set('hmin', 1000);
model.component('comp1').mesh('mesh1').feature('size').set('hgrad', 1.3);
model.component('comp1').mesh('mesh1').create('size2', 'Size');           % refine the heater
model.component('comp1').mesh('mesh1').feature('size2').selection.geom('geom1', 3);
model.component('comp1').mesh('mesh1').feature('size2').selection.set([10]);
model.component('comp1').mesh('mesh1').feature('size2').set('custom', true);
model.component('comp1').mesh('mesh1').feature('size2').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('size2').set('hmax', 2000);
model.component('comp1').mesh('mesh1').create('swe1', 'Sweep');
model.component('comp1').mesh('mesh1').run;
try, st = mphmeshstats(model); fprintf('Mesh elements: %d\n', sum(st.numelem)); catch ME, fprintf('meshstats: %s\n', ME.message); end
if mesh_only
    fprintf('mesh_only mode: stopping before solve.\n');
    try, mphsave(model, 'Plain_Heater_MeshTest.mph'); catch, end
    return;
end

% --- study ---
fprintf('Solving 13 us transient at V = %.4f V...\n', Voltage_Use);
model.study.create('std1'); model.study('std1').create('time', 'Transient');
model.study('std1').feature('time').set('tlist', sprintf('range(0, %.16g, %.16g)', Time_Step, Target_Time));
model.study('std1').feature('time').setSolveFor('/physics/ec', true);
model.study('std1').feature('time').setSolveFor('/physics/ht', true);
model.study('std1').feature('time').setSolveFor('/multiphysics/emh1', true);
model.sol.create('sol1'); model.sol('sol1').attach('std1'); model.sol('sol1').createAutoSequence('std1');
model.sol('sol1').feature('t1').set('tstepsbdf', 'strict');

% --- result template (same datasets / plot groups as the Step3 verification model) ---
try
    model.result.table.create('evl2', 'Table'); model.result.table('evl2').label('Evaluation 2D');
    model.result.dataset.create('cpl1', 'CutPlane'); model.result.dataset('cpl1').set('quickplane', 'xy'); model.result.dataset('cpl1').set('quickz', 110);
    model.result.dataset.create('cln1', 'CutLine3D'); model.result.dataset('cln1').set('genpoints', [centx_val, -35000, 110; centx_val, 75000, 110]);
    model.result.dataset.create('cln2', 'CutLine3D'); model.result.dataset('cln2').set('genpoints', [-35000, centy_val, 110; 95000, centy_val, 110]);
    model.result.create('pg1', 'PlotGroup3D'); model.result('pg1').label('Electric Potential (ec)'); model.result('pg1').create('vol1', 'Volume'); model.result('pg1').feature('vol1').set('colortable', 'Dipole');
    model.result.create('pg2', 'PlotGroup3D'); model.result('pg2').label('Electric Field (ec)'); model.result('pg2').create('mslc1', 'Multislice'); model.result('pg2').feature('mslc1').set('expr', 'ec.normE'); model.result('pg2').feature('mslc1').set('colortable', 'Prism'); model.result('pg2').create('strmsl1', 'StreamlineMultislice'); model.result('pg2').feature('strmsl1').set('expr', {'ec.Ex' 'ec.Ey' 'ec.Ez'}); model.result('pg2').feature('strmsl1').set('color', 'black');
    model.result.create('pg3', 'PlotGroup3D'); model.result('pg3').label('Temperature (ht)'); model.result('pg3').create('vol1', 'Volume'); model.result('pg3').feature('vol1').set('expr', 'T'); model.result('pg3').feature('vol1').set('colortable', 'HeatCameraLight');
    model.result.create('pg4', 'PlotGroup1D'); model.result('pg4').set('xlabel', 'Y-coordinate (nm)'); model.result('pg4').set('ylabel', 'Temperature (K)'); model.result('pg4').create('lngr1', 'LineGraph'); model.result('pg4').feature('lngr1').set('data', 'cln1'); model.result('pg4').feature('lngr1').set('expr', 'T');
    model.result.create('pg5', 'PlotGroup1D'); model.result('pg5').set('xlabel', 'X-coordinate (nm)'); model.result('pg5').set('ylabel', 'Temperature (K)'); model.result('pg5').create('lngr1', 'LineGraph'); model.result('pg5').feature('lngr1').set('data', 'cln2'); model.result('pg5').feature('lngr1').set('expr', 'T');
    model.result.create('pg6', 'PlotGroup2D'); model.result('pg6').create('surf1', 'Surface'); model.result('pg6').feature('surf1').set('data', 'cpl1'); model.result('pg6').feature('surf1').set('expr', 'T');
    fprintf('Result template (cpl1/cln1/cln2 + pg1..pg6) configured.\n');
catch ME
    fprintf('Result template skipped: %s\n', ME.message);
end

model.study('std1').run;

% --- post: total dissipated power + heater temperature uniformity ---
Ptot = NaN;
try, Ptot = mphint2(model, 'ec.Qh', 'volume', 'dataset', 'dset1', 't', Target_Time); catch ME, fprintf('power int (ec.Qh) failed: %s\n', ME.message); end
if isnan(Ptot)
    try, Ptot = mphint2(model, 'ec.Qrh', 'volume', 'dataset', 'dset1', 't', Target_Time); catch, end
end

npx = 42;
xs = linspace(array_x_left + 4000, array_x_right - 4000, npx);   % central 42 um ROI
ys = linspace(array_y_bot  + 4000, array_y_top  - 4000, npx);
[Xg, Yg] = meshgrid(xs, ys);
coords = [Xg(:)'; Yg(:)'; ones(1, numel(Xg)) * 110];
Tg = mphinterp(model, 'T', 'coord', coords, 't', Target_Time);
Tg = Tg(:);
err = Tg - Target_Temp;

fprintf('\n========== PLAIN HEATER RESULT (V = %.4f V, t = %.2f us) ==========\n', Voltage_Use, Target_Time*1e6);
fprintf('Total dissipated power = %.4f W   (inverse-design target = 1.2285 W)\n', Ptot);
if ~isnan(Ptot) && Ptot > 0
    fprintf('  -> voltage to match 1.2285 W: V_new = %.4f V\n', Voltage_Use*sqrt(1.2285/Ptot));
end
fprintf('Central 42x42 um ROI temperature (vs %.0f K target):\n', Target_Temp);
fprintf('  mean = %.2f K,  min = %.2f K,  max = %.2f K\n', mean(Tg), min(Tg), max(Tg));
fprintf('  RMSE = %.2f K,  max+ = %+.2f K,  max- = %+.2f K,  max-abs = %.2f K\n', ...
    sqrt(mean(err.^2)), max(err), min(err), max(abs(err)));

OutFile = sprintf('Plain_Heater_%.2fV.mph', Voltage_Use);
mphsave(model, OutFile);
fprintf('Saved %s\n', OutFile);
end
