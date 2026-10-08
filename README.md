# Transient inverse thermal design: simulation code

MATLAB/COMSOL workflows for inverse design and nonlinear electrothermal
simulation of resistive silicon microheaters. The repository includes three
inverse-designed cases and a uniformly doped reference heater.

## Requirements

- MATLAB R2025b and Optimization Toolbox (`lsqlin`).
- COMSOL Multiphysics 6.4, LiveLink for MATLAB, and licenses for the thermal and
  electrical physics used by the model-building scripts.
- A running COMSOL Multiphysics Server. Default port: 2036.
- Sufficient memory and disk space; solved models can be several gigabytes.

MATLAB and COMSOL are not bundled with this repository.

## Cases

| Case ID | Physical case | Entry point |
|---|---|---|
| `uniform_1um_inverse` | 50 x 50 um heater, 42 x 42 um target, four groups | `run_case('uniform_1um_inverse')` |
| `uniform_1um_plain` | Uniformly doped reference, 20.92 V | `run_case('uniform_1um_plain')` |
| `rectangular_2um_1d1r` | Rectangular 1D1R heater, six groups | `run_case('rectangular_2um_1d1r')` |
| `linear_profile` | Rectangular linear-gradient heater, 66 groups | `run_case('linear_profile')` |

All cases evaluate the temperature field at 13 us. The inverse-designed cases
include response-matrix construction, constrained power optimization,
power-to-linewidth conversion, and nonlinear electrothermal correction.

## Run

Start the COMSOL server separately, using the executable from your installation.
In MATLAB, change to this repository folder:

```matlab
setenv('COMSOL_MLI_PATH', 'C:\Program Files\COMSOL\COMSOL64\Multiphysics\mli');
setenv('COMSOL_SERVER_PORT', '2036');
addpath(pwd);
run_case('uniform_1um_inverse');
```

`COMSOL_MLI_PATH` may be omitted when LiveLink is already on the MATLAB path.
Each entry point writes to a new timestamped result directory. The inverse-design
workflows construct their response matrices and initial layouts before running
nonlinear correction; saved results are not required to start a new calculation.

```matlab
addpath('tools');
smoke_test; % Static parsing and configuration tests; no FEM simulation
```

## Outputs

Run directories contain model files, optimization results, and workflow logs.
`tools/export_temperature_field.m` exports the temperature field and two
cutlines from a solved model using specified coordinates, without modifying or
rerunning the model.

## Attribution and License

Developed by Muzhang Huang with collaborators. The patterned-heater code and
material models draw on earlier work by Khoi Phuong Dao.

Released under the **BSD 3-Clause License**; see `LICENSE`. This license covers
the repository code, not MATLAB or COMSOL.

Please cite this software and the associated paper when available.
Software citation information is provided in `CITATION.cff`.
