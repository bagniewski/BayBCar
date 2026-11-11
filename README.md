# BayBCar

BayBCar is a class of Bayesian proxy system models (PSMs) for benthic δ13C, incorporating the effects of carbonate chemistry, sea ice concentration, and species-specific differences. It enables probabilistic forward and inverse modelling of benthic foraminiferal δ13C and δ13C of DIC, respectively.

Depending on the input variables, one of 8 BayBCar formulations will be selected:
* A3Dp: 3 predictors (δ13C_dic, Depth, Species)
* A3Op: 3 predictors (δ13C_dic, Omega, Species)
* A4DIp: 4 predictors (δ13C_dic, Depth, SIC, Species)
* A4DIp: 4 predictors (δ13C_dic, Omega, SIC, Species)
* W2D: 2 predictors (δ13C_dic, Depth)
* W2O: 2 predictors (δ13C_dic, Omega)
* W3DI: 3 predictors (δ13C_dic, Depth, SIC)
* W3OI: 3 predictors (δ13C_dic, Omega, SIC)  

For a detailed description of BayBCar PSMs please see: Bagniewski, W., Osman, M., Aquino-López, M., Rathburn, A., & Muschitiello, F. (2026). BayBCar: A Bayesian calibration model for δ13C in benthic foraminifera. Paleoceanography and Paleoclimatology.

For a Python version of BayBCar, see [PyBayBCar](https://github.com/bagniewski/PyBayBCar).

## External dependencies

* [CO2SYS](https://zenodo.org/records/7552554) - calculates Omega from available variables if not provided by the user. BayBCar will work without it by defaulting to a formulation that uses depth in place of Omega.

## Files

* **BayBCar.m**: a Matlab function that selects the appropriate BayBCar calibration based on input data provided by the user and performs probabilistic forward and inverse modelling of benthic δ13C. 

* **BayBCar_calibrations.mat**: a data file used by **BayBCar.m**, containing the calibrated BayBCar parameters.

## Usage

### Syntax
* Forward modelling:
```
[d13C_cib,upper,lower] = BayBCar(d13C_dic,Depth,Omega,SIC,Species,CI,DIC,TA,pH,pCO2,T,S,Si,PO4,"forward")
```
* Inverse modelling:
```
[d13C_dic,upper,lower] = BayBCar(d13C_cib,Depth,Omega,SIC,Species,CI,DIC,TA,pH,pCO2,T,S,Si,PO4,"inverse")
```
Use `[]` to skip arguments, for example:
`[d13C_cib,upper,lower] = BayBCar(d13C_dic,Depth,[],SIC,[],0.75)`

Prediction direction (`"forward"` or `"inverse"`) can be specified as the last input argument. When not specified, `"forward"` is used. Example:
`[d13C_dic,upper,lower] = BayBCar(d13C_cib,[],Omega,"inverse")`

### Required Inputs
- `d13C_x`: First input argument, represents either d13C_dic (in forward modelling) or d13C_cib (in inverse modelling):
  - `d13C_dic`: δ13C of dissolved inorganic carbon (permil)
  - `d13C_cib`: δ13C of Cibicidoides foraminifera (permil, C. wuellerstorfi by default)
- At least one of these:
  - `Depth`: Water depth (meters)
  - `Omega`: Calcite saturation state

### Optional Inputs
- `Depth`: Water depth (meters)  
- `Omega`: Calcite saturation state  
- `SIC`: Summer sea ice concentration (decimal)  
- `Species`: Foraminiferal species. 1 = C. spp; 2 (default) = C. wuellerstorfi; 3 = C. kullenbergi; 4 = C. lobatulus; 5 = C. pachyderma; 6 = C. floridanus; 7 = C. mckannai
- `CI`: Confidence interval in decimal format. 0.95 (default) = 95% CI.
- `direction`: Direction of prediction. "forward" (default) = forward modelling to predict d13C_cib; "inverse" = inverse modelling to predict d13C_dic.

### Optional CO2SYS inputs for calculating Omega (at least two of these are required: DIC, TA, pH, pCO2)
- `DIC`: Dissolved inorganic carbon (µmol/kg)  
- `TA`: Total alkalinity (µmol/kg)  
- `pH`: pH  
- `pCO2`: pCO2 (µatm)  
- `T`: Temperature (Celsius)  
- `S`: Salinity  
- `Si`: Silicate (µmol/kg)  
- `PO4`: Phosphate (µmol/kg)  

### Output
- `d13C_pred`: represents either d13C_cib (in forward modelling) or d13C_dic (in inverse modelling)
  - `d13C_cib`: δ13C of Cibicidoides foraminifera (permil, C. wuellerstorfi by default)
  - `d13C_dic`: δ13C of dissolved inorganic carbon (permil)  
- `upper`: upper bound of the confidence interval (95% by default)  
- `lower`: lower bound of the confidence interval (95% by default)

## Examples

Sample data:
```
time     = 1:10;
d13C_dic = [0.43, 0.43, 0.42, 0.40, 0.41, 0.33, 0.33, 0.30, 0.22, 0.31];
Omega    = [0.85, 0.83, 0.81, 0.82, 0.80, 0.70, 0.71, 0.70, 0.68, 0.66];
SIC      = 0.1;
Depth    = 3420;
DIC      = [2255, 2250, 2255, 2251, 2264, 2265, 2270, 2276, 2280, 2278];
TA       = [2365, 2361, 2366, 2361, 2374, 2377, 2380, 2386, 2393, 2397];
% Sample d13C_cib data for inverse modelling:
d13C_cib = [0.33, 0.33, 0.32, 0.30, 0.31, 0.23, 0.23, 0.20, 0.12, 0.21];
```

### Forward modelling

Example 1: W2D model
```
[d13C_cib1,upper1,lower1] = BayBCar(d13C_dic,Depth);
```
Example 2: A4OIp model. Using `BayBCar(d13C_dic,[],Omega,SIC,1)` gives identical result.
```
[d13C_cib2,upper2,lower2] = BayBCar(d13C_dic,Depth,Omega,SIC,1);
```
Example 3: A4OIp model using CO2SYS
```
[d13C_cib3,upper3,lower3] = BayBCar(d13C_dic,Depth,[],SIC,1,[],DIC,TA);
```
Plot the forward modelled d13C_cib from the 3 examples above:
```
figure; 
plot(time,d13C_cib1,'-*r','LineWidth',2); hold on;
fill([time'; flip(time')], [upper1; flip(lower1)], 'r', 'EdgeColor','none','FaceAlpha',0.15);
plot(time,d13C_cib2,'-*b','LineWidth',2);
fill([time'; flip(time')], [upper2; flip(lower2)], 'b', 'EdgeColor','none','FaceAlpha',0.15);
plot(time,d13C_cib3,'-*k','LineWidth',2);
fill([time'; flip(time')], [upper3; flip(lower3)], 'k', 'EdgeColor','none','FaceAlpha',0.15);
xlabel('time'); ylabel('Forward modelled d13C\_cib')
legend('W2D','','A4OIp','','A4OIp (CO2SYS)','')
```

### Inverse modelling

Example 4: A4OIp model
```
[d13C_dic4,upper4,lower4] = BayBCar(d13C_cib,Depth,Omega,SIC,1,"inverse");
```
Plot the inverse modelled d13C_dic:
```
figure; 
plot(time,d13C_dic4,'-*b','LineWidth',2); hold on;
fill([time'; flip(time')], [upper4; flip(lower4)], 'b', 'EdgeColor','none','FaceAlpha',0.15);
xlabel('time'); ylabel('Inverse modelled d13C\_dic')
legend('A4OIp','')
```
