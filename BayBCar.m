function [d13C_pred,upper,lower] = BayBCar(d13C_x,varargin)
%
%% DESCRIPTION
% 
% This function uses BayBCar PSM calibrations to perform probabilistic 
% forward and inverse modelling of benthic δ13C.
%
% Depending on the input variables, one of 8 BayBCar formulations will be
% selected:
%
%    A3Dp  : 3 predictors (d13C_dic, Depth, Species)
%    A3Op  : 3 predictors (d13C_dic, Omega, Species)
%    A4DIp : 4 predictors (d13C_dic, Depth, SIC, Species)
%    A4DIp : 4 predictors (d13C_dic, Omega, SIC, Species)
%    W2D   : 2 predictors (d13C_dic, Depth)
%    W2O   : 2 predictors (d13C_dic, Omega)
%    W3DI  : 3 predictors (d13C_dic, Depth, SIC)
%    W3OI  : 3 predictors (d13C_dic, Omega, SIC)
%
% The function also returns the upper and lower bounds of the confidence
% interval specified by the user.
%
% For a detailed description of BayBCar PSMs please see:
%    Bagniewski, W., Osman, M., Aquino-López, M., Rathburn, A., &
%    Muschitiello, F. (2026). BayBCar: A Bayesian calibration model for
%    δ13C in benthic foraminifera. Paleoceanography and Paleoclimatology.
%
% Author: Witold Bagniewski
% Date: 11/11/2026
%
%% DEPENDENCIES
%
% CO2SYS:
%   When Omega is not provided by the user, the function will attempt to
%   calculate it from other available variables using the CO2SYS toolbox:
%   Sharp, J. D., Pierrot, D., Humphreys, M. P., Epitalon, J.-M., Orr, J. 
%   C., Lewis, E. R., & Wallace, D. W. R. (2023). CO2SYSv3 for MATLAB 
%   (v3.2.1). Zenodo. https://doi.org/10.5281/zenodo.7552554.
%
%% USAGE
%
% Forward modelling syntax:
%   [d13C_cib,upper,lower] = BayBCar(d13C_dic,Depth,Omega,SIC,Species,CI,DIC,TA,pH,pCO2,T,S,Si,PO4,"forward")
%
% Inverse modelling syntax:
%   [d13C_dic,upper,lower] = BayBCar(d13C_cib,Depth,Omega,SIC,Species,CI,DIC,TA,pH,pCO2,T,S,Si,PO4,"inverse")
%
% Use [] to skip arguments, for example:
%   [d13C_cib,upper,lower] = BayBCar(d13C_dic,Depth,[],SIC,[],0.75)
%
% Prediction direction ("forward" or "inverse") can be specified as the last input argument. When not specified, "forward" is used. Example:
%   [d13C_dic,upper,lower] = BayBCar(d13C_cib,[],Omega,"inverse")
%
% Keep BayBCar_calibrations.mat in the same directory as this file
%
% Required Inputs:
%   d13C_x   : First input argument, represents either d13C_dic (in forward modelling) or d13C_cib (in inverse modelling):
%    - d13C_dic : d13C of dissolved inorganic carbon (permil)
%    - d13C_cib : d13C of Cibicidoides foraminifera (permil, C. wuellerstorfi by default)
%   At least one of these:
%    - Depth : Water depth (meters)
%    - Omega : Calcite saturation state
%
% Optional Inputs:
%   Depth     : Water depth (meters)
%   Omega     : Calcite saturation state
%   SIC       : Summer sea ice concentration (decimal)
%   Species   : Foraminiferal species. 1 = C. spp; 2 (default) = C. wuellerstorfi; 3 = C. kullenbergi; 4 = C. lobatulus; 5 = C. pachyderma; 6 = C. floridanus; 7 = C. mckannai
%   CI        : Confidence interval in decimal format. 0.95 (default) = 95% CI.
%   direction : Direction of prediction. "forward" (default) = forward modelling to predict d13C_cib; "inverse" = inverse modelling to predict d13C_dic.
%
% Optional CO2SYS inputs for calculating Omega (at least two of these are required: DIC, TA, pH, pCO2):
%   DIC      : Dissolved inorganic carbon (µmol/kg)
%   TA       : Total alkalinity (µmol/kg)
%   pH       : pH
%   pCO2     : pCO2 (µatm)
%   T        : Temperature (Celsius)
%   S        : Salinity
%   Si       : Silicate (µmol/kg)
%   PO4      : Phosphate (µmol/kg)
%
% Output:
%   d13C_pred : represents either d13C_cib (in forward modelling) or d13C_dic (in inverse modelling)
%    - d13C_cib : d13C of Cibicidoides foraminifera (permil, C. wuellerstorfi by default)
%    - d13C_dic : d13C of dissolved inorganic carbon (permil)
%   upper    : upper bound of the confidence interval (95% by default)
%   lower    : lower bound of the confidence interval (95% by default)
%
%% EXAMPLES
%
% % Sample data:
%   time     = 1:10;
%   d13C_dic = [0.43, 0.43, 0.42, 0.40, 0.41, 0.33, 0.33, 0.30, 0.22, 0.31];
%   Omega    = [0.85, 0.83, 0.81, 0.82, 0.80, 0.70, 0.71, 0.70, 0.68, 0.66];
%   SIC      = 0.1;
%   Depth    = 3420;
%   DIC      = [2255, 2250, 2255, 2251, 2264, 2265, 2270, 2276, 2280, 2278];
%   TA       = [2365, 2361, 2366, 2361, 2374, 2377, 2380, 2386, 2393, 2397];
% % Sample d13C_cib data for inverse modelling:
%   d13C_cib = [0.33, 0.33, 0.32, 0.30, 0.31, 0.23, 0.23, 0.20, 0.12, 0.21];
%
% % Forward modelling 
% % Example 1: W2D model
%   [d13C_cib1,upper1,lower1] = BayBCar(d13C_dic,Depth);
% % Example 2: A4OIp model. Using BayBCar(d13C_dic,[],Omega,SIC,1) gives identical result.
%   [d13C_cib2,upper2,lower2] = BayBCar(d13C_dic,Depth,Omega,SIC,1);
% % Example 3: A4OIp model using CO2SYS
%   [d13C_cib3,upper3,lower3] = BayBCar(d13C_dic,Depth,[],SIC,1,[],DIC,TA);
%
% % Inverse modelling 
%   Example 4: A4OIp model
%   [d13C_dic4,upper4,lower4] = BayBCar(d13C_cib,Depth,Omega,SIC,1,"inverse");

%% Preprocessing the variables

Depth     = [];
Omega     = [];
SIC       = [];

% Error if d13C_x is not provided
if isempty(d13C_x)
    error('Error. d13C_x is not provided')
end

% If the last argument is a string, it specifies the direction of prediction
if ~isempty(varargin) && (isstring(varargin{end}) || ischar(varargin{end}))
    direction = string(varargin{end});
    varargin(end) = [];
else
    direction = "forward"; % forward modelling (if no direction specified)
end

% Error if direction is neither "forward" nor "inverse"
if direction ~= "forward" & direction ~= "inverse"
    error('Error. Direction must be either "forward" or "inverse"')
end

try
    Depth    = varargin{1};
    Omega    = varargin{2};
    SIC      = varargin{3};
    Species  = varargin{4};
    CI       = varargin{5};
end

% Default values
if ~exist('Species','var') | isempty(Species)
    Species = 2;  % C. wuellerstorfi
end
if ~exist('CI','var') | isempty(CI)
    CI = 0.95;    % 95% CI
end

% Make sure all variables are column vectors
if isrow(d13C_x)
    d13C_x = d13C_x';
end
if isrow(Depth)
    Depth = Depth';
end
if isrow(Omega)
    Omega = Omega';
end
if isrow(SIC)
    SIC = SIC';
end
if isrow(Species)
    Species = Species';
end

% Convert scalar variables to vectors (if needed)
if ~isempty(Depth) && isscalar(Depth)
    Depth = Depth * ones(size(d13C_x));
end
if ~isempty(Omega) && isscalar(Omega)
    Omega = Omega * ones(size(d13C_x));
end
if ~isempty(SIC) && isscalar(SIC)
    SIC = SIC * ones(size(d13C_x));
end
if ~isempty(Species) && isscalar(Species)
    Species = Species * ones(size(d13C_x));
end

%% If Omega is not available, attempt to calculate it using CO2SYS toolbox

if isempty(Omega) && ~isempty(Depth)

    if exist('CO2SYS', 'file') == 2

        DIC  = [];
        TA   = [];
        pH   = [];
        pCO2 = [];

        try
            DIC  = varargin{6};
            TA   = varargin{7};
            pH   = varargin{8};
            pCO2 = varargin{9};
            T    = varargin{10};
            S    = varargin{11};
            Si   = varargin{12};
            PO4  = varargin{13};
        end

        % Default values
        if ~exist('T','var') | isempty(T)
            T = 4;
        end
        if ~exist('S','var') | isempty(S)
            S = 35;
        end
        if ~exist('Si','var') | isempty(Si)
            Si = 1;
        end
        if ~exist('PO4','var') | isempty(PO4)
            PO4 = 1;
        end

        if sum([~isempty(DIC),~isempty(TA),~isempty(pH),~isempty(pCO2)]) >= 2

            PRS = Depth;

            % Make sure all variables are row vectors
            if iscolumn(PRS)
                PRS = PRS';
            end
            if iscolumn(DIC)
                DIC = DIC';
            end
            if iscolumn(TA)
                TA = TA';
            end
            if iscolumn(pH)
                pH = pH';
            end
            if iscolumn(pCO2)
                pCO2 = pCO2';
            end
            if iscolumn(T)
                T = T';
            end
            if iscolumn(S)
                S = S';
            end
            if iscolumn(Si)
                Si = Si';
            end
            if iscolumn(PO4)
                PO4 = PO4';
            end

            % Define constants for carbonate system calculations
            SCALE = 1;  % Total pH scale
            K1K2  = 10; % Leuker et al (2000) K1K2
            SO4   = 1;  % Dickson (1990) KSO4
            KF    = 2;  % Perez & Fraga (1987) KF
            BOR   = 2;  % Lee et al (2010) TB

            % Calculate from TA and DIC
            if ~isempty(TA) && ~isempty(DIC)
                Co2s = CO2SYS(TA,DIC,1,2,S,20,T,0,PRS,Si,PO4,0,0,SCALE,K1K2,SO4,KF,BOR);
                Omega = Co2s(:,35);

                % Calculate from TA and pH
            elseif ~isempty(TA) && ~isempty(pH)
                Co2s = CO2SYS(TA,pH,1,3,S,20,T,0,PRS,Si,PO4,0,0,SCALE,K1K2,SO4,KF,BOR);
                Omega = Co2s(:,35);

                % Calculate from DIC and pH
            elseif ~isempty(DIC) && ~isempty(pH)
                Co2s = CO2SYS(DIC,pH,2,3,S,20,T,0,PRS,Si,PO4,0,0,SCALE,K1K2,SO4,KF,BOR);
                Omega = Co2s(:,35);

                % Calculate from TA and pCO2
            elseif ~isempty(TA) && ~isempty(pCO2)
                Co2s = CO2SYS(TA,pCO2,1,4,S,20,T,0,PRS,Si,PO4,0,0,SCALE,K1K2,SO4,KF,BOR);
                Omega = Co2s(:,35);

                % Calculate from DIC and pCO2
            elseif ~isempty(DIC) && ~isempty(pCO2)
                Co2s = CO2SYS(DIC,pCO2,2,4,S,20,T,0,PRS,Si,PO4,0,0,SCALE,K1K2,SO4,KF,BOR);
                Omega = Co2s(:,35);

                % Calculate from pH and pCO2
            elseif ~isempty(pH) && ~isempty(pCO2)
                Co2s = CO2SYS(pH,pCO2,3,4,S,20,T,0,PRS,Si,PO4,0,0,SCALE,K1K2,SO4,KF,BOR);
                Omega = Co2s(:,35);
            end

            if isrow(Omega)
                Omega = Omega';
            end

        end

    else
        disp('CO2SYS toolbox not found - using depth instead of Omega')
    end
end

%% BayBCar modelling

% Output NaNs if neither Omega nor Depth are available
if isempty(Omega) & isempty(Depth)
    d13C_pred = nan(size(d13C_x));
    upper = nan(size(d13C_x));
    lower = nan(size(d13C_x));
    warning("Neither Omega nor Depth are available")
    return
end

% Load calibrated PSM parameters
PSMs = load('BayBCar_calibrations.mat');

if ~isempty(Omega) && ~isempty(SIC) && all(Species == 2)

    PSM_name = 'W3OI';
    mu = PSMs.mu_W3OI;
    Sigma = PSMs.Sigma_W3OI;
    a = PSMs.a_W3OI;
    b = PSMs.b_W3OI;

    X = [ones(size(d13C_x)), d13C_x, log(Omega), SIC];

elseif isempty(Omega) && ~isempty(SIC) && ~isempty(Depth) && all(Species == 2)

    PSM_name = 'W3DI';
    mu = PSMs.mu_W3DI;
    Sigma = PSMs.Sigma_W3DI;
    a = PSMs.a_W3DI;
    b = PSMs.b_W3DI;

    X = [ones(size(d13C_x)), d13C_x, Depth, SIC];

elseif ~isempty(Omega) && isempty(SIC) && all(Species == 2)

    PSM_name = 'W2O';
    mu = PSMs.mu_W2O;
    Sigma = PSMs.Sigma_W2O;
    a = PSMs.a_W2O;
    b = PSMs.b_W2O;

    X = [ones(size(d13C_x)), d13C_x, log(Omega)];

elseif isempty(Omega) && isempty(SIC) && ~isempty(Depth) && all(Species == 2)

    PSM_name = 'W2D';
    mu = PSMs.mu_W2D;
    Sigma = PSMs.Sigma_W2D;
    a = PSMs.a_W2D;
    b = PSMs.b_W2D;

    X = [ones(size(d13C_x)), d13C_x, Depth];

elseif ~isempty(Omega) && ~isempty(SIC) && any(Species ~= 2)

    PSM_name = 'A4OIp';
    mu = PSMs.mu_A4OIp;
    Sigma = PSMs.Sigma_A4OIp;
    a = PSMs.a_A4OIp;
    b = PSMs.b_A4OIp;

    X_Species = zeros(length(d13C_x),7);
    X_Species(sub2ind(size(X_Species), (1:length(Species))', Species)) = 1;
    X = [ones(size(d13C_x)), d13C_x, log(Omega), SIC, X_Species(:, [1,3:end])];

elseif isempty(Omega) && ~isempty(SIC) && ~isempty(Depth) && any(Species ~= 2)

    PSM_name = 'A4DIp';
    mu = PSMs.mu_A4DIp;
    Sigma = PSMs.Sigma_A4DIp;
    a = PSMs.a_A4DIp;
    b = PSMs.b_A4DIp;

    X_Species = zeros(length(d13C_x),7);
    X_Species(sub2ind(size(X_Species), (1:length(Species))', Species)) = 1;
    X = [ones(size(d13C_x)), d13C_x, Depth, SIC, X_Species(:, [1,3:end])];

elseif ~isempty(Omega) && isempty(SIC) && any(Species ~= 2)

    PSM_name = 'A3Op';
    mu = PSMs.mu_A3Op;
    Sigma = PSMs.Sigma_A3Op;
    a = PSMs.a_A3Op;
    b = PSMs.b_A3Op;

    X_Species = zeros(length(d13C_x),7);
    X_Species(sub2ind(size(X_Species), (1:length(Species))', Species)) = 1;
    X = [ones(size(d13C_x)), d13C_x, log(Omega), X_Species(:, [1,3:end])];

elseif isempty(Omega) && isempty(SIC) && ~isempty(Depth) && any(Species ~= 2)

    PSM_name = 'A3Dp';
    mu = PSMs.mu_A3Dp;
    Sigma = PSMs.Sigma_A3Dp;
    a = PSMs.a_A3Dp;
    b = PSMs.b_A3Dp;

    X_Species = zeros(length(d13C_x),7);
    X_Species(sub2ind(size(X_Species), (1:length(Species))', Species)) = 1;
    X = [ones(size(d13C_x)), d13C_x, Depth, X_Species(:, [1,3:end])];

end

if direction == "forward"
    disp(['forward modelling using ' PSM_name]);
    d13C_pred = X * mu;
elseif direction == "inverse"
    disp(['inverse modelling using ' PSM_name]);
    d13C_pred = (d13C_x - X(:,1:end ~= 2) * mu(1:end ~= 2)) / mu(2);
end

%% Confidence interval

if direction == "forward"
    % Calculate variance
    Var = (b / a) * (1 + diag(X * Sigma * X'));

elseif direction == "inverse"
    % Jacobian of the inverse prediction with respect to the coefficients
    J = zeros(size(X));
    J(:,2) = -d13C_pred / mu(2);    % partial derivative for mu(2)
    J(:,1:end ~= 2) = -X(:,1:end ~= 2) / mu(2);  % partial derivatives for other coefficients

    % Propagation of uncertainty
    Var = (b / a) * (1/mu(2)^2 + diag(J * Sigma * J'));
end

% Confidence intervals using Student's t-distribution
alpha = 1 - CI;
tval = tinv(1 - alpha/2, 2*a);
halfCI = tval * sqrt(Var);  % Confidence interval half-widths
upper = d13C_pred + halfCI;
lower = d13C_pred - halfCI;