%% This code calculates temperatures at depth through time for excess ice 
% given the parameters listed in a parameter file required to load in

% This code is based off Bramson et al. (2017, 2019). 
% includes treatment of atmospheric scattering and attenuation

%% 1. Read in parameter files and set up the variables and constants
clear;


% Read in parameter file:
paramFile = 'paramFile_45Lat.txt';
s = readStruct(paramFile);

% Grab some important values from the structure
f = s.impexp;  % Set variable for how explicit vs implicit to numerically integrate
dt = s.dt;
Tref = s.initSurfTemp; % First orbital timestep will use this reference temp, after that will use mean surface temperature from orbital run immediately before

% Load Laskar 2004 solutions (backwards so it is marching forward in time)
load Laskar2004Solutions_backwards2.mat; 
L04Lsperi = L04Lsp - pi; % Laskar calculates at aphelion

% Number of loops it will do in total
numLaskarRuns = ceil((s.Laskar_end+1-s.Laskar_begin)/s.Laskar_step);
eqDepths = zeros(1, numLaskarRuns);
dz_excessIceTimestep = zeros(1, numLaskarRuns);

flatVisSaved = 0; % Will get written over, just to have something to pass into calcOrbitalParams in the initial or flat cases where we haven't calculated a flatVisSaved or flatIRSaved
flatIRSaved = 0; % Will get written over, just to have something to pass into calcOrbitalParams in the initial or flat cases where we haven't calculated a flatVisSaved or flatIRSaved

%% 2. Setting variables
sigmasb = 5.6705119e-8;          % Stefan Boltzmann constant W m^-2 K^-4
avogadro = 6.022140857e23; % molecules/mole
molarMassH2O = 0.01801528; % kg/mole
boltzmann = 1.38064852e-23; % Boltzmann Constant, m2 kg s-2 K-1
densityIce = s.density_ice; % 917 is non-porous, from Dundas et al. 2015, or just use value in param file
PtoRhoFactor = (molarMassH2O / avogadro) / boltzmann; % to make go faster
SF12_a1 = -1.27; % Parameterization from Schorghofer and Forget (2012)
SF12_b1 = 0.139; % Parameterization  from Schorghofer and Forget (2012)
SF12_c1 = -0.00389; % Parameterization from Schorghofer and Forget (2012)
porosity_param = (s.iceDust)*(1/(1-s.icePorosity-s.iceDust))/(1-s.regolithPorosity);

%% 3. Set up file and folder for storing output 
% output - L04 step, time before present, T at surface, then T through
% depth?,

% Will store outputs in folders organized by the date it was ran within an
% Output folder. Each file that is saved will contain the same string in
% its filename with date and time it was run
if exist('Output', 'dir') == 0
    mkdir('Output');
end
cd Output
if exist(datestr(now, 'yyyy-mm-dd'), 'dir') == 0
    mkdir(datestr(now, 'yyyy-mm-dd'));
end
cd ..

strFileLocationParam = strcat('Output/', datestr(now, 'yyyy-mm-dd'), '/', datestr(now, 'yyyy-mm-dd_HHMM'));

if s.saveParamFile == 'y'
    copyfile(paramFile, strcat(strFileLocationParam, '_ParamFile.txt'));
end

fileName = strcat(strFileLocationParam, s.outputSaveFile);   % if want to input file name manually, use something like 'Folder/file.txt'
fileID = fopen(fileName,'w');
% Write a header to the output file with each value that will be saved
% later in a tab delineated text file
fprintf(fileID, ['L04Step\t' 'timebeforepresent\t' ...
        'meanTsurf\t' 'iceTableTemps_mean\t' '\n' ]);
    
formatSpec = repmat('%20.9f\t', 1, 152);
fmt = ['%6.0f\t%12.0f\t' ...
        '%20.9f\t' ...
        formatSpec ...
        '\n'
];

fclose(fileID);

%% 4. Get the numerical grid of layers set up
% Calculate thermophysical properties for dry lag, pf, and excess ice
% This is run just once at the beginning to get vectors of k_input,
% rhoc_input, diurnalSkinDepths, annualSkinDepths and allKappas
layerProperties_setup; % This just run once at beginning

% Prepare vectors for flat cases (generally excess ice or whatever
% stratigraphy is input
k_flat = [s.flat_TI_input(1)*s.flat_TI_input(1)/s.rhoc_flat_input(1); s.flat_TI_input(2)*s.flat_TI_input(2)/s.rhoc_flat_input(2); s.flat_TI_input(3)*s.flat_TI_input(3)/s.rhoc_flat_input(3)];
rhoc_flat = [s.rhoc_flat_input(1); s.rhoc_flat_input(2); s.rhoc_flat_input(3)];
allKappas_flat = k_flat ./ rhoc_flat;
diurnalSkinDepths_flat = sqrt(allKappas.*lengthOfDay/pi); % meters
annualSkinDepths_flat = sqrt(allKappas.*lengthOfYear/pi); % meters

% Get initial lag thickness - assumes no pore filling ice at start
% Either use a user-specific value or performs initial run to determine the
% equilibrium depth of ice for the first timestep and start with a dry lag
% of that thickness on a flat terrain with main albedo
if s.useInputDepth == 'n'
    z_ei = annualSkinDepths(1)*annualLayers; % look within 6 (or whatever annualLayers is) skin depths of dry lag and have no pore filling ice
    z_pf = annualSkinDepths(1)*annualLayers; % no pore-filling ice- set pf value to be same as z_ei
    [k, rhoc, Kappa, dz, z_ei_index, z_pf_index, depthsAtLayerBoundaries, depthsAtMiddleOfLayers, nLayers] = setLayerProperties(z_ei, z_ei, k_input, rhoc_input, diurnalSkinDepths, annualSkinDepths, s.layerGrowth, s.dailyLayers, s.annualLayers);

    ecc = L04ecc(s.Laskar_begin);
    Lsp = L04Lsperi(s.Laskar_begin);
    obl = L04obl(s.Laskar_begin);
    
    oblDeg = rad2deg(obl);
    if oblDeg > 10 && oblDeg < 28
        atmosPressSF12 = exp(SF12_a1 + SF12_b1*(oblDeg-28));
    elseif oblDeg > 28 && oblDeg < 50
        atmosPressSF12 = exp(SF12_a1 + SF12_b1*(oblDeg-28) + SF12_c1*(oblDeg-28)^2);
    end
    atmosPressSF12 = atmosPressSF12 * s.atmosFactor;
    
    slope = 0;
    slope_aspect = 0;
    albedo = s.albedo;
    
    findEquilibriumDepth;
    
    if PFICE == 0
        % case where no equilibrium depth found within several annual skin
        % depths
        fprintf('No equilibrium depth found. Ice not stable at all. Setting excess ice interface at 6 meters depth.');
        z_ei = 6;
        z_pf = 6;
    else 
        z_ei = z_eq; % Will set initial thickness of dry lag to be equilibrium depth
        z_pf = z_eq; % No pore-filling ice to start, excess ice starts at equilibrium depth, z_pf >= z_ei is that condition
    end
    
elseif s.useInputDepth == 'y'
    z_ei = s.initialLagThickness;
    z_pf = s.initialLagThickness;
else
    fprintf('Error in parameter file specifying which depth to use');
end

% This sets up the grids, will be run each time needs to change layer thicknesses
[k, rhoc, Kappa, dz, z_ei_index, z_pf_index, depthsAtLayerBoundaries, depthsAtMiddleOfLayers, nLayers] = setLayerProperties(z_ei, z_pf, k_input, rhoc_input, diurnalSkinDepths, annualSkinDepths, s.layerGrowth, s.dailyLayers, s.annualLayers);

%% 5. Run model throughout orbital timesteps in Laskar's solutions
counterOrbital = 1;
for L04Step = s.Laskar_begin:s.Laskar_step:s.Laskar_end
    
    tic
    ecc = L04ecc(L04Step);
    Lsp = L04Lsperi(L04Step);
    obl = L04obl(L04Step);
    timebeforepresent = L04time(L04Step);
    fprintf('\nTesting %12.0f years ago when obliquity was %12.9f, ecc was %12.9f and Lsp was %12.9f.\n', timebeforepresent, rad2deg(obl), ecc, rad2deg(Lsp));
        
    % Calculate atmospheric water vapor density based on Schorghofer and
    % Forget 2012 scheme, which just uses obliquity. Get the atmospheric
    % partial pressure here and then convert to density once mean surface
    % temperature is calculated
    oblDeg = rad2deg(obl);
    if oblDeg > 10 && oblDeg < 28
        atmosPressSF12 = exp(SF12_a1 + SF12_b1*(oblDeg-28));
    elseif oblDeg > 28 && oblDeg < 50
        atmosPressSF12 = exp(SF12_a1 + SF12_b1*(oblDeg-28) + SF12_c1*(oblDeg-28)^2);
    end
    atmosPressSF12 = atmosPressSF12 * s.atmosFactor;
    
    %% 6a. Calculate flat terrain first if a sloped case
    % It will output values needed to flatIRSaved and flatVisSaved
    if s.slope ~= 0 && s.slope_rerad == 1
        fprintf('\nRunning the no slope case first.\n');
        % For v12, since I know I'm using this for the NPLD, will make a
        % purely icy stratigraphy for calculating surroundings
        z_ei_old = z_ei;
        z_pf_old = z_pf;
        z_ei = 0;
        z_pf = 0;
        
        [k, rhoc, Kappa, dz, z_ei_index, z_pf_index, depthsAtLayerBoundaries, depthsAtMiddleOfLayers, nLayers] = setLayerProperties(z_ei, z_pf, k_flat, rhoc_flat, diurnalSkinDepths_flat, annualSkinDepths_flat, s.layerGrowth, s.dailyLayers, s.annualLayers);
        
        slope = 0;
        slope_aspect = 0;
        albedo = s.flatAlbedo;
        
        % Calculate orbital parameters and insolation values
        [sf, IRdown, visScattered, flatVis, flatIR, sky, CO2FrostPoints, nStepsInYear, lsWrapped, hr, ltst] = calcOrbitalParams(slope, slope_aspect, ecc, obl, Lsp, s, dt, flatVis, flatIR);

        % Run tempCalculator
        [Temps, Tsurf, frostMasses, flatIRSaved, flatVisSaved, iceTableTConverge] = tempCalculator(nLayers, nStepsInYear, s, dt, k, dz, rhoc, f, Kappa, sf, visScattered, flatVis, albedo, IRdown,flatIR, sky, CO2FrostPoints, depthsAtMiddleOfLayers, sigmasb, slope, Tref, z_ei_index);
        mean(Tsurf)
        %reorgDiurnal;
        
        fprintf('step 6a');
        z_ei = z_ei_old;
        z_pf = z_pf_old;
    end
    
    
    % If there is a minimum h0 lag value for this orbital timestep, and the
    % lag is supposed to be below that, set it at that h0 - only do this
    % though if param file has lag allowed to change in some capacity
    if s.growLag == 'y' || s.removeLag > 0
        if h0 > 0  % Is there an h0 value
            if z_ei < h0
                % Set to h0 if there is one and z_ei drops below that
                z_ei = h0;
            end
        else
            if z_ei < s.thinnestLayer % Somewhat arbitrary but want to make sure the code doesn't try to use a super small lag in the cases where there is no minimum h0 value
                % Just set to 0, this should be the cases where convective loss
                % would be pretty small anyways...let's see how it goes
                z_ei = 0;
            end
        end
    end
    
    fprintf('\nRunning the the sloped case now.\n');
    %% 6b. Calculate equilibrium depth, set stratigraphy and calculate temperatures
    if z_ei > 0
        % There is a lag deposit, so calculate these things while solving
        % for the equilibrium depth of pore filling ice in the
        % findEquilibriumDepth script but first set up layers
        [k, rhoc, Kappa, dz, z_ei_index, z_pf_index, depthsAtLayerBoundaries, depthsAtMiddleOfLayers, nLayers] = ...
            setLayerProperties(z_ei, z_pf, k_input, rhoc_input, diurnalSkinDepths, annualSkinDepths, s.layerGrowth, s.dailyLayers, s.annualLayers);
    
        slope = s.slope;
        slope_aspect = s.slope_aspect;
        albedo = s.albedo;
        
        findEquilibriumDepth
        
        %reorgDiurnal;
    else
        % No lag - only have excess ice all the way to the surface
       
        z_ei = 0;
        z_pf = 0;
        PFICE = 0;
        
        [k, rhoc, Kappa, dz, z_ei_index, z_pf_index, depthsAtLayerBoundaries, depthsAtMiddleOfLayers, nLayers] = setLayerProperties(z_ei, z_pf, k_input, rhoc_input, diurnalSkinDepths, annualSkinDepths, s.layerGrowth, s.dailyLayers, s.annualLayers);
        
        % Calculate temperatures!
        slope = s.slope;
        slope_aspect = s.slope_aspect;
        albedo = s.albedo;
        
        % Calculate orbital parameters and insolation values
        [sf, IRdown, visScattered, flatVis, flatIR, sky, CO2FrostPoints, nStepsInYear, lsWrapped, hr, ltst] = calcOrbitalParams(slope, slope_aspect, ecc, obl, Lsp, s, dt, flatVisSaved, flatIRSaved);

        % Run tempCalculator
        [Temps, Tsurf, frostMasses, flatIR, flatVis, iceTableTConverge] = tempCalculator(nLayers, nStepsInYear, s, dt, k, dz, rhoc, f, Kappa, sf, visScattered, flatVisSaved, albedo, IRdown, flatIRSaved, sky, CO2FrostPoints, depthsAtMiddleOfLayers, sigmasb, slope, Tref, z_ei_index);

        %reorgDiurnal;
    end
    
    meanTsurf = mean(Tsurf);
    iceTableTemps_mean = mean(Temps'); % mean temperature at the ice table over the year
    
        
    %% 9. Save outputs
    %{
    fileID = fopen(fileName,'a');
    fprintf(fileID, fmt, [L04Step timebeforepresent ecc oblDeg rad2deg(Lsp) Tsurf_mean iceTableTemps_mean atmosDensitySF12 iceTable_rhov_mean groundicecounter-1 eqDepths(counterOrbital) z_eqConverge PFICE tempTracker1 tempTracker2 tempTracker3 delta_lagThickness old_lagThickness dz_excessIceTimestep(counterOrbital) dz_excessIceCumulate(counterOrbital)]);
    fclose(fileID);
    %}
    fileID = fopen(fileName,'a');
    fprintf(fileID, fmt, [L04Step timebeforepresent ...
        meanTsurf iceTableTemps_mean]);
    fclose(fileID);
    
    fprintf('\nYoure %5.0f / %5.0f of the way there for orbital solutions!\n\n', counterOrbital, numLaskarRuns);
    counterOrbital = counterOrbital + 1;
    
    toc
end