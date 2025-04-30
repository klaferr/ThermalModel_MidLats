% Calculates k and rho*c based on porosities and dust content of regolith
% and excess ice inputted in the param file using endmember values
% OR (new to v11)
% Uses values given in the param file (if set to non-0 values)
% Also in v11, doesn't set up the initial grid anymore with layer
% thicknesses, as that is done in setLayerProperties
% 
% October 29, 2018: ensures linear conductivity for ice and rock
% if no porosity in excess ice and has various options commented out that
% have been tested in the past including the Van Dusen 1929 non-linear
% dependence on density density

au  = 1.4959787061e11;           % Size of one astronomical unit in meters
gc  = 6.67259e-11;               % Gravitational constant
sm  = 1.9891e30;                 % Solar mass
lengthOfYear = 2*pi/sqrt(gc*sm/(s.a*au)^3);   % length of Mars year in seconds using Kepler's 3rd law

if s.rhoc_input(1) == 0
    rhoc_reg = s.density_rock*s.c_rock*(1-s.regolithPorosity);
else
    rhoc_reg = s.rhoc_input(1);
end

if s.rhoc_input(2) == 0
    rhoc_pf = s.density_rock*s.c_rock*(1-s.regolithPorosity) + s.density_ice*s.c_ice*(s.regolithPorosity);
else
    rhoc_pf = s.rhoc_input(2);
end

if s.rhoc_input(3) == 0
    rhoc_excice = s.density_rock*s.c_rock*s.iceDust + s.density_ice*s.c_ice*(1-s.iceDust-s.icePorosity);
else
    rhoc_excice = s.rhoc_input(3);
end

k_reg = s.TI_input(1) * s.TI_input(1) / rhoc_reg; % Set by thermal inertia and rho*c

if s.TI_input(2) == 0
    k_pf = (1-s.regolithPorosity)*s.k_rock + s.regolithPorosity*s.k_ice;
else
    k_pf = s.TI_input(2) * s.TI_input(2) / rhoc_pf; % Set by thermal inertia and rho*c
end
    
if s.TI_input(3) == 0
    if s.icePorosity ~= 0
        % Linear with porosity, for a given dust content
        k_excice = (s.iceDust/(1-s.icePorosity))*s.k_rock + ((1-s.iceDust-s.icePorosity)^2/(1-s.icePorosity))*s.k_ice;
                    
        % Different, non-linear scheme with porosity for a given dust content
        % k_excice = (s.iceDirt/(1-s.icePorosity))*s.k_rock + ((1-s.iceDirt-s.icePorosity)/(1-s.icePorosity))*s.k_ice;
        
        % Van Dusen (1929) equation, non-linear scheme with density
        %density_tmp = s.density_rock*s.iceDirt + s.density_ice.*(1-s.iceDirt-s.icePorosity);
        %k_excice = 2.1e-2 + (4.2e-4).*density_tmp + (2.2e-9).*density_tmp.^3;
    else
        % Linear combination of ice and dirt if no porosity
        k_excice = s.iceDust*s.k_rock + (1-s.iceDust)*s.k_ice;   
    end
else
    k_excice = s.TI_input(3) * s.TI_input(3) / rhoc_excice; % Set by thermal inertia and rho*c
end

k_input = [k_reg; k_pf; k_excice];
rhoc_input = [rhoc_reg; rhoc_pf; rhoc_excice];
allKappas = k_input ./ rhoc_input;

lengthOfDay = s.rot;
layerGrowth = s.layerGrowth;
dailyLayers = s.dailyLayers;
annualLayers = s.annualLayers;

fprintf ('~~~Layer Properties~~~\n');
for ii = 1:size(k_input,1)
    I = sqrt(k_input(ii)*rhoc_input(ii)); % Thermal inertia in J/m^-2 K^-1 s^-1/2
    fprintf ('Thermal Inertia of the %1.0f layer = %5.2f [J/(m^2 K s^1/2)] \n', ii, I);
end

diurnalSkinDepths = sqrt(allKappas.*lengthOfDay/pi); % meters
annualSkinDepths = sqrt(allKappas.*lengthOfYear/pi); % meters
fprintf ('Diurnal thermal skin depth of all layers : %.8f m, %.8f m, %.8f m\n', diurnalSkinDepths(1),  diurnalSkinDepths(2),  diurnalSkinDepths(3));
fprintf ('Annual thermal skin depth of all layers : %.8f m, %.8f m, %.8f m\n', annualSkinDepths(1),  annualSkinDepths(2),  annualSkinDepths(3));