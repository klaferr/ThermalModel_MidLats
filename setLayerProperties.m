function [k_output, rhoc_output, Kappa_output, dz, z_ei_index, z_pf_index, depthsAtLayerBoundaries, depthsAtMiddleOfLayers, nLayers] = setLayerProperties(z_ei, z_pf, k_input, rhoc_input, diurnalSkinDepths, annualSkinDepths, layerGrowth, dailyLayers, annualLayers)
%setLayerProperties Sets up grid of numerical layers for the thermal model
%   Makes a grid of numerical layers, determines the layer thicknesses and
%   sets their thermophysical properties according to the stratigraphy
%   
%   Updated on May 13, 2018 for pore filling ice up to the surface to still use
%   dry lag skin depth for setting first layer thickness, otherwise won't
%   necessarily get layers thin enough if the pore filling ice/lag is on
%   the order of millimeters thick

z_ei
z_pf 

if z_ei <= 0
    % This is the case where there is no lag at all, all layers are set to
    % thermophysical properties of ice. Layering thicknesses are calculated in
    % standard way but with 15 layers in ice's diurnal skin depth
    fprintf('Here 1');
    firstLayerThickness = diurnalSkinDepths(3) / ((1-layerGrowth^dailyLayers)/(1-layerGrowth)); % Thickness of first layer based on diurnal skin depth of excess ice
    nLayers = ceil(log(1-(1-layerGrowth)*(annualLayers*annualSkinDepths(end)/firstLayerThickness))/log(layerGrowth) ); % Number of subsurface layers based on annual skin depth of deepest [excess ice] layer
    dz = (firstLayerThickness * layerGrowth.^[0:nLayers-1])'; % transpose to make column vector
    
    % Set up property vectors
    k_output = dz.*0;           % Thermal conductivities (W m^-1 K^-1)
    rhoc_output = dz.*0;     % Densities of subsurface (kg/m^3) * Specific heats J/(kg K)
    
    % Set all to excess ice values
    k_output(1:end) = k_input(3);
    rhoc_output(1:end) = rhoc_input(3);
    z_ei_index = 0; % Excess ice layer starts at first layer, so basically doesn't make sense to have positive integer index for where it starts
    z_pf_index = -1; % No pore-filling ice
    
    z_ei
    z_ei_index
    z_pf
    z_pf_index
else
    % Every case where there is a lag deposit (z_ei is not 0)
    fprintf('Here 2');
    firstLayerThickness_orig = diurnalSkinDepths(1) / ((1-layerGrowth^dailyLayers)/(1-layerGrowth)); % Thickness of first layer based on diurnal skin depth of surface [porous regolith] material
    
    if z_pf <= 0
        fprintf('Here 3');
        % Case where pore-filling ice is stable up to the surface
        % Differs from above z_ei = 0 where there is no lag deposit so
        % the properties are excess ice all the way up. In this case we'll
        % still have an excess ice interface
        % Don't use skin depth of pore filling ice to calculate dz because
        % skin depth is a lot bigger so may not have a first layer thin
        % enough if the pf ice layer is on the order of millimeters!
        firstLayerThickness = diurnalSkinDepths(1) / ((1-layerGrowth^dailyLayers)/(1-layerGrowth)) % Thickness of first layer based on diurnal skin depth of surface [porous regolith] material
        nLayers = ceil(log(1-(1-layerGrowth)*(annualLayers*annualSkinDepths(3)/firstLayerThickness))/log(layerGrowth) ); % Number of subsurface layers based on annual skin depth of deepest [excess ice] layer
        dz = (firstLayerThickness * layerGrowth.^[0:nLayers-1])'; % transpose to make column vector
        
        % Set up property vectors
        k_output = dz.*0;           % Thermal conductivities (W m^-1 K^-1)
        rhoc_output = dz.*0;     % Densities of subsurface (kg/m^3) * Specific heats J/(kg K)
        
        % SET EXCESS ICE TABLE DEPTH
        % Snaps layer boundary closest to excess ice interface to that exact depth
        depthsAtLayerBoundaries = cumsum(dz);
        [tmp, z_ei_index] = min(abs(depthsAtLayerBoundaries - z_ei))
        
        % Check if z_ei is also really small (but slightly bigger than
        % z_pf, hence why it got to this part of the if/elseif) and the
        % closest index is 1, then change it so the 2nd layer is the
        % pore filling ice
        if z_ei_index == 1
            z_ei_index = 2;
        end
        
        z_pf_index = 0; % pore filling ice all the way to surface
        
        z_ei
        z_ei_index
        z_pf
        z_pf_index
        
        depthsAtLayerBoundaries(z_ei_index) = z_ei;
        
        % Remesh those layers that were affected
        dz(z_ei_index) = abs(depthsAtLayerBoundaries(z_ei_index) - depthsAtLayerBoundaries(z_ei_index-1));
        dz(z_ei_index+1) = abs(depthsAtLayerBoundaries(z_ei_index+1) - depthsAtLayerBoundaries(z_ei_index));
         
        % Set the rest of the thermophysical layer parameters
        k_output(1:z_ei_index) = k_input(2);
        rhoc_output(1:z_ei_index) = rhoc_input(2);
        k_output(z_ei_index+1:end) = k_input(3);
        rhoc_output(z_ei_index+1:end) = rhoc_input(3);
        
    elseif z_pf < firstLayerThickness_orig || z_ei < firstLayerThickness_orig
        fprintf('Here 4');
        % Case where will just have one thin layer of dry regolith and the
        % thickness of this layer is supposed to be thinner than what the
        % first layer thickness is in the standard scheme
        
        % Set first layer equal to thickness of the thinnest layer and grow by
        % 1.03 until 6 annual skin depths of ice
        firstLayerThickness = min(z_pf, z_ei);
        nLayers = ceil(log(1-(1-layerGrowth)*(annualLayers*annualSkinDepths(3)/firstLayerThickness))/log(layerGrowth) ); % Number of subsurface layers based on annual skin depth of deepest [excess ice] layer
        dz = (firstLayerThickness * layerGrowth.^[0:nLayers-1])'; % transpose to make column vector
        
        % Set up property vectors
        k_output = dz.*0;           % Thermal conductivities (W m^-1 K^-1)
        rhoc_output = dz.*0;     % Densities of subsurface (kg/m^3) * Specific heats J/(kg K)
        
        if firstLayerThickness == z_ei
            fprintf('Here 5');
            % This is the case that have one thin layer of dry lag with no pf ice
            
            z_ei_index = 1; % Excess ice layer starts after first layer, depth to start of excess ice will be at depthAtLayerBoundaries(z_ei_index)
            z_pf_index = -1; % No pore-filling ice
            
            z_ei
            z_ei_index
            z_pf
            z_pf_index
            
            k_output(1) = k_input(1);
            rhoc_output(1) = rhoc_input(1);
            k_output(2:end) = k_input(3);
            rhoc_output(2:end) = rhoc_input(3);
            
        elseif firstLayerThickness == z_pf
            fprintf('Here 6');
            % This is the case that have one thin layer of dry lag with pf ice in
            % between the equilibrium depth and excess ice interface z_ei
            k_output(1) = k_input(1);
            rhoc_output(1) = rhoc_input(1);
            z_pf_index = 1; % pore filling ice starts after first layer, first layer is dry lag
            
            % SET EXCESS ICE TABLE DEPTH
            % Snaps layer boundary closest to excess ice interface to that exact depth
            depthsAtLayerBoundaries = cumsum(dz);
            [tmp, z_ei_index] = min(abs(depthsAtLayerBoundaries - z_ei));
            
            % Check if z_ei is also really small (but slightly bigger than
            % z_pf, hence why it got to this part of the if/elseif) and the
            % closest index is 1, then change it so the 2nd layer is the
            % pore filling ice
            if z_ei_index == 1
                z_ei_index = 2;
            end
            
            z_ei
            z_ei_index
            z_pf
            z_pf_index
            
            depthsAtLayerBoundaries(z_ei_index) = z_ei;
            
            % Remesh those layers that were affected
            dz(z_ei_index) = abs(depthsAtLayerBoundaries(z_ei_index) - depthsAtLayerBoundaries(z_ei_index-1));
            dz(z_ei_index+1) = abs(depthsAtLayerBoundaries(z_ei_index+1) - depthsAtLayerBoundaries(z_ei_index));

            % Set the rest of the thermophysical layer parameters
            k_output(z_pf_index+1:z_ei_index) = k_input(2);
            rhoc_output(z_pf_index+1:z_ei_index) = rhoc_input(2);
            k_output(z_ei_index+1:end) = k_input(3);
            rhoc_output(z_ei_index+1:end) = rhoc_input(3);
        else
            fprintf('Shouldn''t have gotten here. Error in calculating first layer thickness.');
        end
    else
        fprintf('Here 7');
        % This is the case where will have more than one layer for the dry
        % lag deposit, as it is thicker than the first layer thickness from
        % the standard scheme. This is normal cases where use standard dz
        % scheme and can just snap pf and ei interfaces to nearest layer
        % boundary.
        
        firstLayerThickness_orig = diurnalSkinDepths(1) / ((1-layerGrowth^dailyLayers)/(1-layerGrowth)); % Thickness of first layer based on diurnal skin depth of surface [porous regolith] material
        
        nLayers = ceil(log(1-(1-layerGrowth)*(annualLayers*annualSkinDepths(3)/firstLayerThickness_orig))/log(layerGrowth) ); % Number of subsurface layers based on annual skin depth of deepest [excess ice] layer
        dz = (firstLayerThickness_orig * layerGrowth.^[0:nLayers-1])'; % transpose to make column vector
        depthsAtLayerBoundaries = cumsum(dz);
        [tmp, z_pf_index] = min(abs(depthsAtLayerBoundaries - z_pf));
        [tmp, z_ei_index] = min(abs(depthsAtLayerBoundaries - z_ei));
        
        z_ei
        z_ei_index
        z_pf
        z_pf_index
        
        % Set up property vectors
        k_output = dz.*0;           % Thermal conductivities (W m^-1 K^-1)
        rhoc_output = dz.*0;     % Densities of subsurface (kg/m^3) * Specific heats J/(kg K)
        
        % Only set layers for pore filling ice if z_pf is shallower than
        % z_ei, otherwise pore-filling ice not stable within the thickness
        % of the lag     
        if z_pf < z_ei
            fprintf('Here 8');
            % Want to check for case where pore-filling ice is barely stable
            % and be within the same layer as the excess ice interface. If that
            % is the case, move the adjacent layer to the excess ice
            % interface so that there is one layer of pore filling ice of the
            % appropriate thickness.
            if z_pf_index == z_ei_index
                z_ei_index = z_ei_index + 1;
            end
            
            % SET PORE FILLING ICE INTERFACE DEPTH
            % Snaps layer boundary closest to pore filling interface to that exact depth
            depthsAtLayerBoundaries(z_pf_index) = z_pf;
            
            z_ei
            z_ei_index
            z_pf
            z_pf_index
            
            % Remesh those layers that were affected by pf layer
            % Check if index is 1, then can't use
            % depthsAtLayerBoundaries(index - 1)
            if z_pf_index > 1
                dz(z_pf_index) = abs(depthsAtLayerBoundaries(z_pf_index) - depthsAtLayerBoundaries(z_pf_index-1));
                dz(z_pf_index+1) = abs(depthsAtLayerBoundaries(z_pf_index+1) - depthsAtLayerBoundaries(z_pf_index)); 
            else
                dz(z_pf_index) = z_pf; % Since in this case, would be subtracting off 0 essentially
                dz(z_pf_index+1) = abs(depthsAtLayerBoundaries(z_pf_index+1) - depthsAtLayerBoundaries(z_pf_index)); 
            end
            depthsAtLayerBoundaries = cumsum(dz);
            
            % SET EXCESS ICE INTERFACE DEPTH
            % Snaps layer boundary closest to excess ice interface to that exact depth
            depthsAtLayerBoundaries(z_ei_index) = z_ei;
            
            % Remesh those layers that were affected by ei layer
            dz(z_ei_index) = abs(depthsAtLayerBoundaries(z_ei_index) - depthsAtLayerBoundaries(z_ei_index-1));
            dz(z_ei_index+1) = abs(depthsAtLayerBoundaries(z_ei_index+1) - depthsAtLayerBoundaries(z_ei_index));
            
            % Set the rest of the thermophysical layer parameters
            k_output(1:z_pf_index) = k_input(1);
            rhoc_output(1:z_pf_index) = rhoc_input(1);
            k_output(z_pf_index+1:z_ei_index) = k_input(2);
            rhoc_output(z_pf_index+1:z_ei_index) = rhoc_input(2);
            k_output(z_ei_index+1:end) = k_input(3);
            rhoc_output(z_ei_index+1:end) = rhoc_input(3);
            
        else
            % If no pore-filling ice, set the index to -1
            fprintf('Here 9');
            z_pf_index = -1;
            
            % SET EXCESS ICE INTERFACE DEPTH
            % Snaps layer boundary closest to excess ice interface to that exact depth
            depthsAtLayerBoundaries(z_ei_index) = z_ei;
            
            z_ei
            z_ei_index
            z_pf
            z_pf_index
            
            % Remesh those layers that were affected by ei layer
            % Check if index is 1, then can't use
            % depthsAtLayerBoundaries(index - 1)
            if z_ei_index > 1
                dz(z_ei_index) = abs(depthsAtLayerBoundaries(z_ei_index) - depthsAtLayerBoundaries(z_ei_index-1));
                dz(z_ei_index+1) = abs(depthsAtLayerBoundaries(z_ei_index+1) - depthsAtLayerBoundaries(z_ei_index)); 
            else
                dz(z_ei_index) = z_ei; % Since in this case, would be subtracting off 0 essentially
                dz(z_ei_index+1) = abs(depthsAtLayerBoundaries(z_ei_index+1) - depthsAtLayerBoundaries(z_ei_index)); 
            end
            
            % Set the rest of the thermophysical layer parameters, with no
            % pore-filling layer
            k_output(1:z_ei_index) = k_input(1);
            rhoc_output(1:z_ei_index) = rhoc_input(1);
            k_output(z_ei_index+1:end) = k_input(3);
            rhoc_output(z_ei_index+1:end) = rhoc_input(3);
        end
        
    end
end


% Layer depths based on dz
depthsAtLayerBoundaries = cumsum(dz);
depthsAtMiddleOfLayers = cumsum(dz) - dz/2;

% Thermal Diffusivity, units m2/s
Kappa_output = k_output ./ rhoc_output;

end


