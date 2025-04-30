% Iterates through, finding the equilibrium depth of ice and setting
% pore-filling ice interface at that depth if it is within the thickness of
% the lag deposit

% Not implemented yet...but thinking about On April 16, 2018 update to change lag thickness so that if the z_eq
% interpolation gives a value below s.thinnestLayer, it sets it equal to
% s.thinnestLayer to avoid running into the possibility of a really thin,
% numerically unstable run.
% May 6, 2018: But, this case never really happens to the
% point that it makes the code freak out, while instead, forcing the lag to
% be a certain thickness can cause issues of Nans/Inf in the rho_v's for
% interpolating to get at z_eq. This caused issues at timestep associated
% with 921000 years ago for lags between 5-8 mm thick.
% So, from May 6 onward using this version of findEquilibriumDepth with the
% lines forcing z_pf to be at minimum s.thinnestLayer commented out for the
% cluster runs...which is why 1 Myr runs didn't work.

% Set variables to iterate through equilibrium pore filling ice
% solutions
FOUNDZeq = 0; % is the condition for breaking out of loop to find z_eq
groundicecounter = 1; % keeps track of how many iterations it takes
PFICE = 1; % tracker of if ground ice is stable (1) or not (0)

% Variable to make sure it checks that no ground ice is stable at least
% twice iterations in a row for finding the equilibrium depth
triedTwice = 0;

while FOUNDZeq == 0
    
    %{
    % if want to make sure the first layer of pore filling ice stays above
    this thickness...so far hasn't crashed the code occasionally having
    fraction of a mm thick pf ice but might consider using this at some
    point - no see above, using this causes more issues
    if z_pf <= s.thinnestLayer
        z_pf = s.thinnestLayer
    end
    %}
    
    %% Generate setup for the model
    [k, rhoc, Kappa, dz, z_ei_index, z_pf_index, depthsAtLayerBoundaries, depthsAtMiddleOfLayers, nLayers] = setLayerProperties(z_ei, z_pf, k_input, rhoc_input, diurnalSkinDepths, annualSkinDepths, s.layerGrowth, s.dailyLayers, s.annualLayers);
    dz_test = dz; % This will be the one we compare against until we settle on a solution
    % fprintf('z_ei = %8.3f, z_ei_index = %3.0f and depth at index is %8.3f\n', z_ei, z_ei_index, depthsAtLayerBoundaries(z_ei_index));
    
    % Last layer before ice table (last regolith/pf ice layer) is at same index
    % as for LayerBoundaries(excessIceTableIndex)
    % since layer boundaries come after the actual layers.
    % Regolith layers are at 1:excessIceTableIndex
    % Excess ice layers are at excessIceTableIndex+1:end
    % The depth of the interface between the layers is at
    % depthsAtLayerBoundaries(excessIceTableIndex).
    % ------------------- layer 1 below this
    % --------1---------- layer 2 below this
    % --------2----------
    % -------...--------- layer for excessIceTableIndex below this
    % ----excess ice ---- so the top excess ice layer is at index +1
    % though the depths is at ...LayerBoundaries(index)
    
    % To use for convergence criteria
    depthsAtLayerBoundaries_test = cumsum(dz_test);
    
    % Run calcOrbitalParams
    [sf, IRdown, visScattered, flatVis, flatIR, sky, CO2FrostPoints, nStepsInYear, lsWrapped, hr, ltst] = ...
        calcOrbitalParams(slope, slope_aspect, ecc, obl, Lsp, s, dt, flatVisSaved, flatIRSaved);

    % Run tempCalculator
    [Temps, Tsurf, frostMasses, flatIRSaved0, flatVisSaved0, iceTableTConverge] = ...
        tempCalculator(nLayers, nStepsInYear, s, dt, k, dz, rhoc, f, Kappa, sf, visScattered, flatVisSaved, albedo, IRdown,flatIRSaved, sky, CO2FrostPoints, depthsAtMiddleOfLayers, sigmasb, slope, Tref, z_ei_index);
        
    
    %% Figure out ground ice in lag deposit
    % Calculate estimates of vapor pressure above ice sheet in
    
    % Calculates vapors densities at layer boundaries to see if all vapor
    % densities below the ice table are higher or lower than the
    % atmospheric value
    onesMatrix = ones(1,size(Temps(2:z_ei_index+1,:),2));
    boundaryTemps = (k(2:z_ei_index+1).*dz(1:z_ei_index)*onesMatrix.*Temps(2:z_ei_index+1,:) + k(1:z_ei_index).*dz(2:z_ei_index+1)*onesMatrix.*Temps(1:z_ei_index,:))./((k(2:z_ei_index+1).*dz(1:z_ei_index) + k(1:z_ei_index).*dz(2:z_ei_index+1))*onesMatrix);
    boundary_Pv = 611 .* exp( (-51058/8.31).*(1./boundaryTemps - 1/273.16) ); % compute vapor pressure at ice table
    boundary_rhov = boundary_Pv .* (0.01801528 / 6.022140857e23) ./ (1.38064852e-23 .* boundaryTemps); % compute vapor densities at ice table
    boundary_rhov_mean = mean(boundary_rhov(:,:),2); % get mean for each layer's boundary's vapor densities
    
    % Calculate mean surface temperature
    Tsurf_mean = mean(Tsurf); % mean surface temperature
    
    % Convert atmospheric pressure to density at the surface
    atmosDensitySF12 = atmosPressSF12 * PtoRhoFactor / Tsurf_mean;
    
    % Check if ground ice stable - it is not if all vapor densities in
    % the lag deposit are greater than the atm vapor density
    if sum(boundary_rhov_mean < atmosDensitySF12) == 0
        % This is the case that ground ice not stable
        
        % Just set it equal to ice sheet table (aka no ground ice)
        z_pf = depthsAtLayerBoundaries(z_ei_index);
        
        z_eqConverge = -1;  % no convergence criteria needed, as it just gets set to excess ice table
        
        fprintf ('Ground ice counter = %2.0f. Saturation vapor densities in the lag deposit are all greater than the atmospheric vapor density. Ground ice not stable. z_pf = z_ei = %3.9f.\n', groundicecounter, z_pf);
        
        if triedTwice > 0
            % break out if z_eq > ice sheet table twice
            FOUNDZeq = 1;
            PFICE = 0;
            
            % No need to track which convergence criteria used in this case
            tempTracker1 = -1;
            tempTracker2 = -1;
            tempTracker3 = -1;
        else
            % Test again if just found a solution where z_eq > excess ice
            % interface to make sure the solution belongs there
            triedTwice = triedTwice + 1;
            FOUNDZeq = 0;
        end
    elseif sum(isinf(boundary_rhov_mean)) > 0
        % This is the case where numerical instability was introduced.
        % Ideally would make dt smaller but then the code will take even 
        % longer to run. For dt = 500 sec, this case can happen if pore filling
        % ice is stable within ~1mm so just going to set z_pf =
        % s.thinnestLayer and then move on since pore filling ice is stable
        % so no excess ice loss will happen in any case
        z_pf = s.thinnestLayer;
        
        z_eqConverge = -1;  % no convergence criteria needed, as it just gets set to excess ice table
        
        fprintf ('Ground ice counter = %2.0f. Saturation vapor densities in the lag deposit are -Inf. Ground ice stable very close to surface. z_pf = %3.9f.\n', groundicecounter, s.thinnestLayer);
        
        if groundicecounter > 1
            % break out if it gets here as long as its not 1st iteration on
            % solution
            FOUNDZeq = 1;
            PFICE = 1;
            
            % No need to track which convergence criteria used in this case
            tempTracker1 = -1;
            tempTracker2 = -1;
            tempTracker3 = -1;
        else
            % Test again to make sure the solution belongs there
            FOUNDZeq = 0;
        end
    else
        % This is the case where ground ice is stable
        
        % Vapor densities at middle of layers
        middle_Pv =  611 .* exp( (-51058/8.31).*(1./Temps(1:z_ei_index,:) - 1/273.16) ); % compute vapor pressure for every layer at every time
        middle_rhov = middle_Pv .* (0.01801528 / 6.022140857e23) ./ (1.38064852e-23 .* Temps(1:z_ei_index,:)); % compute vapor densities
        % middle_Pv_mean = mean(middle_Pv(:, :), 2); % annual mean Pv at every layer - unneeded
        middle_rhov_mean = mean(middle_rhov(:,:),2); % annual mean vapor density at every layer
        
        % Combines into one array to interpolate from
        alldepths = [depthsAtMiddleOfLayers(1:z_ei_index)'; depthsAtLayerBoundaries(1:z_ei_index)'];
        alldepths = reshape(alldepths, [size(depthsAtMiddleOfLayers(1:z_ei_index), 1)+size(depthsAtLayerBoundaries(1:z_ei_index), 1),1]);
        all_rhov = [middle_rhov_mean'; boundary_rhov_mean'];
        all_rhov = reshape(all_rhov, [size(middle_rhov_mean, 1)+size(boundary_rhov_mean, 1),1]);
        
        if groundicecounter > 1
            % Checks if atmospheric value is bigger than first layer's rhov,
            % because if it is, then you would just get a NaN in calculating
            % z_eq. This should be the case where we have surface ice or even
            % accumulation
            if atmosDensitySF12 > all_rhov(1)
                % Everything set to 0 because no interfaces, all ice layers
                % Depending on z_ei would either get pore filling ice all
                % the way to the surface (z_ei > 0) or would get excess ice
                % all the way up to the surface (z_ei = 0) - deal with
                % changing z_ei outside of this equilibrium depth calculation though
                z_eq = 0
                z_pf = 0;
                FOUNDZeq = 1;
                % Set these value just in case want to print as output
                z_eqConverge = -1; 
                tempTracker1 = -1;
                tempTracker2 = -1;
                tempTracker3 = -1;
                fprintf ('Ground ice counter = %2.0f: Ground ice stable up to the surface. z_eq = z_pf = %3.9f\n', groundicecounter, z_eq);
            else
                % Interpolates all rhos and depths
                z_eq = interp1(all_rhov, alldepths, atmosDensitySF12)
                
                % With v11, can handle very thin layers now
                % Keeps track of this to output in case there are hidden bugs
                z_eqConverge = abs(z_pf - z_eq);
                
                % Uses dz/2 as convergence criteria - find layer z_eq would be in (which is the
                % index from the min + 1) and use that dz/2 for convergence
                % uses original dz though (dz_E) as to avoid a constantly
                % changing target. value will be dz_E(tmpindex_pfIce+1)/2
                [tmp, tmp_z_pf_index] = min(abs(depthsAtLayerBoundaries_test - z_eq));
                
                fprintf ('Ground ice counter = %2.0f: Ground ice table testing = %3.9f, z_eq based on Temps = %3.9f, difference is = %3.9f and convergence criteria = %3.9f\n', groundicecounter, z_pf, z_eq, z_eqConverge, min([dz_test(tmp_z_pf_index+1)/2 0.003]));
                
                if z_eqConverge < min([dz_test(tmp_z_pf_index+1)/2 0.003]) && groundicecounter > 1
                    z_pf = z_eq;
                    FOUNDZeq = 1;
                    
                    % To help keep track of which convergence criteria it is
                    % using, will eventually get rid of this if I decide this
                    % method is fine and won't involve any modifications
                    tempTracker1 = dz_test(tmp_z_pf_index+1)/2;
                    tempTracker2 = min([dz_test(tmp_z_pf_index+1)/2 0.003]);
                    if(tempTracker2 == 0.003)
                        % converged based on 3 mm
                        tempTracker3 = 0;
                    else
                        % converged based on dz/2 criteria
                        tempTracker3 = 1;
                    end
                else
                    % stop if it is taking too long...it is close
                    % enough and just oscillating between solutions
                    if groundicecounter > 14
                        FOUNDZeq = 1;
                        fprintf ('Ground ice counter = %2.0f: z_pf = %3.9f, z_eq based on Temps = %3.9f, difference is = %3.9f and convergence criteria = %3.9f. Set to mean of the values: %3.9f \n', groundicecounter, z_pf, z_eq, z_eqConverge, min([dz_test(tmp_z_pf_index+1)/2 0.003]),mean([z_pf z_eq]));
                        z_pf = mean([z_pf z_eq]); 
                        tempTracker1 = -1;
                        tempTracker2 = -1;
                        tempTracker3 = -1;
                    else
                        % this is the case to keep testing....
                        FOUNDZeq = 0;
                        z_pf = z_eq;
                    end
                end
            end
        else
            % The first time through, given we haven't yet calculated
            % a new z_pf, we can't compare to the convergence criteria, and
            % therefore need to check at least twice.
            % Checks to see if the equilibrium depth would be essentially
            % above the surface, suggesting surface ice should be stable
            % If there is a lag deposit , pore-filling ice will be full all
            % the way up
            if atmosDensitySF12 > all_rhov(1)
                z_eq = 0;
                z_pf = 0;
                fprintf ('Ground ice counter = %2.0f, Ice stable up to surface.', groundicecounter);
            else
                z_eq = z_pf;
                fprintf ('Ground ice counter = %2.0f, z_eq being set to incoming z_pf value from last timestep, %3.9f\n', groundicecounter, z_pf);
            end
            
        end
        
    end
    
    groundicecounter = groundicecounter+1;
    
end