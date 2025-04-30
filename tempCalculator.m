function [Temps, Tsurf, frostMasses, flatIRSaved, flatVisSaved, iceTableTConverge] = tempCalculator13(nLayers, nStepsInYear, s, dt, k, dz, rhoc, f, Kappa, sf, visScattered, flatVis, albedo, IRdown, flatIR, sky, CO2FrostPoints, depthsAtMiddleOfLayers, sigmasb, slope, Tref, z_ei_index)
%TEMPCALCULATOR13 Runs the actual thermal model part to invert for the temperature array

% October 30, 2018: Adds in feature for polar night from Kieffer paper
% which uses downwelling IR radiation as either maximum of what we were
% using before OR that from surface CO2 frost

%% Setup numerical grids
% Get ready to run model
Temps = zeros(nLayers, nStepsInYear);
Tsurf = zeros(1,nStepsInYear);
lastTimestepTemps = zeros(nLayers, s.runTime);
oldTemps = zeros(nLayers, 1);
T_e = zeros(nLayers,1);

frostMass = 0;
frostMasses = zeros(1,nStepsInYear);

oldTemps(:) = Tref;

dt_over_rhocdz = dt./(rhoc.*dz);
alpha_u = (2*k.*circshift(k,[1 0])./(k.*circshift(dz,[1 0]) + circshift(k,[1 0]).*dz)).*dt_over_rhocdz;
alpha_u(1) = 0;
alpha_d = (2*k.*circshift(k,[-1 0])./(k.*circshift(dz,[-1 0]) + circshift(k,[-1 0]).*dz)).*dt_over_rhocdz;
alpha_d(end) = 0;

dia_e = zeros(nLayers, 1);
dia_e = 1 - (1-f)*alpha_u - (1-f)*alpha_d;
dia_e(nLayers) = 1 - (1-f)*alpha_u(end);
boundary = zeros(nLayers, 1); % column vector
boundary(end) = dt*s.Q/(rhoc(end)*dz(end));
k1_e = (1-f)*alpha_d;
k3_e = (1-f)*alpha_u;

dia_i = zeros(nLayers, 1);
dia_i = 1 + f*alpha_u + f*alpha_d;
dia_i(nLayers) = 1+f*alpha_u(end);
Amatrix_i = zeros(nLayers, nLayers);
Amatrix_i = diag(dia_i) + diag(-f*alpha_u(2:end),-1) + diag(-f*alpha_d(1:end-1),1);
Amatrix_i = sparse(Amatrix_i);

beta = k(1)*dt/(rhoc(1)*dz(1)*dz(1));

% Total fluxes
%{
Fin = (sf + visScattered.*sky + flatVis.*(1-sky)).*(1-albedo) + (IRdown.*sky + flatIR.*(1-sky)).*s.emis;
Fin_frost = (sf + visScattered.*sky + flatVis.*(1-sky)).*(1-s.albedoFrost) + (IRdown.*sky + flatIR.*(1-sky)).*s.emisFrost;
Fin_i = (circshift(sf, [-1 0]) + circshift(visScattered, [-1 0]).*sky + circshift(flatVis, [-1 0]).*(1-sky)).*(1-albedo) + (circshift(IRdown, [-1 0]).*sky + circshift(flatIR, [-1 0]).*(1-sky)).*s.emis;
%}
% From Keiffer but not in mid-latitude papers. Essentially adds in some
% downwelling IR even during polar night.
DownIRPolarNight = s.downwellingPerc.*sigmasb.*CO2FrostPoints.^4;
maxIRDown = max(IRdown.*sky, DownIRPolarNight);
% Calculate total incoming flux vectors
Fin = (sf + visScattered.*sky + flatVis.*(1-sky)).*(1-albedo) + (IRdown.*sky + flatIR.*(1-sky)).*s.emis;
Fin_frost = (sf + visScattered.*sky + flatVis.*(1-sky)).*(1-s.albedoFrost) + (maxIRDown + flatIR.*(1-sky)).*s.emisFrost;
Fin_i = (circshift(sf, [-1 0]) + circshift(visScattered, [-1 0]).*sky + circshift(flatVis, [-1 0]).*(1-sky)).*(1-albedo) + (circshift(IRdown, [-1 0]).*sky + circshift(flatIR, [-1 0]).*(1-sky)).*s.emis;



% For frost mass calculations
timestepSkinDepth = sqrt(Kappa(1)*dt/pi);
gamma_frost = -(1/s.L_CO2)*(2*k(1)*dt/dz(1));
theta_frost = (dt/s.L_CO2).*(2.*k(1).*CO2FrostPoints/dz(1) - Fin_frost + s.emisFrost.*sigmasb.*CO2FrostPoints.^4);
theta_frost_i = circshift(theta_frost, [-1 0]);
defrosting_decrease = exp(-depthsAtMiddleOfLayers./timestepSkinDepth);

% For making temp calculations in the time loop faster by pulling some
% calculations out of the loop
dzover2k_topLayer = dz(1)/(2*k(1));
ab_denom = 4*s.emis*sigmasb*dzover2k_topLayer;
emissig3 = 3*s.emis*sigmasb;

% Use these if expecting to run into the CO2 frost calculation a lot,
% otherwise it's actually just slower to do them every orbital solution
%top_diae_frosted = 1 - (1-f)*(alpha_d(1)+2*beta);
%top_A_frosted = 1 + f*(alpha_d(1)+2*beta);

% Calculate a and b's for surface temperature calculation
aa = dzover2k_topLayer*(Fin(1) + emissig3*Tref^4)/(1+(ab_denom*Tref^3));
b = 1/(1+(ab_denom*Tref^3));
Tsurf(1) = aa+b*Tref;

%% Loop through time calculating temperatures
for yr=1:s.runTime
    
    for n = 1:nStepsInYear
        
        Tfrost = CO2FrostPoints(n);
        
        if frostMass == 0
            
            % To make faster
            T4 = emissig3*Tref^4;
            T3 = 1+(ab_denom*Tref^3);
            
            % Calculations for new temperatures
            b = 1/(1+(ab_denom*Tref^3));
            a_e = dzover2k_topLayer*(Fin(n) + T4)/T3;
            a_i = dzover2k_topLayer*(Fin_i(n) + T4)/T3;
            boundary(1) = 2*beta*((1-f)*a_e + f*a_i);
            
            % Also for speed
            alpha_b2beta = alpha_d(1)+(2-2*b)*beta;
            
            % Explicit Part
            dia_e(1) = 1 - (1-f)*alpha_b2beta;
            T_e = k3_e.*circshift(oldTemps,[1 0]) + k1_e.*circshift(oldTemps,[-1 0]) + dia_e.*oldTemps + boundary;
            
            % Implicit Part
            Amatrix_i(1,1) = 1 + f*alpha_b2beta;
            Temps(:,n) = Amatrix_i\T_e;
            
            Tsurf(n) = a_i + b*Temps(1,n);  % Uses implicit a with new T calculation- instantanous balance
            
            frostMass = 0;
            
            if Tsurf(n) < Tfrost
                
                deltaTsurf = Tfrost - Tsurf(n);
                frostMass = deltaTsurf*rhoc(1)*timestepSkinDepth/s.L_CO2;
                Temps(:,n) = Temps(:,n) + deltaTsurf.*defrosting_decrease;
                Tsurf(n) = Tfrost;
                
            end
            
        elseif frostMass > 0
            
            % In this case, a = Tfrost and b = 0
            boundary(1) = 2*beta*Tfrost;
            
            % Explicit Part
            dia_e(1) = 1 - (1-f)*(alpha_d(1)+2*beta);
            T_e = k3_e.*circshift(oldTemps,[1 0]) + k1_e.*circshift(oldTemps,[-1 0]) + dia_e.*oldTemps + boundary;
            
            % Implicit Part
            Amatrix_i(1,1) =  1 + f*(alpha_d(1)+2*beta);
            Temps(:,n) = Amatrix_i\T_e;
            
            Tsurf(n) = Tfrost;
            
            % Semi-implicit frost mass calculation
            frostMass = frostMass + (1-f)*(gamma_frost*oldTemps(1) + theta_frost(n)) + f*(gamma_frost*Temps(1,n) + theta_frost_i(n));
            
            if frostMass < 0
                
                shiftedFrostMasses = circshift(frostMasses,[1 0]);
                timeDefrosted = sqrt((0-frostMass)/(shiftedFrostMasses(n)-frostMass));
                deltaTsurf2 = -frostMass*s.L_CO2/(rhoc(1)*timestepSkinDepth*timeDefrosted);
                Tsurf(n) = Tfrost + deltaTsurf2;
                Temps(:,n) = Temps(:,n) + deltaTsurf2.*defrosting_decrease;
                frostMass = 0;
                
            end
        else
            fprintf('Frost mass is negative?! Umm...Houston, we have a problem!');
        end
        
        oldTemps(:) = Temps(:,n);
        Tref = Tsurf(n);
        frostMasses(n) = frostMass;
        
    end
    
    lastTimestepTemps(:,yr) = Temps(:,n);  % To compare for convergence
    
    %fprintf('You''re %2.0f / %2.0f of the way there!\n ', yr, s.runTime);
    
    % If at the end of the windup time, calculate the annual average surface
    % temperature and set all temperatures to that for the real run.
    if yr == s.windupTime
        windupTemp = mean(Tsurf);
        oldTemps(:) = windupTemp;
        %fprintf('Windup done! Setting all temperatures equal to %4.2f K.\n', windupTemp);
    end
    
    % Print out the difference between the temperatures of the last timestep
    % and the year before, and compare to a convergence criteria, though
    % currently it doesn't actually run until reaching this criteria
    if yr == s.runTime
        tempDiffs = lastTimestepTemps(:,s.runTime) - lastTimestepTemps(:,s.runTime-1);
        whichConverge = abs(tempDiffs) < s.convergeT;
        [maxDiffVal, maxDiffIndex] = max(abs(tempDiffs));
        if sum(whichConverge) == size(whichConverge,1)
            fprintf('It converged to within %3.7f K! Congrats!\n', maxDiffVal);
        else
            fprintf('It did not converge to within %4.3f K. You only got to within %3.7f K at layer index of %3.0f corresponding to a depth of %5.3f m.\n', s.convergeT, maxDiffVal, maxDiffIndex, depthsAtMiddleOfLayers(maxDiffIndex));
        end
        %% Print convergence of last timestep at ice table (or 1st layer if no ice table)
        if z_ei_index > 0
            iceTableTConverge = abs(tempDiffs(z_ei_index));
            fprintf('The last timestep at the ice table converged to %3.7f K.\n', iceTableTConverge);
        else
            iceTableTConverge = abs(tempDiffs(1));
            fprintf('The last timestep at the 1st layer converged to %3.7f K.\n', iceTableTConverge);
        end
    end
    
end

%% Save if going to use this output for reradiation

    
if slope == 0 && s.slope_rerad == 1
    
    % Important that these are this size, and not the reverse
    flatIRSaved = zeros(nStepsInYear,1);
    flatVisSaved = zeros(nStepsInYear,1);
    
    for n = 2:nStepsInYear
        if frostMasses(n) > 0
            flatIRSaved(n) = s.emisFrost * sigmasb * Tsurf(n-1)^4;
            flatVisSaved(n) = s.albedoFrost * sf(n);
        else
            flatIRSaved(n) = s.emis * sigmasb * Tsurf(n-1)^4;
            flatVisSaved(n) = albedo * sf(n);
        end
        
    end
    
    if frostMasses(1) > 0
        flatIRSaved(1) = s.emisFrost * sigmasb * Tsurf(end)^4;
        flatVisSaved(1) = s.albedoFrost * sf(1);
    else
        flatIRSaved(1) = s.emis * sigmasb * Tsurf(end)^4;
        flatVisSaved(1) = albedo * sf(1);
    end
    
else
    flatIRSaved = flatIR;
    flatVisSaved = flatVis;
end
end

