function [P_up] = solveLiq_Pup(CdA,T_up,P_d,mdot_lbm,fluid,FlowType)
% Written by: Kevin Dille 04/24/2021
% Edited by: Charlie Black 7/24/2022
% Function determines the upstream pressure required for a fluid to achieve
% a provided mass flow rate through a given area (CdA).
% 
% Modified by: Kevin Dille 03/14/2022
% Function can now handle both venturis and injector-like flow with the
% "FlowType" variable.
% 
% INPUTS:
%   CdA:        Total flow area                         (in^2)
%   T_up:       Upstream total temperature of the fluid (F)
%   P_d:        Downstream pressure of the fluid        (psi)
%   mdot_lbm:   Mass flow rate                          (lbm/s)
%   Fluid:      Fluid name for use with refprop
%   FlowType:   Type of flow measurement (1 or 2), with:
%                   1: Venturi - solves CdA assuming cavitating flow if
%                   downstream pressure is <80% the upstream pressure.
%                   2: Injector - solves CdA for noncaviating flow.
%                   Velocity is taken as the Bernoulli velocity given the
%                   dP across the element.
%
% OUTPUT:
%   P_up:       Upstream total pressure of the fluid    (psi)


%% Iterate on Pressure for Mass Flow Convergance

% Check to ensure venturi size works to begin with
if FlowType == 1 % Only do this for venturis
    P_up_check = P_d/0.7999; % Calculate minimum upstream pressure to remain critical
    Min_possible_mdot = solveLiq_mdot(P_up_check,T_up,P_d,CdA,fluid,FlowType); % Calculate lowest mass flow the critical venturi can provide
    if Min_possible_mdot > mdot_lbm
        Max_possible_CdA = solveLiq_CdA(P_up_check,T_up,P_d,mdot_lbm,fluid,FlowType); % Obtain max possible flow area to provide critical metering of the mass flow requested
        warning('Your selected CdA is not compatible with the requested mdot in order to remain critical at 80% pressure recovery')
        warning('You need to decrease CdA to %.5f in2 or increase mass flow to %.3f lbm/s',Max_possible_CdA,Min_possible_mdot)
        P_up = 0;
        return
    else
    end
end

mdot_tol = 1e-4; %Percent tolerance on mass flow convergnace
mdot_dif = 2*mdot_tol; %Initial difference for entering the loop

P_up_min = P_d; %[psi] initial minimum guess for pressure using bisecting method
P_up_max = 10*P_d; %[psi] initial maximum guess for pressure using biscecting method
P_max_i = P_up_max; %[psi] Loop maximum guess for pressure to determine if this value must increase

while abs(mdot_dif) > mdot_tol
    P_up_guess = (P_up_min+P_up_max)/2; %Guessed value for P_up
    mdot_guess = solveLiq_mdot(P_up_guess,T_up,P_d,CdA,fluid,FlowType); %[lbm/s]
    
    mdot_dif = 100*(mdot_lbm-mdot_guess)/mdot_lbm; %Percent difference in mass flow rates
    if abs(mdot_dif) > mdot_tol
        if P_up_min > (P_max_i*0.99) % Catch if P_max is greater than current maximum
            P_up_max = P_max_i*10;
            P_max_i = P_up_max;
        elseif mdot_dif > 0
            P_up_min = P_up_guess;
        else
            P_up_max = P_up_guess;
        end
    end
end

P_up = P_up_guess; %[psi]
end