function [mdot_lbm] = solveLiq_mdot(P_up,T_up,P_d,CdA,fluid,FlowType)
% Written by: Kevin Dille 04/24/2021
% Function determines the mass flow rate of a gas accounting for real gas
% effects for sonic and subsonic flow.
% 
% Modified by: Kevin Dille 03/14/2022
% Function can now solve mass flow through an injector or cavitating
% venturi through the use of the "FlowType" variable
% 
% INPUTS:
%   P_up:       Upstream total pressure of the fluid    (psi)
%   T_up:       Upstream total temperature of the fluid (F)
%   P_d:        Downstream pressure of the fluid        (psi)
%   CdA:        Total flow area                         (in^2)
%   Fluid:      Fluid name for use with refprop
%   FlowType:   Type of flow measurement (1 or 2), with:
%                   1: Venturi - solves CdA assuming cavitating flow if
%                   downstream pressure is <80% the upstream pressure.
%                   2: Injector - solves CdA for noncavitating flow.
%                   Bernoulli velocity is taken from the dP across the 
%                   element.
%
% OUTPUT:
%   mdot_lbm:   Mass flow rate                          (lbm/s)


%% Convert to Refprop Units
P_up_rp = P_up*6.895; %[kPa]
T_up_rp = (T_up+459.67)*5/9; %[K]
P_d_rp = P_d*6.895; %[kPa]
CdA_m = CdA/(39.37^2); %[m^2]

%% Determine Properties
rho_m = refpropm('D','T',T_up_rp,'P',P_up_rp,fluid); %[kg/m^3]
if FlowType == 1 % Cavitating Venturi
    if P_d_rp < (0.8*P_up_rp) % Flow is cavitating
        P_vap_rp = refpropm('P','T',T_up_rp,'Q',1,fluid); %[kPa]
        vel_inj = sqrt(2*1000*(P_up_rp-P_vap_rp)/rho_m); %[m/s]
    else % Venturi is not cavitating
        vel_inj = 0; %[m/s]
    end
elseif FlowType == 2 % Injector
    vel_inj = sqrt(2*(P_up_rp-P_d_rp)*1000/rho_m); %[m/s]
end

%% Compute Mass Flow
mdot_kg = rho_m*vel_inj*CdA_m; %[kg/s]

mdot_lbm = mdot_kg*2.20462; %[lbm/s]
end