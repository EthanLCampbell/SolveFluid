function [P_up] = solveLiq_CdA(CdA,T_up,P_d,mdot_lbm,fluid,FlowType)
% Written by: Ethan Labianca-Campbell
% Function determines the upstream pressure required for a fluid to achieve
% a provided mass flow rate through a given area (CdA).
% 
% Modified by: Ethan Labianca-Campbell 08/01/2025
% Function now differentiates between cavitating venturis and dP injectors
% with the "FlowType" variable to be able to solve CdA for both injectors 
% and venturis. Venturi assumes 80% pressure recovery.
% 
% INPUTS:
%   P_up:       Upstream total pressure of the fluid    (psi)
%   T_up:       Upstream total temperature of the fluid (F)
%   P_d:        Downstream pressure of the fluid        (psi)
%   mdot_lbm:   Mass flow rate                          (lbm/s)
%   fluid:      Fluid name for use with refprop
%   FlowType:   Type of flow measurement (1 or 2), with:
%                   1: Venturi - solves CdA assuming cavitating flow if
%                   downstream pressure is <80% the upstream pressure.
%                   2: Injector - solves CdA for noncavitating flow.
%                   Bernoulli velocity is taken from the dP across the 
%                   element.
%
% OUTPUT:
%   CdA:    Total flow area                             (in^2)

%% Convert to Refprop Units
P_up_rp = P_up*6.895; %[kPa]
T_up_rp = (T_up+459.67)*5/9; %[K]
P_d_rp = P_d*6.895; %[kPa]
mdot_kg = mdot/2.20462; %[kg/s]

%% Determine Bernoulli Velocity
rho_inj = refpropm('D','T',T_up_rp,'P',P_up_rp,fluid); %[kg/m^3]
if FlowType == 1 % Venturi (cav flow)
    if P_d_rp < (0.8*P_up_rp)
        % Flow is cavitating
        P_vap_rp = refpropm('P','T',T_up_rp,'Q',1,fluid); %[kPa]
        vel_inj = sqrt(2*1000*(P_up_rp-P_vap_rp)/rho_inj); %[m/s]
    else
        % Flow is not cavitating
        vel_inj = 0;
    end
    
elseif FlowType == 2 % Injector
    vel_inj = sqrt(2*1000*(P_up_rp-P_d_rp)/rho_inj); %[m/s]
end
CdA_m = mdot_kg/(rho_inj*vel_inj); %[m^2]
CdA = CdA_m*(39.37^2); %[in^2]


end
