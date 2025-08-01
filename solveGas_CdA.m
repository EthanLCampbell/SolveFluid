function [CdA] = solveGas_CdA(P_up,T_up,P_d,mdot,fluid,FlowType)
% Written by: Ethan Labianca-Campbell
% Function determines the area (CdA including discharge coefficient
% effects) of a fluid for a given mass flow rate and upstream conditions.
% 
% Modified by: Ethan Labianca-Campbell 08/01/2025
% Function works for both injectors and venturis now.
%
% INPUTS:
%   P_up:       Upstream total pressure of the fluid    (psi)
%   T_up:       Upstream total temperature of the fluid (F)
%   P_d:        Downstream pressure of the fluid        (psi)
%   mdot_lbm:   Mass flow rate                          (lbm/s)
%   Fluid:      Fluid name for use with refprop
%   FlowType:   Type of flow measurement, with:
%                   1: Venturi - solves CdA assuming sonic flow if
%                   downstream pressure is <80% the upstream pressure.
%                   2: Injector - solves CdA for sonic or subsonic flow.
%                   Sonic flow is used with the downstream pressure is less
%                   than the critical pressure of the upstream conditions.
%
% OUTPUT:
%   CdA:        Total flow area                         (in^2)

%% Convert to Refprop Units
P_up_rp = P_up*6.895; %[kPa]
T_up_rp = (T_up+459.67)*5/9; %[K]
P_d_rp = P_d*6.895; %[kPa]
mdot_kg = mdot/2.20462; %[kg/s]

%% Determine if flow is sonic or subsonic and compute the mass flow rate
[cstar,Pcrit,~] = refpropm('~','T',T_up_rp,'P',P_up_rp,fluid); %critical flow factor, P and T at throat
Mw = refpropm('M','T',T_up_rp,'P',P_up_rp,fluid); %molecular weight of fluid
R = 8314.5/Mw; %Gas constant

if FlowType == 1 % Venturi CdA Calculation
    if P_d_rp < (0.8*P_up_rp) %FUNCTION ASSUMES 80% PR IS REQUIRED FOR SONIC FLOW
        %Venturi is sonic
        CdA_m = mdot_kg*sqrt(R*T_up_rp)/(cstar*(P_up_rp*1000)); %[m^2] critical flow function solved for area
    else
        % Venturi is unchoked
        CdA_m = 0;
    end
    
elseif FlowType == 2 % Injector/Oriface CdA Calculation
    if P_d_rp < Pcrit
        %Flow is choked at outlet
        CdA_m = mdot_kg*sqrt(R*T_up_rp)/(cstar*(P_up_rp*1000)); %[m^2] critical flow function solved for area
        
    else %Flow is not choked
        [s_0,h_0] = refpropm('SH','T',T_up_rp,'P',P_up_rp,fluid); %Determine entropy and enthalpy for upstream conditions (assumed total properties)
        [h_inj,rho_inj] = refpropm('HD','P',P_d_rp,'S',s_0,fluid); %Enthalpy and density at injection pressure assuming flow is isentropic
        vel_inj = sqrt(2*(h_0-h_inj)); %[m/s] Isentropic velocity
        CdA_m = mdot_kg/(rho_inj*vel_inj); %[m^2]
    end
end

CdA = CdA_m*(39.37^2); %[in^2]



end
