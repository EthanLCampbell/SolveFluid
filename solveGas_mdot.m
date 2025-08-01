function [mdot_lbm] = solveGas_mdot(P_up,T_up,P_d,CdA,fluid,FlowType)
%function [mdot_lbm] = solveGas_mdot(14.7,70,250,.99,oxygen,1);
% Written by: Kevin Dille 04/24/2021
% Function determines the mass flow rate of a gas accounting for real gas
% effects for sonic and subsonic flow.
% 
% Modified by: Kevin Dille 03/14/2022
% Function can now handle both injectors and venturis with the "FlowType"
% variable. Venturi calculations allow for replacing the "Sonicmdot"
% function
% 
% INPUTS:
%   P_up:       Upstream total pressure of the fluid    (psi)
%   T_up:       Upstream total temperature of the fluid (F)
%   P_d:        Downstream pressure of the fluid        (psi)
%   CdA:        Total flow area                         (in^2)
%   Fluid:      Fluid name for use with refprop
%   FlowType:   Type of flow measurement, with:
%                   1: Venturi - solves CdA assuming sonic flow if
%                   downstream pressure is <80% the upstream pressure.
%                   2: Injector - solves CdA for sonic or subsonic flow.
%                   Sonic flow is used with the downstream pressure is less
%                   than the critical pressure of the upstream conditions.
%
% OUTPUT:
%   mdot_lbm:   Mass flow rate                          (lbm/s)

%% Convert to Refprop Units
P_up_rp = P_up*6.895; %[kPa]
T_up_rp = (T_up+459.67)*5/9; %[K]
P_d_rp = P_d*6.895; %[kPa]
CdA_m = CdA/(39.37^2); %[m^2]

%% Determine if flow is sonic or subsonic and compute the mass flow rate
[cstar,Pcrit,~] = refpropm('~','T',T_up_rp,'P',P_up_rp,fluid); %critical flow factor, critical Pressure at throat
Mw = refpropm('M','T',T_up_rp,'P',P_up_rp,fluid); %molecular weight of fluid
R = 8314.5/Mw; %Gas constant

if FlowType == 1 % Venturi Mass Flow rate
    if P_d_rp < (0.8*P_up_rp) % Venturi is choked
        mdot_kg = CdA_m*cstar*(P_up_rp*1000)/sqrt(R*T_up_rp); %[kg/s] critical flow funciton
    else % Venturi is unchoked
        mdot_kg = 0;
    end
    
elseif FlowType == 2 % Injector Mass Flow rate
    if P_d_rp < Pcrit
        %Flow is choked
        mdot_kg = CdA_m*cstar*(P_up_rp*1000)/sqrt(R*T_up_rp); %[kg/s] critical flow funciton
    else %Flow is not choked
        [s_0,h_0] = refpropm('SH','T',T_up_rp,'P',P_up_rp,fluid); %Determine entropy and enthalpy for upstream conditions (assumed total properties)
        [h_inj,rho_inj] = refpropm('HD','P',P_d_rp,'S',s_0,fluid); %Enthalpy and density at injection pressure assuming flow is isentropic
        vel_inj = sqrt(2*(h_0-h_inj)); %[m/s] Isentropic velocity
        mdot_kg = rho_inj*vel_inj*CdA_m; %[kg/s]
    end
end

mdot_lbm = mdot_kg*2.20462; %[lbm/s]
end