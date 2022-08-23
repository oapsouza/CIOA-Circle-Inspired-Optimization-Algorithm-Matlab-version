%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%            CIRCLE-INSPIRED OPTIMIZATION ALGORITHM (CIOA)            %%%
%%%            =============================================            %%%
%%%                                                                     %%%
%%%  Implemented in MATLAB R2015a(8.5.0)                                %%%
%%%                                                                     %%%
%%%  Developed by: Otávio Augusto Peter de Souza¹                       %%%
%%%                Letícia Fleck Fadel Miguel²                          %%%
%%%                                                                     %%%
%%%  ¹otavio.souza@ufrgs.br, ²letffm@ufrgs.br                           %%%
%%%                                                                     %%%
%%%  Federal University of Rio Grande do Sul (UFRGS)                    %%%
%%%  Postgraduate Program in Civil Engineering (PPGEC)                  %%%
%%%                                                                     %%%
%%%  https://doi.org/10.1016/j.softx.2022.101192                        %%%
%%%                                                                     %%%
%%%  Paper Published in: SoftwareX                                      %%%
%%%                      de Souza, O.A.P., and Miguel, L.F.F. (2022).   %%%
%%%                      CIOA: Circle-Inspired Optimization Algorithm,  %%%
%%%                      an Algorithm for Engineering Optimization.     %%%
%%%                      SoftwareX, Volume 19, July 2022, 101192        %%%
%%%                      https://doi.org/10.1016/j.softx.2022.101192    %%%
%%%                                                                     %%%
%%%  Version: V01, Updated in: January 2022                             %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sol, obj, rg, rh, nug, nuh] = Objective_Function_CIOA(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%  OBJECTIVE FUNCTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%
k1 = 0.09755988;
k2 = 0.99*k1;
k3 = 0.0391908;
k4 = 0.9*k3;

obj = -x(4);       % Objective Function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CONSTRAINTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nug = 1;           % Number of inequality constraints
nuh = 4;           % Number of equality constraints

rg = zeros(nug,1);
rh = zeros(nuh,1);

% Inequality Constraints:
rg(1) = x(5)^0.5+x(6)^0.5-4;

% Equality Constraints:
rh(1) = x(1)+k1*x(2)*x(5)-1;
rh(2) = x(2)-x(1)+k2*x(2)*x(6);
rh(3) = x(3)+x(1)+k3*x(3)*x(5)-1;
rh(4) = x(4)-x(3)+x(2)-x(1)+k4*x(4)*x(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PENALTIES  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Weight = 10^5;     % Define a value multiplied by the penalty

Pen = 0;
for icong=1:nug;
    Pen = Pen+(max(rg(icong),0));              % Penalty for inequality constraints
end

for iconh = 1:nuh;
    Pen = Pen+(max(abs(rh(iconh))-0.0001,0));  % Penalty for equality constraints
end
Pen = Weight*Pen;

sol = obj + Pen;   % Objective function penalized

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CITE AS:
% de Souza, O.A.P., and Miguel, L.F.F. CIOA: Circle-Inspired Optimization
% Algorithm, an algorithm for engineering optimization. SoftwareX 2022;19:101192
% https://doi.org/10.1016/j.softx.2022.101192
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
