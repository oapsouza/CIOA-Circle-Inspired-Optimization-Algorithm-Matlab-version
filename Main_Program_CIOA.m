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

close all
clear all
clc
format long
tic


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  INPUT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nvar = 6;                   % Number of variables of the objective function
Nag = 250;                  % Number of search agents
Nit = 400;                  % Number of iterations
Ub = [1 1 1 1 16 16];       % Upper bound of the design variables
Lb = [0 0 0 0 1e-5 1e-5];   % Lower bound of the design variables

ThetaAngle = 17;            % Theta angle in degrees (suggestions: 13, 17 or 19°)
GlobIt = 0.85;              % Proportion of iterations for global search (sugg.: 0.75 to 0.95)


%%%%%%%%%%%%%%%%%%%%%%%  ALGORITHM  INITIALIZATION  %%%%%%%%%%%%%%%%%%%%%%%

Angle = (ThetaAngle*pi/180);

Convergence = zeros(1,Nit);
BestVar = zeros(Nit,Nvar);              % Initializing solution vectors
x = ones(1,Nvar);
Bestobj = zeros(Nit,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%  RADIUS CALCULATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%

raux = 1*(sqrt(max(Ub)-min(Lb)))/Nag;   % Radius auxiliary variable
r = zeros(1,Nag);                       % Initializing radius vector

h = Nag;
while h >= 1;
    r(h) = raux*(h)^2/Nag;              % Completing radius vector
    h = h-1;
end


%%%%%%%%%%%%%%%%%%%%%  RANDOM SOLUTION DETERMINATION  %%%%%%%%%%%%%%%%%%%%%

i0 = 1;
Var0 = zeros(Nag,Nvar);
obj0 = zeros(Nag,1);
while i0 <= Nag;
    Var0(i0,:) = Lb+(Ub-Lb).*rand(1,Nvar); % Random variables generation
    i0 = i0+1;
end

Sol0 = zeros(1,Nag);
i1 = 1;
while i1 <= Nag;
    x = Var0(i1,:);                     % Assigns the random variables to the design variables
    [sol] = Objective_Function_CIOA(x); % Calls the objective function of the auxiliary file
    Sol0(i1) = sol;                     % Random solutions for each search agent 
    i1 = i1+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%  GLOBAL SEARCH LOOP  %%%%%%%%%%%%%%%%%%%%%%%%%%%

It = 1;
while It <= floor(Nit*GlobIt);
    [class,ord] = sort(Sol0);           % Solution vector classification
    i2 = 1;
    while i2 <= Nag;
        i3 = 1;
        while i3 <= Nvar; % Update of the position of agents in the global search
            if mod(i3,2)==0;
                Var0(i2,i3) = (Var0(ord(i2),i3))-((0.0+(1-0.0)*rand())*r(i2)*sin((It*Angle)-Angle))+((0.0+(1-0.0)*rand())*r(i2)*sin(It*Angle));
            else
                Var0(i2,i3) = (Var0(ord(i2),i3))-((0.0+(1-0.0)*rand())*r(i2)*cos((It*Angle)-Angle))+((0.0+(1-0.0)*rand())*r(i2)*cos(It*Angle));
            end
            i3 = i3+1;
        end
        i2 = i2+1;
    end

    i4 = 1;
    while i4 <= Nag;
        i5 = 1;
        while i5 <= Nvar; % Update of variables that exceeded the limits
            if Var0(1,i5) > Ub(i5);
                Var0(1,i5) = Ub(i5);
            else if Var0(1,i5) < Lb(i5);
                    Var0(1,i5) = Lb(i5);
                end
            end
            if Var0(i4,i5) > Ub(i5);
                Var0(i4,i5) = Var0(1,i5);
            else if Var0(i4,i5) < Lb(i5);
                    Var0(i4,i5) = Var0(1,i5);
                end
            end
            i5 = i5+1;
        end
        i4 = i4+1;
    end

    i6 = 1;
    while i6 <= Nag;
        x = Var0(i6,:);                          % Assigns the new particle positions to the design variables
        [sol, obj] = Objective_Function_CIOA(x); % Calls the objective function of the auxiliary file
        Sol0(i6) = sol;                          % Solutions according to the new positions for each search agent
        obj0(i6) = obj;
        i6 = i6+1;
    end

    [Best,Position] = min(Sol0);        % Finding the best solution among agents
    BestSolution = Best;                % Saves the best solution of the iteration
    BestVar(It,:) = Var0(Position,:);   % Variables corresponding to the best solution
    Bestobj(It) = obj0(Position);
    Convergence(It) = BestSolution;

    if mod(It,floor(360/ThetaAngle))==0;
        r = r.*0.99;                    % Radius reduction every 360º (full lap) of the algorithm
    end
    It = It+1;
end


%%%%%%%%%%%%%%%%%%%%%%  LOCAL SEARCH INITIALIZATION  %%%%%%%%%%%%%%%%%%%%%%

Lb1 = BestVar(It-1,:)-(ones(1,Nvar).*((Ub-Lb)/10000)); % Reducing the limits of the design variables
Ub1 = BestVar(It-1,:)+(ones(1,Nvar).*((Ub-Lb)/10000));

i7 = 1;
Var1 = zeros(Nag,Nvar);
while i7 <= Nag;
    Var1(i7,:) = BestVar(It-1,:);    % All search agents will begin their search in the promising region
    i7 = i7+1;
end

x = Var1(1,:);
[sol] = Objective_Function_CIOA(x);  % Calls the objective function of the auxiliary file
Sol1 = zeros(1, Nag);

i8 = 1;
while i8 <= Nag;
    Sol1(i8) = sol;
    i8 = i8+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%  LOCAL SEARCH LOOP  %%%%%%%%%%%%%%%%%%%%%%%%%%%

while It <= Nit;
    [class,ord] = sort(Sol1);           % Solution vector classification
    i9 = 1;
    while i9 <= Nag;
        i10 = 1;
        while i10 <= Nvar;
            if mod(i10,2)==0;
                Var1(i9,i10) = (Var1(ord(i9),i10))-((0.0+(1-0.0)*rand())*r(i9)*sin((It*Angle)-Angle))+((0.0+(1-0.0)*rand())*r(i9)*sin(It*Angle));
            else
                Var1(i9,i10) = (Var1(ord(i9),i10))-((0.0+(1-0.0)*rand())*r(i9)*cos((It*Angle)-Angle))+((0.0+(1-0.0)*rand())*r(i9)*cos(It*Angle));
            end
            i10 = i10+1;
        end
        i9 = i9+1;
    end

    i11 = 1;
    while i11 <= Nag;
        i12 = 1;
        while i12 <= Nvar;
            if Var1(1,i12) > Ub(i12);
                Var1(1,i12) = Ub(i12);
            else if Var1(1,i12) < Lb(i12);
                    Var1(1,i12) = Lb(i12);
                end
            end
            if Var1(i11,i12) > Ub1;
                Var1(i11,i12) = Var1(1,i12);
            else if Var1(i11,i12) > Ub(i12);
                    Var1(i11,i12) = Ub(i12);
                end
            end
            if Var1(i11,i12) < Lb1;
                Var1(i11,i12) = Var1(1,i12);
            else if Var1(i11,i12) < Lb(i12);
                    Var1(i11,i12) = Lb(i12);
                end
            end            
            i12 = i12+1;
        end        
        i11 = i11+1;
    end

    i13 = 1;
    while i13 <= Nag;
        x = Var1(i13,:);
        [sol, obj] = Objective_Function_CIOA(x); % Calls the objective function of the auxiliary file
        Sol1(i13) = sol;
        obj0(i13) = obj;
        i13 = i13+1;
    end

    [Best,Position] = min(Sol1);
    BestSolution = Best;
    BestVar(It,:) = Var1(Position,:);
    Bestobj(It) = obj0(Position);
    Convergence(It) = BestSolution;

    if mod(It,floor(360/ThetaAngle))==0;
        r = r.*0.99;
    end
    It = It+1;
end


%%%%%%%%%%  CONVERGENCE CURVE CONSTRUCTION AND RESULTS DISPLAY  %%%%%%%%%%%

ConvergenceCurve = zeros(1,Nit);
ConvergencePlot = zeros(1,Nit);
ConvergenceCurve(1) = Convergence(1);
ConvergencePlot(1) = Bestobj(1);

i14 = 2;
while i14 <= Nit;
    if Convergence(i14) < ConvergenceCurve(i14-1);
        ConvergenceCurve(i14) = Convergence(i14);
        ConvergencePlot(i14) = Bestobj(i14);
    else
        ConvergenceCurve(i14) = ConvergenceCurve(i14-1);
        ConvergencePlot(i14) = ConvergencePlot(i14-1);
    end
    i14 = i14+1;
end

plot(ConvergencePlot','b-','LineWidth',2)  % Convergence curve plot
xlabel('\fontsize{11}\bf Iteration Number')
ylabel('\fontsize{11}\bf Objective Function Value')
grid on

[BestSol,PositionSol] = min(Convergence);
BestSolution = ConvergencePlot(Nit)        % Display of best solution
BestVariables = BestVar(PositionSol,:)     % Display of the design variables that constitute the best solution


%%%%%%%%%%%%%%%%%%%%%%%%%  CONSTRAINT VIOLATIONS  %%%%%%%%%%%%%%%%%%%%%%%%%

x = BestVariables;
[sol, obj, rg, rh, nug, nuh] = Objective_Function_CIOA(x);

Violg = 0;
for icong=1:nug;                        % Violation of inequality constraints
    Violg = Violg+max(rg(icong),0);
end

Violh = 0;
for iconh=1:nuh;                        % Violation of equality constraints
    Violh = Violh+(max(abs(rh(iconh))-0.0001,0));
end

if (nug+nuh)==0;
    Violation = (Violg+Violh)/1;        % Total violation
else
    Violation = (Violg+Violh)/(nug+nuh);
end

Violation                               % Display violation

toc

save Results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CITE AS:
% de Souza, O.A.P., and Miguel, L.F.F. CIOA: Circle-Inspired Optimization
% Algorithm, an algorithm for engineering optimization. SoftwareX 2022;19:101192
% https://doi.org/10.1016/j.softx.2022.101192
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
