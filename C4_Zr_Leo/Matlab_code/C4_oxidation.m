%% C4_oxidation
%%Main code computing all the oxidation-related parameters (oxide thickness, oxygen concentration) using finite difference

%%Author : Leo Borrel
%%Email : borrel@wisc.edu

%%Last updated : 09/03/2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Creation of the matrices used in the code

% Create the mesh in the metal and in the oxide
x = zeros(K + 2, N + 1);                % Matrix containing the N+1 (time) discretization of the space of the metal [cm]
x_ox = zeros(K_ox + 1, N + 1);          % Matrix containing the N+1 (time) discretization of the space in the oxide [cm]

% Create the concentration matrices in the metal (1 column = 1 time, 1 line = 1 position)
Co = Co_nat * ones(K + 2, N + 1);       % Matrix containing the N+1 (time) oxygen concentration at the K+2 (space) position [/cm^3]
xo = xo_nat * ones(K + 2, N + 1);       % Matrix containing the N+1 (time) oxygen fraction at the K+2 (space) position [/cm^3]
Ch = Ch_nat * ones(1, N + 1);           % Matrix containing the N+1 (time) hydrogen concentration at the K+2 (space) position [wt ppm]

% Create the concentration matrices in the oxide (1 column = 1 time, 1 line = 1 position)
Cv_ox = zeros(K_ox + 1, N + 1);         % Matrix containing the N+1 (time) vacancy concentration at the K_ox+1 (space) position in the oxide [/cm^3]
Ce_ox = zeros(K_ox + 1, N + 1);         % Matrix containing the N+1 (time) elecrons concentration at the K_ox+1 (space) position in the oxide [/cm^3]
Ch_ox = zeros(K_ox + 1, N + 1);         % Matrix containing the N+1 (time) hydrogen concentration at the K_ox+1 (space) position in the oxide [/cm^3]
Co_ox = zeros(K_ox + 1, N + 1);         % Matrix containing the N+1 (time) oxygen concentration at the K_ox+1 (space) position in the oxide [/cm^3]
xo_ox = zeros(K_ox + 1, N + 1);         % Matrix containing the N+1 (time) oxygen fraction at the K_ox+1 (space) position in the oxide [/cm^3]
Nb_needed = zeros(K_ox + 1, N + 1);     % Matrix containing the N+1 (time) amount of niobium needed to compensate the space charge at the K_ox+1 (space) position in the oxide [wt fraction]

% Create the result vector and the interface position vectors
res = zeros(K + 4, 1);                  % Time-evolving column vector containing the K+2 (space) oxygen concentration + the position of the oxide/alpha-phase interface + the position of the alpha-phase/beta-phase interface
res_save = zeros(K + 2, N + 1);         % Save the raw concentration results of the vector res without the interpolation
xi_ox = e_Zr * ones(1, N + 1);          % Vector containing the N+1 (time) position of the oxide/alpha-phase interface [cm]
xi_ab = e_Zr * ones(1, N + 1);          % Vector containing the N+1 (time) position of the alpha-phase/beta-phase interface [cm]

% Create the layer thickness vectors
d_protect = zeros(1, N + 1);            % Vector containing the N+1 (time) thickness of the protective oxide layer [cm]
d_total = zeros(1, N + 1);              % Vector containing the N+1 (time) thickness of the total oxide layer (protective + not protective) [cm]
duct_b = zeros(1, N+1);                 % Vector containing the N+1 (time) thickness of the ductile beta-phase layer for high-temperature corrosion [cm]

% Create the fluxes vectors
Jv = zeros(1, N + 1);                   % Vector containing the N+1 (time) oxygen flux in the oxide [/cm^2/s]
Je = zeros(1, N + 1);                   % Vector containing the N+1 (time) electron flux in the oxide [/cm^2/s]
Jh = zeros(1, N + 1);                   % Vector containing the N+1 (time) hydrogen flux in the oxide [/cm^2/s]
Jo = zeros(1, N + 1);                   % Vector containing the N+1 (time) oxygen flux through the oxide/alpha-phase interface [/cm^2/s]
%%%%%Jxo = zeros(1, N + 1);                  % Vector containing the N+1 (time) oxygen flux through the oxide/alpha-phase interface [/cm^2/s]
E = zeros(1, N + 1);                    % Vector containing the N+1 (time) homogeneous electric field in the oxide [V/cm]

% Create mass, weight gain and resistivity vectors
m = m_0 * ones(1, N + 1);               % Vector containing the N+1 (time) mass pf the sample [mg]
wg = zeros(1, N + 1);                   % Vector containing the N+1 (time) weight gain [mg/cm^2]
wg_int = zeros(1, N + 1);               % Vector containing the N+1 (time) weight gain based on integration of oxygen concentration [mg/cm^2]
rho = zeros(1, N + 1);                  % Vector containing the N+1 (time) resistivity of the oxide layer [ohm*cm]

% Create the hydrogen pickup fraction vectors
Dh = zeros(1, N + 1);                   % Vector containing the N+1 (time) hydrogen diffusion coefficient in the oxide [cm^2/s]
fh_inst = 0.15 * ones(1, N + 1);        % Vector containing the N+1 (time) instantaneous hydrogen pickup fraction []
fh_tot = zeros(1, N + 1);               % Vector containing the N+1 (time) total hydrogen pickup fraction []

% Initiate the transition vector
N_transition = 1;                       % Iteration at which the transition occurs (default at 1 to add 0 when the transition hasn't occured yet)
d_transition = 0;                       % Thickness of the oxide at the last transition [cm]


%% Initialisation : t = 0

dx_fixed = e_fixed / K_fixed;               % Spatial step size in the fixed mesh zone [cm]
dx_moving = (e_Zr - e_fixed) / K_moving;    % Spatial step size in the moving mesh zone [cm]

% Create the position vector of the fixed zone for all the times
for k = 1:K_fixed
    x(k, :) = (k - 1) * dx_fixed;
end
% Create the position vector of the moving zone for the first time
for k = 1:K_b + 1
    x(K_fixed + k, 1) = e_fixed + (k - 1) * dx_moving;
end
for k = 1:K_a + 1
    x(K_fixed + K_b + k + 1, 1) = x(K_fixed + K_b + 1, 1) + (k - 1) * dx_moving;
end
% (To prepare the arrival of the 2 phases in the moving mesh, the point K_fixed + K_b + 1 is the same as the point K_fixed + K_b + 2)

xi_ox(1) = e_Zr;    % Initial position of the oxide/alpha-phase interface
xi_ab(1) = e_Zr;    % Initial position of the alpha-/beta-phase interface

dx_a = 0;           % Artificially set the step size in the alpha-phase layer at 0
dx_b = 0;           % Artificially set the step size in the beta-phase layer at 0

% Set to defaults all the parameters related to oxygen ingress into the metals
if strcmp(model, 'C4') || strcmp(model, 'C4-H')
    Da = 0;                                                 % Oxygen diffusion coefficient in the alpha-phase [/cm^2/s]
    Db = 0;                                                 % Oxygen diffusion coefficient in the beta-phase [/cm^2/s]
end

%% First iteration : use isotherm method because the oxygen flux would be infinite : oxygen ingress into the metal is ignored for this first step

n = 1;      % Loop counter for the time step

C4_parameter_T;             % Compute the temperature-dependent parameters

% Case without preoxide: the current oxide thickness is 0
if preoxide == 0
    
    C4_ox_flux;                                     % Run the flux computation code
    
    e = sqrt(Vm * gammav * dt) / PBR;               % thickness of metal turned into oxide (<> from oxide thickness because of the PBR)
    
    xi_ox(2) = xi_ox(1) - e;                        % Update the position of the oxide/alpha-phase interface at the next time step
    xi_ab(2) = xi_ox(2);                            % The position of the alpha-/beta-phase interface is the same as the oxide/alpha-phase interface because there is no alpha phase yet
    
    d_protect(2) = PBR * (e_Zr - xi_ox(2));         % The thickness of the oxide layer is the size of the metal turned into oxide times the Pilling-Bedworth Ratio (PBR)
    d_total(2) = d_protect(2);                      % Total size of the oxide layer (no transition yet)
    dx_moving = (xi_ox(2) - e_fixed) / K_moving;    % Redefine the step size of the moving zone mesh taking into account the new position of the oxide/alpha-phase interface
    
    % Use this oxide thickness estimation to compute the different fluxes (used for the following step)
    n = 2;  % Second time step
    
    C4_parameter_T;             % Compute the temperature-dependent parameters
    
    Co(K + 2, n) = Co_a_ox;     % The concentration at the last position is fixed: it is the concentration in the oxide at the oxide/alpha-phase interface
    
    C4_ox_flux;                 % Run the flux computation code
    
    % Compute the species fluxes by dividing the constants calculated in C4_ox_flux by the current oxide thickness
    
    Jv(n) = gammav / d_protect(n);
    Je(n) = gammae / d_protect(n);
    Jh(n) = gammah / d_protect(n);
%     Jh(n) = 2 * fh_inst(n) * Jv(n);         % If we already know the hydrogen pickup (from experiment)
    E(n) = potential / d_protect(n);
    
    % put a preoxide to simulate previous corrosion
    % We are back to the first iteration n = 1
else
    
    d_protect(n) = preoxide;        % The current oxide thickness is the size of preoxide
    d_total(n) = d_protect(n);      % No transition yet
    
    Co(K + 2, n) = Co_a_ox;         % The concentration at the last position is fixed: it is the concentration in the oxide at the oxide/alpha-phase interface
    
    C4_ox_flux;                     % Run the flux computation code
    
    % Compute the species fluxes by dividing the constants calculated in C4_ox_flux by the current oxide thickness
    
    Jv(n) = gammav / d_protect(n);
    Je(n) = gammae / d_protect(n);
    Jh(n) = gammah / d_protect(n);
%     Jh(n) = 2 * 0.15 * Jv(n);       % If we already know the hydrogen pickup (from experiment)
    E(n) = potential / d_protect(n);
    
end         % End of the preoxide choice


% Compute the concentrations in the oxide

x_ox(:, n) = linspace(0, d_protect(n), K_ox + 1);   % Update the mesh in the oxide according to the current oxide thickness

for k_ox = 1:K_ox + 1
    Cv_ox(k_ox, n) = Cv_ox_a * exp(Zv * E(n) * x_ox(k_ox, n) / (kb * T(n))) + Jv(n) / (mu_v * E(n)) * (1 - exp(Zv * E(n) * x_ox(k_ox, n) / (kb * T(n))));
%     Cv_ox(k_ox, n) = 0;
    Co_ox(k_ox, n) = 2 * RhoZi * Na / MZi - Cv_ox(k_ox, n);
    Ce_ox(k_ox, n) = Ce_ox_a * exp(Ze * E(n) * x_ox(k_ox, n) / (kb * T(n))) + Je(n) / (mu_e * E(n)) * (1 - exp(Ze * E(n) * x_ox(k_ox, n) / (kb * T(n))));
    if strcmp(model, 'C4-H') || strcmp(model, 'C4-OH')
        Ch_ox(k_ox, n) = Ch_ox_a * exp(Zh * E(n) * x_ox(k_ox, n) / (kb * T(n))) - Jh(n) / (mu_h * E(n)) * (1 - exp(Zh * E(n) * x_ox(k_ox, n) / (kb * T(n))));     % - Jh because Jh is in the other direction
    end
    Nb_needed(k_ox, n) = MNb * (2 * Cv_ox(k_ox, n) - Ce_ox(k_ox, n) + Ch_ox(k_ox, n)) / ((4 - Nb_ox_state) * Na * RhoZi);
end

% Compute the weight gain from oxygen and hydrogen pickup
wg_o = Jv(n-1) * dt * Mo / Na * 1000;     % *1000 to convert it in mg
wg_h = Jh(n-1) * dt * Mh / Na * 1000;

%%% COMMENT MORE

wg(n) = wg_o + wg_h;        % There is no weight gain before
m(n) = m_0 + wg(n) * Area;
Ch(n) = 1e6 * (Ch_nat * 1e-6 * m_0 + wg_h * Area) / m_0;
fh_tot(n) = 1e-6 * (Ch(n) * m(n) - Ch_nat * m_0) / (2 * (m(n) - m_0) / 16);     % Compute the hydrogen pickup fraction
rho(n) = E(n) / (Je(n) * c_e);

% Compute the size of the ductile beta phase
if strcmp(model, 'C4-O') || strcmp(model, 'C4-OH')
    idx_duct_b = find(Co(:, n) > Co_duct_b, 1);
    if isempty(idx_duct_b)
        duct_b(n) = x(K+2, n);
    else
        duct_b(n) = x(idx_duct_b, n);
    end
end

% Remesh the moving zone of the mesh, accounting for the different oxide thickness

for k = 1:K_b + 1       % In the beta-phase
    x(K_fixed + k, n) = e_fixed + (k - 1) * dx_moving;
end
for k = 1:K_a + 1       % In the alpha-phase
    x(K_fixed + K_b + k + 1, n) = x(K_fixed + K_b + 1, n) + (k - 1) * dx_moving;
end
for k = 1:K + 2             % Update the result vector with the new concentrations
    res(k, 1) = Co(k, n);
end
res_save(:, n) = res(1:K + 2, 1);

%%%%%%%%%%%

wg_int(n) = (trapz(x(:,n),Co(:,n)) + trapz(x_ox(:,n),Co_ox(:,n)) - trapz(x(:,1),Co(:,1)) - trapz(x_ox(:,1),Co_ox(:,1))) * Mo / Na * 1000;
wg(n) = trapz(x_ox(:,n),Co_ox(:,n)) * Mo / Na * 1000;      % There is no computed flux for the first time step (n = 2), so we use the integration of the oxygen in the oxide as a weight gain (no diffusion in the metal yet)

% Update the result vector with the new interfaces positions
res(K + 3, 1) = xi_ox(n);
res(K + 4, 1) = xi_ab(n);

init = n + 1;       % Custom starting iteration depending if we used a preoxide or not


%% Iteration

for n = init:N+1    % Starting iteration depends if we used a preoxide or not
        
        C4_parameter_T;             % Compute the temperature-dependent parameters
        
        %%%%%%%%%%%%%%%%%%%%%%
        
        % Set the oxygen diffusion coefficient depending on the temperature
        %%{
        if (strcmp(model, 'C4-O') || strcmp(model, 'C4-OH')) && strcmp(mode, 'D HT optimization') == 0
            if T(n) == 633.15
                Da = 1.36e-15;
                Db = 4.8e-12;
            elseif T(n) == 1223.15
                Da = 4.3e-9;
                Db = 2.4e-7;
            elseif T(n) == 1273.15
%                 Da = 8.63516356504005e-09;  % LS optimized
%                 Db = 4.67287494996384e-06;  % LS optimized
%                 Da = 7.96256467581779e-09;  % LS 90%-alpha optimized
%                 Db = 5.41930221393628e-07;  % LS 90%-alpha optimized
%                 Da = 3.76533294987210e-09;  % LS alpha optimized
%                 Db = 3.79072975947571e-07;  % LS alpha optimized
                Da = 4.5e-9;                % Chosen parameter for publication
                Db = 3.8e-7;                % Chosen parameter for publication
            elseif T(n) == 1373.15
%                 Da = 2.91841589411482e-08;  % LS optimized
%                 Db = 1.57330557652387e-06;  % LS optimized
%                 Da = 2.38142736959768e-08;  % LS 90%-alpha optimized
%                 Db = 4.30821059739715e-07;  % LS 90%-alpha optimized
%                 Da = 2.34027770549067e-08;  % LS alpha optimized
%                 Db = 3.63968085565743e-07;  % LS alpha optimized
                Da = 2.6e-8;                % Chosen parameter for publication
                Db = 3.7e-7;                % Chosen parameter for publication
%                 Da = 2.48e-8;               % Corvolan
%                 Db = 8.44e-7;               % Corvolan
            elseif T(n) == 1473.15
%                 Da = 7.71767743625447e-08;  % LS optimized
%                 Db = 1.48345237506789e-07;  % LS optimized
%                 Da = 8.17305090977462e-08;  % LS alpha optimized
%                 Db = 1.56064975154785e-07;  % LS alpha optimized
%                 Da = 1.6e-7;                % Best d + wg
%                 Db = 1e-6;                  % Best d + wg
                Da = 1e-7;                  % Chosen parameter for publication
                Db = 6e-7;                  % Chosen parameter for publication
%                 Da = 8.5e-8;                % Corvolan
%                 Db = 1.55e-6;               % Corvolan
            elseif T(n) == 1573.15
%                 Da = 3.35180760866770e-07;  % LS optimized
%                 Db = 5.96370019459895e-06;  % LS optimized
                Da = 2.61152343691429e-07;  % LS alpha optimized
                Db = 2.01070471839367e-06;  % LS alpha optimized
            elseif T(n) == 1673.15
%                 Da = 9.67247347395458e-07;  % LS optimized
%                 Db = 9.56883554131739e-06;  % LS optimized
                Da = 8.72249840061976e-07;  % LS alpha optimized
                Db = 8.56187744202073e-06;  % LS alpha optimized
            elseif T(n) == 1773.15
%                 Da = 1.62917175797368e-06;  % LS optimized
%                 Db = 1.01412268582896e-05;  % LS optimized
                Da = 1.73131961527490e-06;  % LS alpha optimized
                Db = 9.27853860803972e-06;  % LS alpha optimized
            end
            %}
            
%             Da = 15 * exp(-24.91e3/(RR*T(n)));             % Zanella
%             Da = 0.196 * exp(-41000./(RR*T(n)));           % Mallett
%             Da = 0.196 * exp(-42500./(RR*T(n)));           % Mallett lower uncertainty bound
%             Da = 5.2 * exp(-50800./(RR*T(n)));             % Pemsler (400C --> 1500C)
%             Da = 7.3 * exp(-53330./(RR*T(n)));             % UW-Madison (1000C --> 1500C)
%             
%             Db = 0.0453 * exp(-28200./(RR*T(n)));          % Mallett
%             Db = 0.0263 * exp(-28200 ./ (RR * T(n)));      % Pawel
%             Db = 2.48e-2 * exp(-28200./(RR*T(n)));         % Perkins
%             Db = 1.71e-2 * exp(-29300./(RR*T(n)));         % Perkins lower uncertainty bound
%             Db = 0.11 * exp(-33390./(RR*T(n)));            % UW-Madison (1000C --> 1500C)
        end
    
    if dx_a == 0
        %% Temperature below the phase transition: single alpha-phase
        %See thesis p. 33-37
        
        % Coefficient used for the matrix creation
        A_fixed = Da * dt / dx_fixed^2;
        A_moving = Da * dt / dx_moving^2;
        B_moving = Da * dt / dx_moving / (PBR * Co_ox_a - Co_a_ox);
        C_mat = dt * Jv(n - 1) / (PBR * Co_ox_a - Co_a_ox);
        
        % Creation of the diagonal of the matrix
        diagon = horzcat((1 + 2 * A_fixed) * ones(1, K_fixed), (1 + 2 * A_moving) * ones(1, K_moving + 2), ones(1, 2));
        diagon(1) = 1;
        diagon(K_fixed + K_b + 1) = 1;
        diagon(K_fixed + K_b + 2) = 1;
        diagon(K + 2) = 1;
        
        low_diag = horzcat(-A_fixed * ones(1, K_fixed), -A_moving * ones(1, K_moving + 2), zeros(1));
        low_diag(K_fixed + K_b) = 0;
        low_diag(K_fixed + K_b + 1) = 0;
        low_diag(K + 1) = 0;
        low_diag(K + 2) = 0;
        
        up_diag = horzcat(-A_fixed * ones(1, K_fixed), -A_moving * ones(1, K_moving + 2), zeros(1));
        up_diag(1) = 0;
        up_diag(K_fixed + K_b + 1) = 0;
        up_diag(K_fixed + K_b + 2) = 0;
        up_diag(K + 2) = 0;
        
        % Assemble the diagonals to form the matrix
        M1 = diag(diagon, 0) + diag(low_diag, -1) + diag(up_diag, 1);
        M1(K_fixed + 1, K_fixed) = -2 * Da * dt / (dx_fixed * (dx_fixed + dx_moving));              % Correction for the term at the boundary between the fixed and the moving mesh
        M1(K_fixed + 1, K_fixed + 1) = 1 + 2 * Da * dt / (dx_fixed * dx_moving);                    % Correction for the term at the boundary between the fixed and the moving mesh
        M1(K_fixed + 1, K_fixed + 2) = -2 * Da * dt / (dx_moving * (dx_fixed + dx_moving));         % Correction for the term at the boundary between the fixed and the moving mesh
        M1(K + 3, K + 1) = B_moving;
        M1(K + 3, K + 2) = -B_moving;
        
        M2 = eye(K + 4);
        
        M3 = zeros(K + 4, 1);
        M3(K + 3, 1) = -C_mat;
        
        % Solve the system: M1 * y = M2 * x + M3
        res = M1 \ (M2 * res + M3);
        res_save(:, n) = res(1:K+2, 1);
        
        % Update the position of the interfaces and the oxide thickness
        xi_ox(n) = min(res(K + 3, 1), e_Zr);        % prevent any negative oxide thickness
        xi_ab(n) = xi_ox(n);
        d_protect(n) = preoxide + PBR * (e_Zr - xi_ox(n)) - d_transition;
        
        %% Transition
        
        if T(n) < T_transf
            %% Phase transition temperature not reached yet
            
            if d_protect(n) > transition(length(N_transition))
                %% Transition occurs
                
                d_transition = d_transition + transition(length(N_transition));     % Update the thickness of the oxide at the transition
                N_transition(end + 1) = n - 1;                                      % update the vector containing the iterations at which the transition occurs
                
                C4_ox_flux;                                                         % Run the flux computation code
                
                e = sqrt(Vm * gammav * dt) / PBR;                                   % thickness of metal turned into oxide (<> from oxide thickness because of the PBR)
                
                xi_ox(n) = e_Zr - d_transition / PBR - e;                           % Update the position of the oxide/alpha-phase interface
                d_protect(n) = PBR * (e_Zr - xi_ox(n)) - d_transition;              % The thickness of the oxide layer is the size of the metal turned into oxide times the Pilling-Bedworth Ratio (PBR)
                
            else
                %% No transition
                
                C4_ox_flux;                         % Run the flux computation code
            end
            
            % Compute the species fluxes by dividing the constants calculated in C4_ox_flux by the current oxide thickness
            
            Jv(n) = gammav / d_protect(n);
            Jo(n) = Da * (res(K + 2, 1) - res(K + 1, 1)) / dx_moving;
            Je(n) = gammae / d_protect(n);
            Jh(n) = gammah / d_protect(n);
%             Jh(n) = 2 * fh_inst(n) * Jv(n);         % If we already know the hydrogen pickup (from experiment)
            E(n) = potential / d_protect(n);
            
            
            dx_moving = (xi_ox(n) - e_fixed) / K_moving;        % Redefine the step size of the moving zone mesh taking into account the new position of the oxide/alpha-phase interface
            
            % Remesh the moving zone of the mesh, accounting for the different oxide thickness
            
            for k = 1:K_b + 1       % In the beta-phase
                x(K_fixed + k, n) = e_fixed + (k - 1) * dx_moving;
            end
            for k = 1:K_a + 1       % In the alpha-phase
                x(K_fixed + K_b + k + 1, n) = x(K_fixed + K_b + 1, n) + (k-1) * dx_moving;
            end
            
            % Update the concentration vector by taking the values from the result vector
            for k = 2:K_fixed
                Co(k, n) = res(k, 1);
            end
            
            % As the mesh is moving, the concentration on the new mesh are linearly interpolated from the old mesh
            for k = K_fixed + 1:K + 1                           % Go through all the points in the moving zone
                l = K_fixed;
                while l <= K + 2 && x(l, n - 1) < x(k, n)       % Find where is the point located in the old mesh
                    l = l + 1;
                end
                if l == K + 3                                   % If the point is after the last one in the old mesh, the concentration is the one at the interface
                    Co(k, n) = Co_a_ox;
                else
                    Co(k, n) = res(l - 1, 1) + (x(k, n) - x(l - 1, n - 1)) / (x(l, n - 1) - x(l - 1, n - 1)) * (res(l, 1) - res(l - 1, 1));     % Linear interpolation
                end
            end
            
            for k = K_fixed + 1:K + 1
                res(k, 1) = Co(k, n);                           % Update the result vector
            end
            
            Co(K + 2, n) = Co_a_ox;                             % The last point has the concentration at the interface oxide/alpha-phase
            res(K + 3, 1) = xi_ox(n);                           % Update the result vector
            
        else
            %% Phase transition: creation of the alpha and beta phase
            
            Co_interface = (Co_a_b + Co_b_a) / 2;       % We use the mean concentration between the two sides of the alpha-/beta-phase interface as the initial location of the interface
            
            % Find the position of the initial alpha-/beta-phase interface
            l = 1;
            while (l <= K + 2) && res(l, 1) < Co_interface
                l = l + 1;
            end
            
            if l == K + 3       % If no point in the metal as the concentration required to become an alpha-phase, the position of the alpha-/beta-phase is the same as the position of the oxide/alpha-phase
                xi_ab(n) = xi_ox(n);
            else
                xi_ab(n) = x(l - 1, n - 1) + (Co_interface - res(l - 1, 1)) / (res(l, 1) - res(l - 1, 1)) * (x(l, n - 1) - x(l - 1, n - 1));        % Linear interpolation to find the initial location of the interface
            end
            
            C4_ox_flux;         % Run the flux computation code
            
            % Compute the species fluxes by dividing the constants calculated in C4_ox_flux by the current oxide thickness
            
            if d_protect(n) ~= 0
                Jv(n) = gammav / d_protect(n);
                Jo(n) = Da * (res(K + 2, 1) - res(K + 1, 1)) / dx_a;
                Je(n) = gammae / d_protect(n);
                Jh(n) = gammah / d_protect(n);
%                 Jh(n) = 2 * fh_inst(n) * Jv(n);         % If we already know the hydrogen pickup (from experiment)
                E(n) = potential / d_protect(n);
            else
                Jv(n) = 0;
                Jo(n) = Da * (res(K + 2, 1) - res(K + 1, 1)) / dx_a;
                Je(n) = 0;
                Jh(n) = 0;
%                 Jh(n) = 2 * fh_inst(n) * Jv(n);         % If we already know the hydrogen pickup (from experiment)
                E(n) = 0;
            end
            
            
            % Create the mesh with new alpha and beta phases
            
            dx_a = (xi_ox(n) - xi_ab(n)) / K_a;     % Spatial step size in the alpha-phase
            dx_b = (xi_ab(n) - e_fixed) / K_b;      % Spatial step size in the beta-phase
            
            for k = 1:K_b + 1       % In the beta-phase
                x(K_fixed + k, n) = e_fixed + (k - 1) * dx_b;
            end
            for k = 1:K_a + 1       % In the alpha-phase
                x(K_fixed + K_b + k + 1, n) =xi_ab(n) + (k - 1) * dx_a;
            end
            
            % Update the concentration vector by taking the values from the result vector
            for k = 2:K_fixed
                Co(k, n) = res(k, 1);
            end
            
            %%%%%%%%%%%%%%%%% CHECK THE CHANGE IN THE INTERFACE CONCENTRATION WHEN PHASE CREATION
            
            % As the mesh is moving, the concentrations on the new mesh are linearly interpolated from the old mesh
            for k = K_fixed + 1:K_fixed + K_b + 1                       % Go through all the points in the moving zone
                l = K_fixed;
                while l <= K_fixed + K_b + 2 && x(l, n - 1) < x(k, n)   % Find where is the point located in the old mesh
                    l = l + 1;
                end
                if l == K_fixed + K_b + 3                               % If the point is after the last one in the old mesh, the concentration is the one at the interface
                    Co(k, n) = Co_nat;
                elseif l == K_fixed                                     % If the point is after the first one in the old mesh, the concentration is the one at the boundary
                    Co(k, n) = Co(K_fixed, n);
                else
                    Co(k, n) = res(l - 1, 1) + (x(k, n) - x(l - 1, n - 1)) / (x(l, n - 1) - x(l - 1, n - 1)) * (res(l, 1) - res(l - 1, 1));     % Linear interpolation
                end
            end
            
            % Redefine the concentration at the alpha-/beta-phase interface
            Co(K_fixed + K_b + 1, n) = Co_b_a;
            Co(K_fixed + K_b + 2, n) = Co_a_b;
            
            for k = K_fixed + K_b + 3:K + 1
                l = K_fixed + K_b + 2;
                while l <= K + 2 && x(l, n - 1) < x(k, n)
                    l = l + 1;
                end
                if l == K + 3                                           % If the point is after the last one in the old mesh, the concentration is the one at the interface
                    Co(k, n) = Co_a_ox;
                elseif l == K_fixed + K_b + 2                           % If the point is after the first one in the old mesh, the concentration is the one at the interface
                    Co(k, n) = Co_a_b;
                else
                    Co(k, n) = res(l - 1, 1) + (x(k, n) - x(l - 1, n - 1)) / (x(l, n - 1) - x(l - 1, n - 1)) * (res(l, 1) - res(l - 1, 1));     % Linear interpolation
                end
            end
            
            Co(K + 2, n) = Co_a_ox;                 % The last point has the concentration at the interface oxide/alpha-phase
            
            for k = 1:K + 2
                res(k, 1) = Co(k, n);               % Update the result vector
            end
            
            % Update the result vector
            res(K + 3, 1) = xi_ox(n);
            res(K + 4, 1) = xi_ab(n);
            
            
        end
        
    else
        if T(n) > T_transf
            %% Temperature above phase transition: alpha and beta phases already existing
            %See thesis p. 67-70
            
            % If during a fast tempererature transient the oxide layer gets totally dissolved into the metal
            if d_protect(n-1) == 0
                C4_ox_flux;                                     % Run the flux computation code
                
                e = sqrt(Vm * gammav * dt) / PBR;               % thickness of metal turned into oxide (<> from oxide thickness because of the PBR)
                
                xi_ox(n) = xi_ox(n-1) - e;                      % Update the position of the oxide/alpha-phase interface at the next time step
                xi_ab(n) = xi_ab(n-1);                          % The position of the alpha-/beta-phase interface is the same as the oxide/alpha-phase interface because there is no alpha phase yet
                
                d_protect(n) = PBR * (e_Zr - xi_ox(n));         % The thickness of the oxide layer is the size of the metal turned into oxide times the Pilling-Bedworth Ratio (PBR)
                
            else
                
                % Coefficient used for the matrix creation
                A_a = Da * dt / dx_a^2;
                A_b = Db * dt / dx_b^2;
                B_ox = Da * dt / (dx_a * (PBR * Co_ox_a - Co_a_ox));
                C_mat = dt * Jv(n - 1) / (PBR * Co_ox_a - Co_a_ox);
                B_a = Da * dt / (dx_a * (Co_a_b - Co_b_a));
                B_b = Db * dt / (dx_b * (Co_a_b - Co_b_a));
                
                % Creation of the diagonal of the matrix
                diagon = horzcat((1 + 2 * A_fixed) * ones(1, K_fixed), (1 + 2 * A_b) * ones(1, K_b + 1), (1 + 2 * A_a) * ones(1, K_a + 1), ones(1, 2));
                diagon(1) = 1;
                diagon(K_fixed + K_b + 1) = 1;
                diagon(K_fixed + K_b + 2) = 1;
                diagon(K + 2) = 1;
                
                low_diag = horzcat(-A_fixed * ones(1, K_fixed), -A_b * ones(1, K_b + 1), -A_a * ones(1, K_a + 1), zeros(1));
                low_diag(K_fixed + K_b) = 0;
                low_diag(K_fixed + K_b + 1) = 0;
                low_diag(K + 1) = 0;
                low_diag(K + 2) = 0;
                
                up_diag = horzcat(-A_fixed * ones(1, K_fixed), -A_b * ones(1, K_b + 1), -A_a * ones(1, K_a + 1), zeros(1));
                up_diag(1) = 0;
                up_diag(K_fixed + K_b + 1) = 0;
                up_diag(K_fixed + K_b + 2) = 0;
                up_diag(K + 2) = 0;
                
                % Assemble the diagonals to form the matrix
                M1 = diag(diagon, 0) + diag(low_diag, -1) + diag(up_diag, 1);
                M1(K_fixed + 1, K_fixed) = -2 * Db * dt / (dx_fixed * (dx_fixed + dx_b));               % Correction for the term at the boundary between the fixed and the moving mesh
                M1(K_fixed + 1, K_fixed + 1) = 1 + 2 * Db * dt / (dx_fixed * dx_b);                     % Correction for the term at the boundary between the fixed and the moving mesh
                M1(K_fixed + 1, K_fixed + 2) = -2 * Db * dt / (dx_b * (dx_fixed + dx_b));               % Correction for the term at the boundary between the fixed and the moving mesh
                M1(K + 3, K + 1) = B_ox;
                M1(K + 3, K + 2) = -B_ox;
                M1(K + 4, K_fixed + K_b) = B_b;
                M1(K + 4, K_fixed + K_b + 1) = -B_b;
                M1(K + 4, K_fixed + K_b + 2) = -B_a;
                M1(K + 4, K_fixed + K_b + 3) = B_a;
                
                M2 = eye(K + 4);
                
                M3 = zeros(K + 4, 1);
                M3(K + 3, 1) = -C_mat;
                
                % Solve the system
                res = M1 \ (M2 * res + M3);
                res_save(:, n) = res(1:K+2, 1);
                
                % Update the position of the interfaces and the oxide thickness
                xi_ox(n) = min(res(K + 3, 1), e_Zr);
                xi_ab(n) = min(res(K + 4, 1), xi_ox(n));
                d_protect(n) = preoxide + PBR * (e_Zr - xi_ox(n));
                
                C4_ox_flux;             % Run the flux computation code
                
                % Compute the species fluxes by dividing the constants calulated in C4_ox_flux by the current oxide thickness
                if d_protect(n) ~= 0
                    Jv(n) = gammav / d_protect(n);
                    Jo(n) = Da * (res(K + 2, 1) - res(K + 1, 1)) / dx_a;
                    Je(n) = gammae / d_protect(n);
                    Jh(n) = gammah / d_protect(n);
%                     Jh(n) = 2 * fh_inst(n) * Jv(n);         % If we already know the hydrogen pickup (from experiment)
                    E(n) = potential / d_protect(n);
                else
                    Jv(n) = 0;
                    Jo(n) = Da * (res(K + 2, 1) - res(K + 1, 1)) / dx_a;
                    Je(n) = 0;
                    Jh(n) = 0;
%                     Jh(n) = 2 * fh_inst(n) * Jv(n);         % If we already know the hydrogen pickup (from experiment)
                    E(n) = 0;
                end
            
            end
            
            % Update the mesh size in the alpha- and beta-phases
            dx_a = (xi_ox(n) - xi_ab(n)) / K_a;
            dx_b = (xi_ab(n) - e_fixed) / K_b;
            
            for k = 1:K_b + 1       % In the beta-phase
                x(K_fixed + k, n) = e_fixed + (k - 1) * dx_b;
            end
            for k = 1:K_a + 1       % In the alpha-phase
                x(K_fixed + K_b + k + 1, n) =xi_ab(n) + (k - 1) * dx_a;
            end
            
            % Update the concentration vector by taking the values from the result vector
            for k = 2:K_fixed
                Co(k, n) = res(k, 1);
            end
            
            % As the mesh is moving, the concentrations on the new mesh are linearly interpolated from the old mesh
            for k = K_fixed + 1:K_fixed + K_b + 1                       % Go through all the points in the moving zone
                l = K_fixed;
                while l <= K_fixed + K_b + 2 && x(l, n - 1) < x(k, n)   % Find where is the point located in the old mesh
                    l = l + 1;
                end
                if l == K_fixed + K_b + 3                               % If the point is after the last one in the old mesh, the concentration is the one at the interface
                    Co(k, n) = Co_b_a;
                elseif l == K_fixed                                     % If the point is after the first one in the old mesh, the concentration is the one at the boundary
                    Co(k, n) = Co(K_fixed, n);
                else
                    Co(k, n) = res(l - 1, 1) + (x(k, n) - x(l - 1, n - 1)) / (x(l, n - 1) - x(l - 1, n - 1)) * (res(l, 1) - res(l - 1, 1));     % Linear interpolation
                end
            end
            
            % Redefine the concentration at the alpha-/beta-phase interface
            Co(K_fixed + K_b + 1, n) = Co_b_a;
            Co(K_fixed + K_b + 2, n) = Co_a_b;
            
            for k = K_fixed + K_b + 3:K + 1
                l = K_fixed + K_b + 2;
                while l <= K + 2 && x(l, n - 1) < x(k, n)
                    l = l + 1;
                end
                if l == K + 3                                           % If the point is after the last one in the old mesh, the concentration is the one at the interface
                    Co(k, n) = Co_a_ox;
                elseif l == K_fixed + K_b + 2                           % If the point is after the first one in the old mesh, the concentration is the one at the interface
                    Co(k, n) = Co_a_b;
                else
                    Co(k, n) = res(l - 1, 1) + (x(k, n) - x(l - 1, n - 1)) / (x(l, n - 1) - x(l - 1, n - 1)) * (res(l, 1) - res(l - 1, 1));     % Linear interpolation
                end
            end
            
            Co(K + 2, n) = Co_a_ox;         % The last point has the concentration at the interface oxide/alpha-phase
            
            for k = 1:K + 2
                res(k, 1) = Co(k, n);       % Update the result vector
            end
            
            % Update the result vector
            res(K + 3, 1) = xi_ox(n);
            res(K + 4, 1) = xi_ab(n);
            
        else
            %% Temperature below phase transition during cooling down: alpha and beta phases already existing, but alpha/beta interface not moving
            
            % The oxygen concentration at the alpha/beta interface are kept 
            xo_a_b = (45.86e-3*(T_transf-1136.15) - 44.77e-6*(T_transf-1136.15)^2 + 17.40e-9*(T_transf-1136.15)^3) * 1e-2;
            Co_a_b = xo_a_b * RhoZr * Na / (1 - xo_a_b) / MZr;                  % Concentration of oxygen in the alpha-phase at the alpha-/beta-phase interface [at/cm^3]
            xo_b_a = (9.59e-3*(T_transf-1136.15) + 4.72e-6*(T_transf-1136.15)^2 - 4.35e-9*(T_transf-1136.15)^3) * 1e-2;
            Co_b_a = xo_b_a * RhoZr * Na / (1 - xo_b_a) / MZr;                  % Concentration of oxygen in the beta-phase at the alpha-/beta-phase interface [at/cm^3]
        
            % Coefficient used for the matrix creation
            A_a = Da * dt / dx_a^2;
            A_b = Da * dt / dx_b^2;     % The oxygen diffusion coefficient in the alpha phase is also used for the prior-beta phase
            B_ox = Da * dt / (dx_a * (PBR * Co_ox_a - Co_a_ox));
            C_mat = dt * Jv(n - 1) / (PBR * Co_ox_a - Co_a_ox);
            
            % Creation of the diagonal of the matrix
            diagon = horzcat((1 + 2 * A_fixed) * ones(1, K_fixed), (1 + 2 * A_b) * ones(1, K_b + 1), (1 + 2 * A_a) * ones(1, K_a + 1), ones(1, 2));
            diagon(1) = 1;
            diagon(K_fixed + K_b + 1) = 1;
            diagon(K_fixed + K_b + 2) = 1;
            diagon(K + 2) = 1;
            
            low_diag = horzcat(-A_fixed * ones(1, K_fixed), -A_b * ones(1, K_b + 1), -A_a * ones(1, K_a + 1), zeros(1));
            low_diag(K_fixed + K_b) = 0;
            low_diag(K_fixed + K_b + 1) = 0;
            low_diag(K + 1) = 0;
            low_diag(K + 2) = 0;
            
            up_diag = horzcat(-A_fixed * ones(1, K_fixed), -A_b * ones(1, K_b + 1), -A_a * ones(1, K_a + 1), zeros(1));
            up_diag(1) = 0;
            up_diag(K_fixed + K_b + 1) = 0;
            up_diag(K_fixed + K_b + 2) = 0;
            up_diag(K + 2) = 0;
            
            % Assemble the diagonals to form the matrix
            M1 = diag(diagon, 0) + diag(low_diag, -1) + diag(up_diag, 1);
            M1(K_fixed + 1, K_fixed) = -2 * Db * dt / (dx_fixed * (dx_fixed + dx_b));               % Correction for the term at the boundary between the fixed and the moving mesh
            M1(K_fixed + 1, K_fixed + 1) = 1 + 2 * Db * dt / (dx_fixed * dx_b);                     % Correction for the term at the boundary between the fixed and the moving mesh
            M1(K_fixed + 1, K_fixed + 2) = -2 * Db * dt / (dx_b * (dx_fixed + dx_b));               % Correction for the term at the boundary between the fixed and the moving mesh
            M1(K + 3, K + 1) = B_ox;
            M1(K + 3, K + 2) = -B_ox;
            
            M2 = eye(K + 4);
            
            M3 = zeros(K + 4, 1);
            M3(K + 3, 1) = -C_mat;
            
            % Solve the system
            res = M1 \ (M2 * res + M3);
            res_save(:, n) = res(1:K+2, 1);
            
            % Update the position of the interfaces and the oxide thickness
            xi_ox(n) = min(res(K + 3, 1), e_Zr);
            xi_ab(n) = min(res(K + 4, 1), xi_ox(n));
            d_protect(n) = preoxide + PBR * (e_Zr - xi_ox(n));
            
            C4_ox_flux;             % Run the flux computation code
            
            % Compute the species fluxes by dividing the constants calulated in C4_ox_flux by the current oxide thickness
            Jv(n) = gammav / d_protect(n);
            Jo(n) = Da * (res(K + 2, 1) - res(K + 1, 1)) / dx_a;
            Je(n) = gammae / d_protect(n);
            Jh(n) = gammah / d_protect(n);
            %         Jh(n) = 2 * fh_inst(n) * Jv(n);         % If we already know the hydrogen pickup (from experiment)
            E(n) = potential / d_protect(n);
            
            % Update the mesh size in the alpha- and beta-phases
            dx_a = (xi_ox(n) - xi_ab(n)) / K_a;
            dx_b = (xi_ab(n) - e_fixed) / K_b;
            
            for k = 1:K_b + 1       % In the beta-phase
                x(K_fixed + k, n) = e_fixed + (k - 1) * dx_b;
            end
            for k = 1:K_a + 1       % In the alpha-phase
                x(K_fixed + K_b + k + 1, n) =xi_ab(n) + (k - 1) * dx_a;
            end
            
            % Update the concentration vector by taking the values from the result vector
            for k = 2:K_fixed
                Co(k, n) = res(k, 1);
            end
            
            % As the mesh is moving, the concentrations on the new mesh are linearly interpolated from the old mesh
            for k = K_fixed + 1:K_fixed + K_b + 1                       % Go through all the points in the moving zone
                l = K_fixed;
                while l <= K_fixed + K_b + 2 && x(l, n - 1) < x(k, n)   % Find where is the point located in the old mesh
                    l = l + 1;
                end
                if l == K_fixed + K_b + 3                               % If the point is after the last one in the old mesh, the concentration is the one at the interface
                    Co(k, n) = Co_b_a;
                elseif l == K_fixed                                     % If the point is after the first one in the old mesh, the concentration is the one at the boundary
                    Co(k, n) = Co(K_fixed, n);
                else
                    Co(k, n) = res(l - 1, 1) + (x(k, n) - x(l - 1, n - 1)) / (x(l, n - 1) - x(l - 1, n - 1)) * (res(l, 1) - res(l - 1, 1));     % Linear interpolation
                end
            end
            
            % Redefine the concentration at the alpha-/beta-phase interface
            Co(K_fixed + K_b + 1, n) = Co_b_a;
            Co(K_fixed + K_b + 2, n) = Co_a_b;
            
            for k = K_fixed + K_b + 3:K + 1
                l = K_fixed + K_b + 2;
                while l <= K + 2 && x(l, n - 1) < x(k, n)
                    l = l + 1;
                end
                if l == K + 3                                           % If the point is after the last one in the old mesh, the concentration is the one at the interface
                    Co(k, n) = Co_a_ox;
                elseif l == K_fixed + K_b + 2                           % If the point is after the first one in the old mesh, the concentration is the one at the interface
                    Co(k, n) = Co_a_b;
                else
                    Co(k, n) = res(l - 1, 1) + (x(k, n) - x(l - 1, n - 1)) / (x(l, n - 1) - x(l - 1, n - 1)) * (res(l, 1) - res(l - 1, 1));     % Linear interpolation
                end
            end
            
            Co(K + 2, n) = Co_a_ox;         % The last point has the concentration at the interface oxide/alpha-phase
            
            for k = 1:K + 2
                res(k, 1) = Co(k, n);       % Update the result vector
            end
            
            % Update the result vector
            res(K + 3, 1) = xi_ox(n);
            res(K + 4, 1) = xi_ab(n);
        end
    end
    
    
    d_total(n) = d_transition + d_protect(n);
    
    % Compute the concentration in the oxide
    
    x_ox(:, n) = linspace(0, d_protect(n), K_ox+1);     % Update the mesh in the oxide according to the current oxide thickness
    
    for k_ox = 1:K_ox+1
        Cv_ox(k_ox, n) = Cv_ox_a * exp(Zv * E(n) * x_ox(k_ox, n) / (kb * T(n))) + Jv(n) / (mu_v * E(n)) * (1 - exp(Zv * E(n) * x_ox(k_ox, n) / (kb * T(n))));
%         Cv_ox(k_ox, n) = 0;
        Co_ox(k_ox, n) = 2 * RhoZi * Na / MZi - Cv_ox(k_ox, n);
        Ce_ox(k_ox, n) = Ce_ox_a * exp(Ze * E(n) * x_ox(k_ox, n) / (kb * T(n))) + Je(n) / (mu_e * E(n)) * (1 - exp(Ze * E(n) * x_ox(k_ox, n) / (kb* T(n))));
        if strcmp(model, 'C4-H') || strcmp(model, 'C4-OH')
            Ch_ox(k_ox, n) = Ch_ox_a * exp(Zh * E(n) * x_ox(k_ox, n) / (kb * T(n))) - Jh(n) / (mu_h * E(n)) * (1 - exp(Zh * E(n) * x_ox(k_ox, n) / (kb * T(n))));     % - Jh because Jh is in the other direction
        end
        Nb_needed(k_ox, n) = MNb * (2 * Cv_ox(k_ox, n) - Ce_ox(k_ox, n) + Ch_ox(k_ox, n)) / ((4 - Nb_ox_state) * Na * RhoZi);
    end
    
    % Compute the weight gain from oxygen and hydrogen pickup
    wg_o = Jv(n-1) * dt * Mo / Na * 1e3;      % *1e3 to convert it in mg
    wg_h = Jh(n-1) * dt * Mh / Na * 1e3;
    
    wg(n) = wg(n - 1) + wg_o + wg_h;
    wg_int(n) = (trapz(x(:,n), Co(:,n)) + trapz(x_ox(:,n), Co_ox(:,n)) + Co_ox(end,n)*(d_total(n)-d_protect(n)) - trapz(x(:,1),Co(:,1)) - trapz(x_ox(:,1),Co_ox(:,1))) * Mo / Na * 1000;    % Add the contribution of the non-protective oxide layer to the weight gain
    m(n) = m_0 + wg(n) * Area;        % wg * surface area in cm^2
    
    % Compute the hydrogen concentration, hydrogen pickup fraction and the resistivity of the oxide layer
    Ch(n) = 1e6 * (Ch(n - 1) * 1e-6 * m(n - 1) + wg_h * Area) / m(n);
    fh_tot(n) = 1e-6 * (Ch(n) * m(n) - Ch_nat * m_0) / (2 * (m(n) - m_0) / 16);
    rho(n) = E(n) / (Je(n) * c_e);
    
    % Compute the size of the ductile beta phase
    if strcmp(model, 'C4-O') || strcmp(model, 'C4-OH')
        idx_duct_b = find(Co(:, n) > Co_duct_b, 1);
        if isempty(idx_duct_b)
            duct_b(n) = x(K+2, n);
        else
            duct_b(n) = x(idx_duct_b, n);
        end
    end
end         % End of the iterations over the time

xo = 1 ./ (1 + RhoZr * Na ./ (MZr .* Co));          % Convert oxygen concentration into oxygen fraction
res_plot = 1 ./ (1 + RhoZr * Na ./ (MZr .* res_save));
xo_ox = 1 ./ (1 + RhoZi * Na ./ (MZi .* Co_ox));    % Convert oxygen concentration into oxygen fraction
alpha = xi_ox - xi_ab;                              % Size of the alpha phase
wg_approx = d_total * RhoZi * 2*Mo / MZi * 1e3;     % Approximate weight gain (if neglect oxygen ingress into the metal) [mg/cm^2] ; d in [cm]
fh_inst = 0.5 * Jh ./ Jv;                           % Instantaneous hydrogen pickup fraction


% Save the data in the output file
save(strcat('output_file/', file_name), '-regexp', '^(?!(progress)$).')
