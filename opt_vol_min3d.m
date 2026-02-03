% ---------------------------------------------------------------------------------------
% To run the code use: opt_vol_min3d(80, 50, 4, 2.3, 3.0, 1.5)
% 3D Topology Optimization: Minimize Volume under Compliance Constraint
% This code is written by Md Zakirul et. al. for educational purpose only
% ---------------------------------------------------------------------------------------
function opt_vol_min3d(nelx,nely,nelz,Ctarget,penal,rmin)
%% Initialize User Configuration
maxIter = 500;       % Maximum iterations number
tol = 0.01;          % Convergence tolerance
show_Plot = 1;       % Enable real-time visualization
%% Define Material Properties
E0 = 1;              % Young's modulus (solid)
Emin = 1e-9;         % Young's modulus (void)
nu = 0.3;            % Poisson's ratio
%% Define Degrees of Freedom (Loads)
mid_j_index = floor((nely + 1)/2); % Middle node index
jl = mid_j_index - 1;              % Convert to 0-based index
[il,jl,kl] = meshgrid(nelx, jl, 0:nelz);                % Node coordinates (x, y, z)
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Node IDs
loaddof = 3*loadnid(:) - 1;                             % Extract Y-direction DOFs
%% Define Fixed Degrees of Freedom (Supports)
[iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                % Node coordinates of fixed face
fixednid = kf*(nelx+1)*(nely+1) + iif*(nely+1) + (nely+1-jf); % 1D node IDs of fixed nodes
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % X, Y, Z DOFs for each fixed node
%% Initialize Finite Element Model
nele = nelx * nely * nelz;                   % Total number of elements
ndof = 3 * (nelx+1) * (nely+1) * (nelz+1);   % Total degrees of freedom
F = sparse(loaddof, 1, -1, ndof, 1);         % External force vector (unit load in -Y)
U = zeros(ndof, 1);                          % Initialize displacement vector
freedofs = setdiff(1:ndof, fixeddof);        % Identify free (unconstrained) DOFs
KE = lk_H8(nu);                              % Element stiffness matrix for an 8-node hexahedral solid
%% Construct Global Element Connectivity
nodegrd = reshape(1:(nely+1)*(nelx+1), nely+1, nelx+1); % 2D node grid (XY-plane, layer z=0)
nodeids = reshape(nodegrd(1:end-1, 1:end-1), nely*nelx, 1); % Bottom-left node of each XY element
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);   % Offset for each layer in Z
nodeids = repmat(nodeids, 1, numel(nodeidz)) + ...
          repmat(nodeidz(:).', numel(nodeids), 1);      % Replicate XY pattern across Z-layers
edofVec = 3 * nodeids(:) + 1;                           % Base DOF index (X-component) for each element
%% Define Local DOF Ordering for A Single Hexahedral Element (24 Entries)
edof_pattern = [0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
                3*(nely+1)*(nelx+1) + ...
                [0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]];
%% Expand to All Elements
edofMat = repmat(edofVec, 1, 24) + repmat(edof_pattern, nele, 1);
%% Prepare Sparse Matrix Indexing for Global Stiffness Assembly
iK = reshape(kron(edofMat, ones(24,1))', 24*24*nele, 1); % Row indices
jK = reshape(kron(edofMat, ones(1,24))', 24*24*nele, 1); % Column indices
%% Sensitivity Filter
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for k1 = 1:nelz
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (k1-1)*nelx*nely + (i1-1)*nely + j1;
            for k2 = max(k1-(ceil(rmin)-1),1):min(k1+(ceil(rmin)-1),nelz)
                for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
                    for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                        e2 = (k2-1)*nelx*nely + (i2-1)*nely + j2;
                        k = k + 1;
                        iH(k) = e1;
                        jH(k) = e2;
                        sH(k) = max(0, rmin - sqrt((i1-i2)^2 + (j1-j2)^2 + (k1-k2)^2));
                    end
                end
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% Compute Initial Compliance C0
x = ones(nely,nelx,nelz);   % Start with full material
xPhys = x;
sK = reshape(KE(:)*(Emin + xPhys(:)'.^penal*(E0-Emin)), 24*24*nele, 1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
ce = reshape(sum((U(edofMat)*KE) .* U(edofMat), 2), [nely,nelx,nelz]);
C0 = sum(sum(sum((Emin + xPhys.^penal*(E0-Emin)) .* ce)));
Cmax = C0 * Ctarget;        % Maximum target compliance
fprintf('Initial compliance C0: %.4f\n', C0);
fprintf('Target compliance: %.4f (%.2fx C0)\n', Cmax, Ctarget);
%% Optimization Initialization
x = ones(nely,nelx,nelz) * 0.5;  % Start with 50% density
xPhys = x; loop = 0; change = 1;
%% Iteration History
C_ratio_hist = zeros(maxIter, 1);
V_hist = zeros(maxIter, 1);
%% Start Iteration Loop
while change > tol && loop < maxIter
    loop = loop + 1;
    %% FE-Analyisis
    sK = reshape(KE(:)*(Emin + xPhys(:)'.^penal*(E0-Emin)), 24*24*nele, 1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);
    %% Compute Compliance & Strain Energy for Sensitivity Analysis
    ce = reshape(sum((U(edofMat)*KE) .* U(edofMat), 2), [nely,nelx,nelz]);
    c = sum(sum(sum((Emin + xPhys.^penal*(E0-Emin)) .* ce)));
    %% Sensitivities
    dv = ones(nely,nelx,nelz);                    % Volume sensitivity
    dc = -penal*(E0-Emin)*xPhys.^(penal-1) .* ce; % Compliance sensitivity
    %% Filtering of Sensitivites
    dv(:) = H*(dv(:)./Hs);
    dc(:) = H*(dc(:)./Hs);
    %% Density-Based Design Update with Compliance Constraint
    l1 = 0; l2 = 1e9; move = 0.2;
    while (l2 - l1)/(l1 + l2) > 1e-3
        lmid = 0.5*(l2 + l1);
        %% Update Rule: x_new = x * sqrt(-lambda*dc/dv)
        xnew = max(0, ...
                   max(x - move, ...
                       min(1, ...
                           min(x + move, ...
                               x .* sqrt(-lmid * dc ./ dv)))));
        xPhys(:) = (H * xnew(:)) ./ Hs;
        %% Check Compliance Constraint
        sK_temp = reshape(KE(:)*(Emin + xPhys(:)'.^penal*(E0-Emin)), 24*24*nele, 1);
        K_temp = sparse(iK,jK,sK_temp); K_temp = (K_temp + K_temp')/2;
        U_temp = zeros(ndof,1);
        U_temp(freedofs,:) = K_temp(freedofs,freedofs) \ F(freedofs,:);
        ce_temp = reshape(sum((U_temp(edofMat)*KE) .* U_temp(edofMat), 2), [nely,nelx,nelz]);
        c_temp = sum(sum(sum((Emin + xPhys.^penal*(E0-Emin)) .* ce_temp)));
        %% Adjust Lagrange Multiplier Based on Constraint
        if c_temp > Cmax
            l1 = lmid;
        else
            l2 = lmid;
        end
    end
    change = max(abs(xnew(:) - x(:)));
    x = xnew;
    %% Compute Volume Fraction and Compliance Ratio
    vol = mean(xPhys(:));
    c_ratio = c / C0;
    %% Record History
    C_ratio_hist(loop) = c_ratio;
    V_hist(loop) = vol;
    %% Print Results
    fprintf(' It.:%5i Vol.:%7.3f Comp.:%11.4f C/C0:%.3f (target:%.2f) ch.:%.4f\n', ...
            loop, vol, c, c_ratio, Ctarget, change);
    %% Plot Densities
    if show_Plot
        clf;
        display_3D(xPhys);
    end
end
%% Final Visualization
clf; display_3D(xPhys);
fprintf('\nOptimization completed:\n');
fprintf('Final volume fraction: %.3f\n', mean(xPhys(:)));
fprintf('Final compliance: %.4f\n', c);
fprintf('Compliance ratio C/C0: %.3f (target: %.2f)\n', c/C0, Ctarget);
%% Historical Plot
iter = 1:loop;
figure; clf;
set(gcf,'Color','w','Units','centimeters','Position',[3 3 18 12]);
yyaxis left
h1 = plot(iter, C_ratio_hist(1:loop), '-', 'Color', [0.00 0.45 0.74], 'MarkerFaceColor', [0.00 0.45 0.74], ...
    'MarkerEdgeColor', [0.00 0.45 0.74], 'MarkerSize', 6, 'LineWidth', 2);
hold on;
yline_left = yline(Ctarget, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 2);
set(gca,'YScale','log');
ylabel('Compliance Ratio', 'FontSize', 30, 'FontWeight','bold');
set(gca, 'FontSize', 30, 'LineWidth', 2, 'TickLength', [0.02 0.02], 'YColor', [0.00 0.45 0.74], 'XColor', 'k');
yyaxis right
h2 = plot(iter, V_hist(1:loop), '-', 'Color', [0.49 0.18 0.56], 'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', [0.49 0.18 0.56], 'MarkerSize', 6, 'LineWidth', 2);
ylabel('Volume Fraction', 'FontSize', 30, 'FontWeight','bold');
set(gca, 'FontSize', 30, 'LineWidth', 2, 'TickLength', [0.02 0.02], 'YColor', [0.49 0.18 0.56]);
xlabel('Iteration', 'FontSize', 30, 'FontWeight','bold');
title('Convergence History', 'FontSize', 30, 'FontWeight','bold');
legend([h1, h2, yline_left], ...
    {'Compliance Ratio $C/C_0$','Volume Fraction','Target Compliance'}, 'Location','best', 'FontSize', 24, 'Box','on','Interpreter','latex');
cmin = min(C_ratio_hist(1:loop));
cmax = max(C_ratio_hist(1:loop));
vmin = min(V_hist(1:loop));
vmax = max(V_hist(1:loop));
yyaxis left
ylim([cmin*0.9, cmax*1.1]);
yyaxis right
ylim([vmin - 0.05*(vmax-vmin), vmax + 0.05*(vmax-vmin)]);
grid on;
yyaxis left
set(gca, 'XMinorGrid','off', 'YMinorGrid','off', 'XMinorTick','off', 'YMinorTick','off', 'Box','on', 'Layer','top');
yyaxis right
set(gca, 'XMinorGrid','off', 'YMinorGrid','off', 'XMinorTick','off', 'YMinorTick','off', 'Box','on', 'Layer','top');
drawnow;
end
%% Element Stiffness Matrix (8-node Hex)
function [KE] = lk_H8(nu)
    A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
        -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
    k = (A' * [1; nu]) / 144;
    K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
          k(2) k(1) k(2) k(4) k(6) k(7);
          k(2) k(2) k(1) k(4) k(7) k(6);
          k(3) k(4) k(4) k(1) k(8) k(8);
          k(5) k(6) k(7) k(8) k(1) k(2);
          k(5) k(7) k(6) k(8) k(2) k(1)];
    K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
          k(8)  k(9)  k(12) k(5)  k(3)  k(5);
          k(10) k(10) k(13) k(7)  k(4)  k(6);
          k(6)  k(5)  k(11) k(9)  k(2)  k(10);
          k(4)  k(3)  k(5)  k(2)  k(9)  k(12);
          k(11) k(4)  k(6)  k(12) k(10) k(13)];
    K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
          k(7)  k(6)  k(4)  k(10) k(13) k(10);
          k(5)  k(5)  k(3)  k(8)  k(12) k(9);
          k(9)  k(10) k(2)  k(6)  k(11) k(5);
          k(12) k(13) k(10) k(11) k(6)  k(4);
          k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
    K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
          k(11) k(14) k(11) k(12) k(9)  k(8);
          k(11) k(11) k(14) k(12) k(8)  k(9);
          k(13) k(12) k(12) k(14) k(7)  k(7);
          k(10) k(9)  k(8)  k(7)  k(14) k(11);
          k(10) k(8)  k(9)  k(7)  k(11) k(14)];
    K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
          k(2) k(1)  k(8)  k(4) k(6)  k(11);
          k(8) k(8)  k(1)  k(5) k(11) k(6);
          k(3) k(4)  k(5)  k(1) k(8)  k(2);
          k(5) k(6)  k(11) k(8) k(1)  k(8);
          k(4) k(11) k(6)  k(2) k(8)  k(1)];
    K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
          k(11) k(14) k(7)  k(12) k(9)  k(2);
          k(7)  k(7)  k(14) k(10) k(2)  k(9);
          k(13) k(12) k(10) k(14) k(7)  k(11);
          k(10) k(9)  k(2)  k(7)  k(14) k(7);
          k(12) k(2)  k(9)  k(11) k(7)  k(14)];
    KE = (1/((nu+1)*(1-2*nu))) * ...
        [K1,  K2,  K3,  K4;
         K2', K5,  K6,  K3';
         K3', K6,  K5', K2';
         K4,  K3,  K2,  K1'];
end
%% Display 3D Topology Layout Isometric View
function display_3D(rho)
[nely,nelx,nelz] = size(rho);
hx = 1; hy = 1; hz = 1;  % Element size
face = [1 2 3 4; ...
        2 6 7 3; ...
        4 3 7 8; ...
        1 5 8 4; ...
        1 2 6 5; ...
        5 6 7 8];           % Cubic face
set(gcf,'Name','Volume Minimization with Compliance Constraint','NumberTitle','off');
%% Loop Over All Elements In The 3D Grid
for k = 1:nelz
    for i = 1:nelx
        for j = 1:nely
            if rho(j,i,k) > 0.5
                x = (i-1)*hx; y = nely*hy - (j-1)*hy; z = (k-1)*hz;
                vert = [x y z; ...
                        x y-hy z; ...
                        x+hx y-hy z; ...
                        x+hx y z; ...
                        x y z+hz; ...
                        x y-hy z+hz; ...
                        x+hx y-hy z+hz; ...
                        x+hx y z+hz];
                vert(:,[2 3]) = vert(:,[3 2]); vert(:,2) = -vert(:,2);  % Reorient axes
                patch('Faces',face,'Vertices',vert, ...
                      'FaceColor',[0.2+0.8*(1-rho(j,i,k)), ...
                                   0.2+0.8*(1-rho(j,i,k)), ...
                                   0.2+0.8*(1-rho(j,i,k))]);
                hold on;
            end
        end
    end
end
%% Finalize View Settings
axis equal;          % Equal scaling on all axes
axis tight;          % Tighten axis limits to data
axis off;            % Hide axes for cleaner visualization
box on;              % Show bounding box
view([30, 30]);      % Isometric viewing angle
pause(1e-6);         % Force graphics refresh (useful during animation)
end