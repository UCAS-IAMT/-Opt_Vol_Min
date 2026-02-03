% ---------------------------------------------------------------------------------------
% To run the code use: opt_vol_min2d(200, 100, 1.9552, 3.0, 2)
% 2D Topology Optimization: Minimize Volume under Compliance Constraint
% This code is written by Md Zakirul et. al. for educational purpose only
% ---------------------------------------------------------------------------------------
function opt_vol_min2d(nelx,nely,Ctarget,penal,rmin)
%% Initialize User Configuration
maxIter = 500;    % Maximum iterations number
tol = 0.01;      % Convergence tolerance
show_Plot = 1;  % Enable real-time visualization
%% Define Material Properties
E0 = 1;           % Young's modulus (solid)
Emin = 1e-9;      % Young's modulus (void)
nu = 0.3;         % Poisson's ratio
%% Define Degrees of Freedom (Loads)
[il,jl] = meshgrid(nelx, 0);        % Coordinates
loadnid = il*(nely+1)+(nely+1-jl);  % Convert node IDs
loaddof = 2*loadnid(:);             % Extract Y-direction DOFs
%% Define Fixed Degrees of Freedom (Supports)
[iif,jf] = meshgrid(0, 0:nely);                     % Node coordinates
fixednid = iif*(nely+1)+(nely+1-jf);                % Node IDs
fixeddof = [2*fixednid(:); 2*fixednid(:)-1];        % X, Y DOFs for each fixed node
%% Initialize Finite Element Model
nele = nelx*nely;            % Total number of elements
ndof = 2*(nely+1)*(nelx+1);  % Total degrees of freedom
load_dof_x = 2*loadnid(:) - 1;   % x-DOFs
load_dof_y = 2*loadnid(:);       % y-DOFs
F = sparse([load_dof_x; load_dof_y], 1, [-0.5; -1], ndof, 1);  % External force vector
U = zeros(ndof,1);                    % Initialize displacement vector
freedofs = setdiff(1:ndof,fixeddof);  % Identify free (unconstrained) DOFs
KE = lk_Q4(nu);          % Element stiffness matrix for an a 4-node quad element
%% Construct Global Element Connectivity
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);    % Create 2D grid of node numbers
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1); % Extract the bottom-left node index
edofVec = 2*nodeids(:)+1; % Global DOF index of the x-component
%% Expand to All Elements
edofMat = repmat(edofVec,1,8)+ ...
    repmat([0 1 2*nely+[2 3 0 1] -2 -1],nele,1);
%% Prepare Sparse Matrix Indexing for Global Stiffness Assembly
iK = reshape(kron(edofMat,ones(8,1))',8*8*nele,1);
jK = reshape(kron(edofMat,ones(1,8))',8*8*nele,1);
%% Sensitivity Filter
iH = ones(nele*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% Compute Initial Compliance C0
x = ones(nely,nelx);  % Start with full material
xPhys = x;
sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),8*8*nele,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx]);
C0 = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
Cmax = C0 * Ctarget;  % Maximum target compliance
fprintf('Initial compliance C0: %.4f\n', C0);
fprintf('Target compliance: %.4f (%.2fx C0)\n', Cmax, Ctarget);

%% Record Initial State (Iteration 0)
C_ratio_hist = zeros(maxIter + 1, 1);  % Compliance history (incl. iteration 0)
V_hist = zeros(maxIter + 1, 1);        % Volume fraction history (incl. iteration 0)
C_ratio_hist(1) = 1.0;                 % C/C0 = 1 for solid initial design
V_hist(1) = 1.0;                       % Initial volume fraction = 1 (fully dense)

%% Optimization Initialization
x = ones(nely,nelx) * 0.5;  % Start with 50% density
xPhys = x; 
loop = 0; 
change = 1;
%% Start Iteration Loop
while change > tol && loop < maxIter
    loop = loop+1;
    %% FE-Analyisis
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nele,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);
    %% Compute Compliance & Strain Energy for Sensitivity Analysis
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),[nely,nelx]);
    c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
    %% Sensitivities
    dv = ones(nely,nelx);                        % Volume sensitivity
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;  % Compliance sensitivity
    %% Filtering of Sensitivites
    dv(:) = H*(dv(:)./Hs);
    dc(:) = H*(dc(:)./Hs);
    %% Density-Based Design Update with Compliance Constraint
    l1 = 0; l2 = 1e9; move = 0.1;
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        %% Update Rule: x_new = x * sqrt(-lambda*dc/dv)
        xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-lmid*dc./dv)))));
        xPhys(:) = (H*xnew(:))./Hs;
        %% Check Compliance Constraint
        sK_temp = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nele,1);
        K_temp = sparse(iK,jK,sK_temp); K_temp = (K_temp+K_temp')/2;
        U_temp = zeros(ndof,1);
        U_temp(freedofs,:) = K_temp(freedofs,freedofs)\F(freedofs,:);
        ce_temp = reshape(sum((U_temp(edofMat)*KE).*U_temp(edofMat),2),[nely,nelx]);
        c_temp = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce_temp));
        %% Adjust Lagrange Multiplier Based on Constraint
        if c_temp > Cmax, l1 = lmid; else l2 = lmid; end
    end
    change = max(abs(xnew(:)-x(:)));
    x = xnew;
    %% Compute Volume Fraction and Compliance Ratio
    vol = mean(xPhys(:));
    c_ratio = c/C0;
    %% Record History
    C_ratio_hist(loop + 1) = c_ratio;
    V_hist(loop + 1) = vol;
    %% Print Results
    fprintf(' It.:%5i Vol.:%7.3f Comp.:%11.4f C/C0:%.3f (target:%.2f) ch.:%.4f\n',...
            loop, vol, c, c_ratio, Ctarget, change);
    %% Plot Densities
    if show_Plot
        clf;
        display_2D(xPhys);
    end
end
%% Final Visualization
clf; display_2D(xPhys);
fprintf('\nOptimization completed:\n');
fprintf('Final volume fraction: %.3f\n', mean(xPhys(:)));
fprintf('Final compliance: %.4f\n', c);
fprintf('Compliance ratio C/C0: %.3f (target: %.2f)\n', c/C0, Ctarget);
%% Historical Plot
iter = 0:loop;
figure; clf;
set(gcf,'Color','w','Units','centimeters','Position',[3 3 18 12]);
yyaxis left
h1 = plot(iter, C_ratio_hist(1:loop+1), '-', 'Color', [0.00 0.45 0.74], 'MarkerFaceColor', [0.00 0.45 0.74], ...
    'MarkerEdgeColor', [0.00 0.45 0.74], 'MarkerSize', 6, 'LineWidth', 2);
hold on;
yline_left = yline(Ctarget, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 2);
set(gca,'YScale','log');
ylabel('Compliance Ratio', 'FontSize', 30, 'FontWeight','bold');
set(gca, 'FontSize', 30, 'LineWidth', 2, 'TickLength', [0.02 0.02], 'YColor', [0.00 0.45 0.74], 'XColor', 'k');
yyaxis right
h2 = plot(iter, V_hist(1:loop+1), '-', ...
    'Color', [0.49 0.18 0.56], 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.49 0.18 0.56], ...
    'MarkerSize', 6, 'LineWidth', 2);

ylabel('Volume Fraction', 'FontSize', 30, 'FontWeight','bold');

set(gca, 'FontSize', 30, 'LineWidth', 2, 'TickLength', [0.02 0.02], 'YColor', [0.49 0.18 0.56]);
xlabel('Iteration', 'FontSize', 30, 'FontWeight','bold');
legend([h1, h2, yline_left], ...
    {'Compliance Ratio $C/C_0$','Volume Fraction','Target Compliance'}, 'Location','best', 'FontSize', 24, 'Box','on','Interpreter','latex');
grid on;
yyaxis left
set(gca, 'XMinorGrid','off', 'YMinorGrid','off', 'XMinorTick','off', 'YMinorTick','off', 'Box','on', 'Layer','top');
yyaxis right
set(gca, 'XMinorGrid','off', 'YMinorGrid','off', 'XMinorTick','off', 'YMinorTick','off', 'Box','on', 'Layer','top');
drawnow;
end
%% 4-node Quadrilateral Element Stiffness Matrix
function [KE] = lk_Q4(nu)
k = [1/2-nu/6   1/8+nu/8   -1/4-nu/12 -1/8+3*nu/8 ...
     -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
KE = 1/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
                  k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
                  k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
                  k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
                  k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
                  k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
                  k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
                  k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
end
%% Display 2D Topology Layout
function display_2D(rho)
%% Apply Threshold for Sharper Visualization
rho_display = rho;
rho_display(rho < 0.5) = 0;  % Set low densities to void
rho_display(rho >= 0.5) = 1; % Set high densities to solid
colormap([1 1 1; 0 0 1]); imagesc(rho_display); caxis([0 1]); axis equal; axis off; drawnow;
end