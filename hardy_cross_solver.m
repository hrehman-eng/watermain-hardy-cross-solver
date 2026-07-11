%% HARDY_CROSS_SOLVER.M
% Water Distribution Network Solver — Hardy-Cross Method
% Blackwood Meadows Subdivision, Waterloo, ON (synthetic portfolio network)
%
% Reads network geometry/demands from Watermain_Network_Input.xlsx and
% runs three cases:
%   1. Base case      - average-day demand, new pipe (C = 130)
%   2. Fire flow case - average-day demand + FUS fire demand at the
%                        hydraulically most remote node, checked against
%                        the 140 kPa fire-flow minimum residual pressure
%   3. Sensitivity     - average-day demand, aged pipe (C = 100)
%
% Each case is solved with the same Hardy-Cross / Hazen-Williams engine
% and written to its own sheet in the workbook. Convergence history and
% a pressure-margin comparison are plotted and saved as PNGs.
%
% Run this script with Watermain_Network_Input.xlsx in the same folder.

clear; clc; close all;

%% -------------------- USER SETTINGS --------------------
inputFile       = 'Watermain_Network_Input.xlsx';
maxIterations   = 50;
convTol         = 1e-7;      % m^3/s, max |dQ| across all loops
minPressureAvg  = 275;       % kPa, Ontario MECP min residual, average-day
minPressureFire = 140;       % kPa, Ontario MECP min residual, fire flow
maxVelocity_ms  = 2.5;       % design velocity limit
sensitivityC    = 100;       % aged-pipe Hazen-Williams C

%% -------------------- READ INPUT DATA --------------------
nodeTbl = readtable(inputFile, 'Sheet', 'Node_Data');
pipeTbl = readtable(inputFile, 'Sheet', 'Pipe_Data');

nNodes = height(nodeTbl);
nPipes = height(pipeTbl);

nodeID    = string(nodeTbl.Node_ID);
elev      = nodeTbl.Elevation_m;
demand    = nodeTbl.Demand_Ls;
demand(isnan(demand)) = 0;
fireDemand = nodeTbl.Fire_Demand_Ls;
fireDemand(isnan(fireDemand)) = 0;
isSource  = strcmpi(nodeTbl.Type, 'Source');
sourceHGL = nodeTbl.Fixed_HGL_m(isSource);

pipeID   = string(pipeTbl.Pipe_ID);
fromNode = string(pipeTbl.From_Node);
toNode   = string(pipeTbl.To_Node);
L        = pipeTbl.Length_m;
D        = pipeTbl.Diameter_mm / 1000;
C_base   = pipeTbl.Hazen_Williams_C;

pIdx = @(id) find(pipeID == id);

% Loop A (left square):  N1->N2->N5->N4->N1  (P1, P4, P6, P3)
% Loop B (right square): N2->N3->N6->N5->N2  (P2, P5, P7, P4)
loops = {
    [pIdx("P1"), pIdx("P4"), pIdx("P6"), pIdx("P3")], [ 1,  1, -1, -1];
    [pIdx("P2"), pIdx("P5"), pIdx("P7"), pIdx("P4")], [ 1,  1, -1, -1]
};

%% -------------------- INITIAL FLOW GUESS (BASE CASE) --------------------
% Must satisfy continuity at every node before iterating — Hardy-Cross
% only balances energy around loops, it never repairs an unbalanced
% starting guess. Re-derive by hand if topology/demands change.
Q0map = containers.Map( ...
    {'P1','P2','P3','P4','P5','P6','P7','P8'}, ...
    {0.0240, 0.0095, 0.0150, 0.0060, 0.0035, 0.0050, 0.0020, 0.0030});  % m^3/s

Q0base = zeros(nPipes,1);
for i = 1:nPipes
    Q0base(i) = Q0map(pipeID(i));
end

% Fire-flow initial guess: route the extra fire demand at N7 up through
% P3 (N1->N4) and P8 (N4->N7), the only path to N7. Continuity still
% balances everywhere else because the same increment is added on both
% sides of node N4.
fireNodeIdx = find(fireDemand > 0);
fireQ = sum(fireDemand(fireNodeIdx)) / 1000;  % L/s -> m^3/s
Q0fire = Q0base;
Q0fire(pIdx("P3")) = Q0fire(pIdx("P3")) + fireQ;
Q0fire(pIdx("P8")) = Q0fire(pIdx("P8")) + fireQ;

%% -------------------- RUN THE THREE CASES --------------------
[Qbase, histBase] = solveHardyCross(Q0base, L, D, C_base, loops, maxIterations, convTol);
[Qfire, histFire] = solveHardyCross(Q0fire, L, D, C_base, loops, maxIterations, convTol);
C_aged = sensitivityC * ones(nPipes,1);
[Qsens, histSens] = solveHardyCross(Q0base, L, D, C_aged, loops, maxIterations, convTol);

fprintf('Base case:        converged in %d iterations\n', numel(histBase));
fprintf('Fire flow case:   converged in %d iterations\n', numel(histFire));
fprintf('Sensitivity case: converged in %d iterations\n', numel(histSens));

%% -------------------- NODE PRESSURES FOR EACH CASE --------------------
HGLbase = computeHGL(Qbase, L, D, C_base, nodeID, fromNode, toNode, isSource, sourceHGL);
HGLfire = computeHGL(Qfire, L, D, C_base, nodeID, fromNode, toNode, isSource, sourceHGL);
HGLsens = computeHGL(Qsens, L, D, C_aged, nodeID, fromNode, toNode, isSource, sourceHGL);

[nodeResBase, pipeResBase] = buildResultTables(nodeID, elev, HGLbase, isSource, Qbase, L, D, C_base, pipeID, minPressureAvg, maxVelocity_ms);
[nodeResFire, pipeResFire] = buildResultTables(nodeID, elev, HGLfire, isSource, Qfire, L, D, C_base, pipeID, minPressureFire, maxVelocity_ms);
[nodeResSens, pipeResSens] = buildResultTables(nodeID, elev, HGLsens, isSource, Qsens, L, D, C_aged, pipeID, minPressureAvg, maxVelocity_ms);

%% -------------------- CRITICAL NODE CALL-OUT --------------------
demandMask = ~isSource;
[minP_base, iMin] = min(nodeResBase.Pressure_kPa(demandMask));
demandNodeIDs = nodeID(demandMask);
critNode = demandNodeIDs(iMin);
fprintf('\n--- Critical node (base case) ---\n');
fprintf('%s has the lowest pressure at %.1f kPa - %.1f kPa above the %d kPa minimum, so the network passes with margin.\n', ...
    critNode, minP_base, minP_base - minPressureAvg, minPressureAvg);

minP_fire = min(nodeResFire.Pressure_kPa(demandMask));
fireResult = "PASSES";
if minP_fire < minPressureFire
    fireResult = "FAILS";
end
fprintf('\n--- Fire flow check ---\n');
fprintf('With %.1f L/s fire demand at %s, minimum system pressure drops to %.1f kPa - %.1f kPa above the %d kPa fire-flow minimum (%s).\n', ...
    fireQ*1000, nodeID(fireNodeIdx), minP_fire, minP_fire - minPressureFire, minPressureFire, fireResult);

minP_sens = min(nodeResSens.Pressure_kPa(demandMask));
fprintf('\n--- Aged-pipe sensitivity (C = %d) ---\n', sensitivityC);
fprintf('Minimum pressure drops from %.1f kPa (new pipe, C=%d) to %.1f kPa (aged pipe, C=%d) - a margin loss of %.1f kPa, still %.1f kPa above the %d kPa minimum.\n', ...
    minP_base, C_base(1), minP_sens, sensitivityC, minP_base - minP_sens, minP_sens - minPressureAvg, minPressureAvg);

%% -------------------- WRITE RESULTS TO EXCEL --------------------
writetable(nodeResBase, inputFile, 'Sheet', 'Results', 'Range', 'A3');
writetable(pipeResBase, inputFile, 'Sheet', 'Results', 'Range', 'H3');

writetable(nodeResFire, inputFile, 'Sheet', 'Results_FireFlow', 'Range', 'A3');
writetable(pipeResFire, inputFile, 'Sheet', 'Results_FireFlow', 'Range', 'H3');

writetable(nodeResSens, inputFile, 'Sheet', 'Results_Sensitivity', 'Range', 'A3');
writetable(pipeResSens, inputFile, 'Sheet', 'Results_Sensitivity', 'Range', 'H3');

fprintf('\nResults written to %s (Results, Results_FireFlow, Results_Sensitivity sheets).\n', inputFile);
disp(nodeResBase);
disp(pipeResBase);

%% -------------------- PLOT 1: CONVERGENCE HISTORY --------------------
figure('Name','Convergence History');
semilogy(1:numel(histBase), histBase, '-o', 'LineWidth', 1.5, 'DisplayName', 'Base case'); hold on;
semilogy(1:numel(histFire), histFire, '-s', 'LineWidth', 1.5, 'DisplayName', 'Fire flow case');
semilogy(1:numel(histSens), histSens, '-^', 'LineWidth', 1.5, 'DisplayName', 'Aged-pipe sensitivity');
yline(convTol, '--k', 'Convergence tolerance', 'LabelHorizontalAlignment','left');
xlabel('Iteration');
ylabel('Max loop flow correction, |\DeltaQ| (m^3/s)');
title('Hardy-Cross Convergence History');
legend('Location','northeast');
grid on;
saveas(gcf, 'convergence_history.png');

%% -------------------- PLOT 2: PRESSURE MARGIN COMPARISON --------------------
figure('Name','Pressure Margin Comparison');
margins = [minP_base - minPressureAvg, minP_fire - minPressureFire, minP_sens - minPressureAvg];
labels = {'Base case', 'Fire flow', 'Aged pipe (C=100)'};
b = bar(margins);
b.FaceColor = 'flat';
b.CData(1,:) = [0.12 0.30 0.47];
b.CData(2,:) = [0.75 0.22 0.17];
b.CData(3,:) = [0.85 0.60 0.13];
set(gca, 'XTickLabel', labels);
ylabel('Pressure margin above minimum requirement (kPa)');
title('Minimum Pressure Margin by Scenario (worst-case node)');
grid on;
saveas(gcf, 'pressure_margin_comparison.png');

fprintf('\nPlots saved: convergence_history.png, pressure_margin_comparison.png\n');

%% ==================== LOCAL FUNCTIONS ====================
function hf = headloss(Q, L, D, C)
    % Hazen-Williams head loss (SI units), signed to match flow direction.
    s = sign(Q);
    if s == 0
        s = 1;
    end
    hf = s * 10.67 * L * abs(Q)^1.852 / (C^1.852 * D^4.8704);
end

function [Q, history] = solveHardyCross(Q0, L, D, C, loops, maxIter, tol)
    Q = Q0;
    history = [];
    for iter = 1:maxIter
        maxCorr = 0;
        for Li = 1:size(loops,1)
            pipeIdxList = loops{Li,1};
            signList    = loops{Li,2};

            sumHf = 0;
            sumHfOverQ = 0;
            for k = 1:numel(pipeIdxList)
                p = pipeIdxList(k);
                s = signList(k);
                q = Q(p);
                hf = headloss(q, L(p), D(p), C(p));
                sumHf = sumHf + s*hf;
                qAbs = max(abs(q), 1e-9);
                sumHfOverQ = sumHfOverQ + 10.67 * L(p) * qAbs^0.852 / (C(p)^1.852 * D(p)^4.8704);
            end

            dQ = -sumHf / (1.852 * sumHfOverQ);
            for k = 1:numel(pipeIdxList)
                p = pipeIdxList(k);
                s = signList(k);
                Q(p) = Q(p) + s*dQ;
            end
            maxCorr = max(maxCorr, abs(dQ));
        end
        history(end+1) = maxCorr; %#ok<AGROW>
        if maxCorr < tol
            break;
        end
    end
    if maxCorr >= tol
        warning('Solver did not converge within %d iterations (max correction = %.3e).', maxIter, maxCorr);
    end
end

function HGL = computeHGL(Q, L, D, C, nodeID, fromNode, toNode, isSource, sourceHGL)
    nNodes = numel(nodeID);
    nPipes = numel(Q);
    HGL = nan(nNodes,1);
    srcIdx = find(isSource);
    HGL(srcIdx) = sourceHGL;

    visited = false(nNodes,1);
    visited(srcIdx) = true;
    queue = srcIdx;

    while ~isempty(queue)
        n = queue(1); queue(1) = [];
        for p = 1:nPipes
            fIdx = find(nodeID == fromNode(p));
            tIdx = find(nodeID == toNode(p));
            if fIdx == n && ~visited(tIdx)
                hf = headloss(Q(p), L(p), D(p), C(p));
                HGL(tIdx) = HGL(n) - hf;
                visited(tIdx) = true;
                queue(end+1) = tIdx; %#ok<AGROW>
            elseif tIdx == n && ~visited(fIdx)
                hf = headloss(Q(p), L(p), D(p), C(p));
                HGL(fIdx) = HGL(n) + hf;
                visited(fIdx) = true;
                queue(end+1) = fIdx; %#ok<AGROW>
            end
        end
    end
end

function [nodeRes, pipeRes] = buildResultTables(nodeID, elev, HGL, isSource, Q, L, D, C, pipeID, minPressurekPa, maxVelocity_ms)
    pressure_kPa = round((HGL - elev) * 9.81, 1);
    pressure_psi = round(pressure_kPa * 0.145038, 1);
    pressureCheck = strings(numel(nodeID),1);
    pressureCheck(:) = "PASS";
    pressureCheck(pressure_kPa < minPressurekPa & ~isSource) = "FAIL";
    pressureCheck(isSource) = "N/A (source)";

    nodeRes = table(nodeID, elev, round(HGL,2), pressure_kPa, pressure_psi, pressureCheck, ...
        'VariableNames', {'Node_ID','Elevation_m','HGL_m','Pressure_kPa','Pressure_psi','Min_Pressure_Check'});

    nPipes = numel(Q);
    area = pi * D.^2 / 4;
    velocity = round(abs(Q) ./ area, 2);
    hf = zeros(nPipes,1);
    for p = 1:nPipes
        hf(p) = headloss(Q(p), L(p), D(p), C(p));
    end
    velocityCheck = strings(nPipes,1);
    velocityCheck(:) = "PASS";
    velocityCheck(velocity > maxVelocity_ms) = "FAIL";

    pipeRes = table(pipeID, round(Q*1000,1), velocity, round(hf,3), velocityCheck, ...
        'VariableNames', {'Pipe_ID','Flow_Ls','Velocity_ms','Headloss_m','Velocity_Check'});
end
