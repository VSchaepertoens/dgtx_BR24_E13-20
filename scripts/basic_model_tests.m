% Basic stats
if isfield(model, 'description') && ischar(model.description)
    fprintf('Model ID/description: %s\n', model.description);
else
    fprintf('Model ID/description: no description field found\n');
end
fprintf('Reactions: %d\nMets: %d\nGenes: %d\n', length(model.rxns), length(model.mets), length(model.genes));

% Model stoichiometry dimensions
disp('Stoichiometry matrix size:');
disp(size(model.S));

% % --- Check mass & charge balance --- not giving any info for now
% fprintf('\nChecking mass and charge balance...\n');
% [unbalancedRxns, balancedRxns, elementalImbalance, hasFormula] = checkMassChargeBalance(model);
% 
% % Ensure unbalancedRxns is full logical
% unbalancedRxnsFull = full(unbalancedRxns);   % convert sparse to full logical vector
% 
% % Count number of unbalanced reactions
% numUnbalanced = sum(unbalancedRxnsFull);     % scalar now
% fprintf('Number of unbalanced reactions: %d\n', numUnbalanced);
% 
% % Get numeric indices for display
% rxnIdxList = find(unbalancedRxnsFull);
% 
% if ~isempty(rxnIdxList)
%     disp('Example of unbalanced reactions:');
%     nShow = min(5, numel(rxnIdxList));
%     for i = 1:nShow
%         rxnIdx = rxnIdxList(i);           % numeric index of unbalanced reaction
%         fprintf('%d: %s\n', rxnIdx, model.rxns{rxnIdx});
%     end
% else
%     disp('No unbalanced reactions found.');
% end

% --- Optimize model with default objective --- 
FBAsol = optimizeCbModel(model);
fprintf('FBA objective value: %g\n', FBAsol.f);

% Show top fluxes (absolute value)
[~, idx] = sort(abs(FBAsol.x), 'descend');
topN = 10;
fprintf('Top %d fluxes:\n', topN);
for i = 1:topN
    fprintf('%s: %g\n', model.rxns{idx(i)}, FBAsol.x(idx(i)));
end

% --- FVA --- 
% Compute FVA at 100% of optimum (can take a few seconds)
[minFlux, maxFlux] = fluxVariability(model, 100, 'max', model.rxns);

% Show example of first 5 reactions
fprintf('FVA ranges for first 5 reactions:\n');
for i = 1:5
    fprintf('%s: [%g, %g]\n', model.rxns{i}, minFlux(i), maxFlux(i));
end