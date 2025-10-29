%% ==========================================
% Initialize COBRA Toolbox and load CHO model
% ===========================================

% --- Adjust these paths to your setup ---
projectDir = 'C:\Users\b1095820\OneDrive - Universität Salzburg\Documents\DGTX\BR24_E13-20';  % your project folder
cobraDir   = 'C:\Users\b1095820\Documents\MATLAB\cobratoolbox';  % path to COBRA Toolbox

% --- Add COBRA Toolbox to MATLAB path ---
addpath(genpath(cobraDir));

% --- Initialize COBRA Toolbox ---
initCobraToolbox(false);   % false = do not update from GitHub

% --- Change to your project directory ---
cd(projectDir);

% --- Define model path ---
modelPath = fullfile(projectDir, 'cho_gems', 'iCHOv1_final.xml');  % adjust file name if needed

% --- Load model ---
model = readCbModel(modelPath);

% --- Confirm successful loading ---
disp('✅ COBRA Toolbox initialized.');
disp(['✅ Model loaded: ', model.description]);
disp(['Number of reactions: ', num2str(length(model.rxns))]);
disp(['Number of metabolites: ', num2str(length(model.mets))]);
disp(['Number of genes: ', num2str(length(model.genes))]);
