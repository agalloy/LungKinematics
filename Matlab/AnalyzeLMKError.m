%% Initialize Matlab
clear
clc

%% Data parameters
% Important directories
disp_dir = 'Y:\Documents\BioMOST_Research\Lung_FE\FEBio\FEBio_Runs\TLCtoFRC_PenaltyStep';
mesh_dir = 'Y:\Documents\BioMOST_Research\Lung_FE\FEBio\Meshes_v3';
lmk_dir = 'Z:\Lung\Dissertation\';

% Important file patterns
mesh_pattern = '${SUBJECT}_${SIDE}Lung_Lobes_Mesh_v3.mat';
disp_pattern = '${SUBJECT}_${SIDE}Lung_Lobes_f${FC}_ndata.txt';
lmk_pattern = '${SUBJECT}_${SIDE}_landmarks-disp.txt';

% Subjects to analyze
subjects = "all"; % To analyze all subjects found in disp_dir use "all"
exclude = ""; % May be empty

% Set up the models to compare
% Model parameters (1 x P string array)
model_params = ["${SIDE}","${FC}"];
% Values to set those parameters to (M x P string array)
% M is the number of models, P is the number of parameters
% If ${SIDE} is a parameter, then consecutive models will be combined
% after loading data assuming the first is "Left" and the second is "Right"
model_values = ["Left","1.5"
                "Right","1.5"
                "Left","0"
                "Right","0"];
            
% Initial and final step for each model
t_int = {[0,1],[0,1],[0,1],[0,1]};

% Entries for each model in output table
names = {'Left: non-sliding', 'Right: non-sliding', 'Left: sliding', 'Right: sliding'};

% Region names (Same order as disp data, left lung first then right)
%region_names = ["LL","LLL","LUL"];
region_names = ["LL","LLL","LUL","RL","RLL","RML","RUL"];
              
% Use the same mesh for models sharing a subject?
same_mesh = false;

%% Create a subject list
if strcmp(subjects(1),"all")
    subject_pattern = replace( disp_pattern, ["${SUBJECT}",model_params], repmat("*",1,1+numel(model_params)) );
    disp_dir_files = dir(fullfile(disp_dir,subject_pattern));
    subject_list = strings(length(disp_dir_files),1);
    % Get all the subject names in folder
    for i = 1:length(disp_dir_files)
        subject = disp_dir_files(i).name;
        subject = extractBefore(subject,'_');
        subject_list(i) = subject;
    end
    % Remove redundancies
    subject_list = unique(subject_list);
else
    subject_list = subjects';
end

% Remove the entries in subjects that are in the exclude list
subject_list = subject_list( ~ismember(subject_list,exclude) );  

num_models = size(model_values,1);
num_subjects = size(subject_list,1);

%% Load mesh data for each subject
if same_mesh
    num_mesh = 1;
else
    num_mesh = num_models;
end
NodeCells = cell(num_subjects,num_mesh);
ElementCells = cell(num_subjects,num_mesh);
eIDCells = cell(num_subjects,num_mesh);

tic
fprintf('\nLoading mesh data...\n')
for i = 1:num_subjects
    for j = 1:num_mesh
        subject = char(subject_list(i));

        % Get mesh data
        mesh_name = replace( mesh_pattern, ["${SUBJECT}",model_params], [subject,model_values(j,:)] );
        mesh_file = fullfile( mesh_dir, mesh_name );
        load( mesh_file, "NodeArray", "ElementArray", "nID", "eID" );

        % Fill cell arrays
        NodeCells{i,j} = NodeArray;
        ElementCells{i,j} = ElementArray;
        eIDCells{i,j} = eID;
    end
end
fprintf('Mesh data loaded!\n')
toc

%% Loop through each model and get average landmark error and displacement
% Initialize Arrays
u_lmk_avg = cell(num_subjects,num_models);
e_lmk_avg = cell(num_subjects,num_models);

tic
fprintf('\nLoading landmark data...\n')
for i = 1:num_subjects
    subject = char(subject_list(i));
    
    for j = 1:num_models
        % Get subject's mesh info
        if same_mesh
            NodeArray = NodeCells{i};
            ElementArray = ElementCells{i};
            eID = eIDCells{i};
        else
            NodeArray = NodeCells{i,j};
            ElementArray = ElementCells{i,j};
            eID = eIDCells{i,j};
        end
    
        % Get disp data
        disp_name = replace( disp_pattern, ["${SUBJECT}",model_params], [subject,model_values(j,:)] );
        disp_file = fullfile(disp_dir,disp_name);
        t_get = t_int{j};
        try
            [ n_data, ~ ] = GetNodeData( disp_file, t_get );
        catch GND_exc
            fprintf('\nFailed to read disp_file. Subject: %s Model: %i \n',subject,j)
            continue
        end
        
        % Set reference geometry to the desired time step
        if t_get(1) > 0
            DispOffset = n_data{1}(:,2:4);
            NodeArray_Ref = NodeArray + DispOffset;
            % Offset final displacements accordingly
            DispArray = n_data{2}(:,2:4) - DispOffset;
        else
            NodeArray_Ref = NodeArray;
            DispArray = n_data{end}(:,2:4);
        end 
        
        % Get landmark error for the subject if available
        lmk_name = replace( lmk_pattern, ["${SUBJECT}",model_params], [subject,model_values(j,:)] );
        lmk_file = fullfile(lmk_dir,lmk_name);
        try
            lmk_struct = GetLMKError(NodeArray_Ref, ElementArray(eID~=1,:), DispArray, lmk_file);
        catch GND_exc
            fprintf('\nFailed to read lmk_file. Subject: %s Model: %i \n',subject,j)
            continue
        end            
        u_lmk = lmk_struct.u_lmk;
        u_fea = lmk_struct.u_fea;
        e_lmk = u_fea - u_lmk;
        
        % Get the mean landmark displacement and error
        u_lmk_mag = vecnorm( u_lmk, 2, 2 );
        u_lmk_avg{i,j} = mean( u_lmk_mag );
        e_lmk_mag = vecnorm( e_lmk, 2, 2 );
        e_lmk_avg{i,j} = mean( e_lmk_mag );        
    end
end
fprintf('Landmark data loaded!\n')
toc

%% Display results
disp('Mean landmark displacements')
u_table = [cellstr(subject_list), u_lmk_avg];
u_table = [['Subjects',names];u_table];
disp(u_table)
disp('Mean landmark errors')
e_table = [cellstr(subject_list), e_lmk_avg];
e_table = [['Subjects',names];e_table];
disp(e_table)