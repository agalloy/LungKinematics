%% Initialize Matlab
clear
clc

%% Data parameters
% Important directories
disp_dir = 'Y:\Documents\BioMOST_Research\Lung_FE\FEBio\FEBio_Runs\TLCtoFRC_PenaltyStep';
mesh_dir = 'Y:\Documents\BioMOST_Research\Lung_FE\FEBio\Meshes_v3';

% Important file patterns
mesh_pattern = '${SUBJECT}_${SIDE}Lung_Lobes_Mesh_v3.mat';
disp_pattern = '${SUBJECT}_${SIDE}Lung_Lobes_f${fc}_ndata.txt';

% Subjects to analyze
subjects = "all"; % To analyze all subjects found in disp_dir use "all"
exclude = ""; % May be empty

% Set up the models to compare
% Model parameters (1 x P string array)
model_params = ["${SIDE}","${fc}"];
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

% Entries for each model in legend
legend_entries = ["Non-Sliding"
                  "Sliding"];

% Region names (Same order as disp data, left lung first then right)
region_names = ["Left Lung","LLL","LUL","Right Lung","RLL","RML","RUL"];
% Which regions to plot? If empty will default to all. Use the index of the above regions.
to_plot = [];
              
% Use the same mesh for models sharing a subject?
same_mesh = false;

%% Additional Analyses
plot_ADIdiff = true;
    % Scatter plot colors
    sc = {[1,0,0],[1,0,0]};
    % Scatter plot shapes
    sm = {'x','o'};
    % Scatter plot tick offsets and spacing
    tick_width = 0.8;
    tick_spacing = 1.5;

stats_ADIdiff = true;

maxShear_analysis = true;
    % Max shear voxel grid spacing
    spacing = [5,5,5];

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

%% Loop through each model and calculate mean ADI in lung and lobes
% Initialize Arrays
mean_ADI = cell(num_subjects,num_models);
mean_MaxShear = cell(num_subjects,num_models);
ref_vol = cell(num_subjects,num_models);
def_vol = cell(num_subjects,num_models);
vol_change = cell(num_subjects,num_models);

tic
fprintf('\nLoading deformation data...\n')
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
        
        % Using new reference configuration and displacements calculate ADI
        [lambda,vol] = Disp2Stretches( NodeArray_Ref, ElementArray(eID~=1,:), DispArray );
        vol = vol/1000; % convert vol to mL
        ADI = sqrt( ((lambda(:,1)-lambda(:,2))./lambda(:,2)).^2 + ((lambda(:,2)-lambda(:,3))./lambda(:,3)).^2 );
        J = lambda(:,1) .* lambda(:,2) .* lambda(:,3);
        
        % Get mean ADI in whole lung and each lobe
        eID_solid = eID(eID~=1);
        mean_ADI_reg = nan(max(eID),1);
        mean_ADI_reg(1) = (vol' * ADI) / sum(vol);        
        for k = 2:max(eID)
            % Find elements in current region
            in_region = eID_solid == k;
            mean_ADI_reg(k) = (vol(in_region)' * ADI(in_region)) / sum(vol(in_region));
        end
        mean_ADI{i,j} = mean_ADI_reg;
        
        % Calculate deformed and reference lung volumes (in mL)
        ref_vol_reg = nan(max(eID),1);
        def_vol_reg = nan(max(eID),1);
        ref_vol_reg(1) = sum(vol);     
        def_vol_reg(1) = sum( vol .* J );
        for k = 2:max(eID)
            % Find elements in current region
            in_region = eID_solid == k;
            % Reference volume
            ref_vol_reg(k) = sum(vol(in_region));
            def_vol_reg(k) = sum( vol(in_region) .* J(in_region) );
        end
        ref_vol{i,j} = ref_vol_reg;
        def_vol{i,j} = def_vol_reg;
        vol_change{i,j} = ref_vol{i,j} - def_vol{i,j};
        
        % Get mean max shear along the fissures
        if maxShear_analysis
            [mean_MaxShear{i,j}, ~] = ShearContourV2( NodeArray_Ref, ElementArray, eID, DispArray, spacing );
        end
    end
end
fprintf('Deformation data loaded!\n')
toc

%% Combine left and right lung data if necessary
num_regions = numel(region_names);
if ismember("${SIDE}",model_params)
    % Combine paired models
    mean_ADI_comb = cellfun(@(x,y) [x;y], mean_ADI(:,1:2:end), mean_ADI(:,2:2:end), 'UniformOutput', false );
    num_comb_models = num_models / 2;
else
    mean_ADI_comb = mean_ADI;
    num_comb_models = num_models;
end

%% Plot results
if plot_ADIdiff
    offsets = linspace( -tick_width/2, tick_width/2, num_comb_models );

    if isempty(to_plot)
        plot_regions = 1:num_regions;
    else
        plot_regions = to_plot;
    end
    mean_ADI_plot = cellfun( @(x) x(plot_regions), mean_ADI_comb, 'UniformOutput', false);
    num_pregions = numel(plot_regions);

    figure()
    hold on
    % Generate scatter plot
    for i = 1:num_comb_models
        offset = offsets(i);
        scatter_y = cell2mat( mean_ADI_plot(:,i) );
        scatter_x = repmat( (1:num_pregions)', 1, num_subjects )*tick_spacing + offset;
        scatter_x = reshape( scatter_x, [], 1 );
        plot( scatter_x, scatter_y, 'o', ...
             'Marker', sm{i}, 'Color', sc{i}, 'MarkerSize', 10, 'LineWidth', 3 )
    end
    % Generate lines
    for i = 1:num_pregions
        line_x = i*ones(num_subjects,num_comb_models)*tick_spacing + offsets;
        line_y = cellfun(@(x) x(i), mean_ADI_plot);
        line(line_x',line_y','Color','k','LineStyle','-')
    end

    %Format Axes
    box on
    xlim([0,tick_spacing*(num_pregions+1)])
    %ylim([0.04 0.18])
    ax = gca;
    ax.XTick = tick_spacing:tick_spacing:tick_spacing*num_pregions;
    %ax.YTick = 0.0:0.05:0.50;
    ax.XTickLabels = region_names(plot_regions);
    ax.FontSize = 13;
    ylabel('Mean ADI','FontSize',20)
    %Format Legend
    legend(legend_entries)
    hold off
end

%% Paired statistical comparison of the 1st and 2nd (combined) model
if stats_ADIdiff
    p = nan(1,num_regions);
    % Get lists of values for each model for paired comparison
    stat_data1 = reshape( cell2mat(mean_ADI_comb(:,1)), num_regions, num_subjects )';
    stat_data2 = reshape( cell2mat(mean_ADI_comb(:,2)), num_regions, num_subjects )';

    % Perform paired Sum of Signed Rank test
    for i = 1:num_regions
        p(i) = signrank( stat_data1(:,i), stat_data2(:,i), 'tail', 'both');
    end
    sig = p < 0.05;

    % Compute absolute and percent differences 
    diff_data = stat_data1 - stat_data2;
    pdiff_data = (stat_data1 - stat_data2) ./ stat_data2 * 100;

    % Display results in table
    out_table = table( 'Size', [4,num_regions], 'VariableTypes', repmat("double",1,num_regions) );
    out_table.Properties.VariableNames = region_names;
    out_table.Properties.RowNames = {'p-values', 'Significance', 'Median Absolute Difference', 'Median Percent Difference'};
    out_table{1,:} = p;
    out_table{2,:} = sig;
    out_table{3,:} = median(diff_data);
    out_table{4,:} = median(pdiff_data);
    fprintf('\nMean ADI: %s vs. %s\n',legend_entries(2),legend_entries(1))
    disp(out_table)
end