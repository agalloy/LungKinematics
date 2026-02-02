% Script for examining ADI vs distance from the fissure and comparing
% between sliding and non-sliding models

%% Initialize Matlab
clear
clc
addpath('Z:\AdamGalloy\TecPlotTools')

%% User Inputs: File Processing
% FEBio results directory
results_dir = 'Z:\AdamGalloy\Lung FE\TecPlot\FEBioRuns\TLCtoFRC_PenaltyStep';
% Results file pattern
results_pattern = '${SUBJECT}_LeftLung_Lobes_${FC}.txt';
% Subject(s) to run
% subjects = ["H5972","H5974","H5978","H5983","H6012","H6019"];
subjects = "H5972";
% Lung to analyze (left or right)
side = 'left';
% Parameters to for each model (for each subject)
model_params = ["${FC}"];
% Values for each model (size: num_models X num_params)
model_values = ["f1.5"
                "f0"];
% Name for each model (for use in figures)
model_names = ["Non-Sliding"
               "Sliding"];

%% Plot settings
% Include additional debugging plots
debug_plots = false;

if strcmpi(side,'left')
    lobe_names = ["LLL";"LUL"];
    lineColor = {[0.28,0.07,0.39],[0.48,0.82,0.32]};  % presentation green [24,222,43]/255
elseif strcmpi(side,'right')
    lobe_names = ["RLL";"RML";"RUL"];
    lineColor = {[0.28,0.07,0.39],[0.93,0.69,0.13],[0.48,0.82,0.32]};  % presentation green [24,222,43]/255
else
    error('Pick a correct side!')
end

num_models = size(model_values,1);
num_regions = numel(lobe_names);

% Legend text
legend_template = "${LOBE}: ${MODEL}";
legend_text = repmat(legend_template, num_regions, num_models);
legend_text = strrep( legend_text, repmat("${LOBE}", num_regions, num_models), repmat(lobe_names,1,num_models) );
legend_text = strrep( legend_text, repmat("${MODEL}", num_regions, num_models), repmat(model_names',num_regions,1) );
legend_text = reshape( legend_text, [], 1 );

% Plot settings
lineType = {'-','--'};
lineWidth = 1.5;

%% Loop through each subject and each model and load results
num_subjects = size(subjects,2);

% Loop through each subject
results = cell(num_subjects,num_models);
for i = 1:num_subjects
    subject = char(subjects(i));
    % Loop through each model for this subject and collect results data
    for j = 1:num_models
        % Get the appropriate file name from the pattern
        results_name = replace(results_pattern,["${SUBJECT}",model_params],[subject,model_values(j)]);
        results_file = fullfile(results_dir,results_name);        
        % Read the data for the current model
        results{i,j} = ReadTecPlotData(results_file);
    end
end

%% Plot ADI vs Fissure distance for each model and compare differences
fissure_area = nan(num_subjects,num_regions);
lobe_volume = nan(num_subjects,num_regions);

for i = 1:num_subjects
    subject = char(subjects(i));
    % Get each lobar surface
    NodePos = results{i,1}.NodeArray(2:end);
    ElementArray = results{i,1}.ElementArray(2:end);
    FaceArray = cellfun( @(x) FESurface(x), results{i,1}.ElementArray(2:end), 'UniformOutput', false );
    num_lobes = numel(NodePos);
    
    % Get each pair of surfaces
    % WARNING only works for sets of 2 or 3!!!
    lobe_pair1 = repmat(1:num_lobes,num_lobes,1);
    lobe_pair1(1 : num_lobes+1 : end) = [];
    lobe_pair1 = lobe_pair1';
    lobe_pair2 = repmat(1:num_lobes,num_lobes,1)';
    lobe_pair2(1 : num_lobes+1 : end) = [];
    lobe_pair2 = lobe_pair2';
    lobe_pairs = [lobe_pair1,lobe_pair2];
    num_pairs = size(lobe_pairs,1);
    
    % Find the fissure surfaces for each lobe
    cppStructL = arrayfun( @(x,y) ClosestPointTriSurfV2(FaceArray{x},NodePos{x},NodePos{y}), lobe_pair2, lobe_pair1, 'UniformOutput', false );
    f_nodes = cellfun( @(x) x.dist < 1 | x.inside, cppStructL, 'UniformOutput', false );
    f_nodes = reshape( f_nodes, num_lobes-1, num_lobes )';
    f_nodes = arrayfun( @(x) find(any(cell2mat(f_nodes(x,:)),2)), (1:num_lobes)', 'UniformOutput', false );
    f_faces = arrayfun( @(x) any(ismember(FaceArray{x},f_nodes{x}),2), (1:num_lobes)', 'UniformOutput', false ); 
    
    % Get distance map from fissure each node
    cppStructF = arrayfun( @(x) ClosestPointTriSurfV2(FaceArray{x}(f_faces{x},:),NodePos{x},NodePos{x}), (1:num_lobes)', 'UniformOutput', false );
    dist_n = cellfun( @(x) x.dist, cppStructF, 'UniformOutput', false );

    % Get sparse E X N adjacency matrix from ElementArray
    adj = cellfun(@(x) sparse(repmat((1:size(x,1))',1,4), x, ones(size(x))), ElementArray, 'UniformOutput', false);
    % Get average distance from fissure in a given element
    dist_e = arrayfun(@(x) adj{x}/4 * dist_n{x}, (1:num_lobes)', 'UniformOutput', false);
    
    % Divide distance into bins
    bin_edges = cellfun( @(x) 0:10:max(x), dist_e, 'UniformOutput', false); 
    bin_centers = cellfun( @(x) diff(x)/2 + x(1:end-1), bin_edges, 'UniformOutput', false);
    
    % Figure out which elements are in which bins
    in_bin = cellfun( @(x,y) x <= y(2:end) & x > y(1:end-1), dist_e, bin_edges, 'UniformOutput', false );

    % Compute fissure surface areas and volumes for each lobe
    fissure_area(i,:) = cellfun( @(x,y,z) TriSurfArea(x(z,:),y), FaceArray, NodePos, f_faces );
    lobe_volume(i,:) = cellfun( @(x,y) TetVolume(x,y), ElementArray, NodePos);

    % Plot ADI against distance for each model
    figure()
    hold on
    for j = 1:num_models
        % Load ADI data
        ADI_index = 2;
        ADI = cellfun( @(x) x(:,ADI_index), results{i,j}.ElementData(2:end), 'UniformOutput', false ); 
        
        % Calculate mean ADI for each bin
        meanADI = cellfun(@(x,y) x'*y ./ sum(y,1), ADI, in_bin, 'UniformOutput', false );  
        
        %cellfun(@(x,y,z) plot(x,y,'LineStyle',lineType{j},'Color',z,'LineWidth',lineWidth), bin_centers, meanADI, lineColor' )
        cellfun(@(x,y,z) stairs(x,[y,y(end)],'LineStyle',lineType{j},'Color',z,'LineWidth',lineWidth), bin_edges, meanADI, lineColor' )
    end    
    xlabel('Distance from fissure (mm)','FontSize',20)
    ylabel('ADI','FontSize',20)
    ax = gca;
    ax.FontSize = 13;
    legend(legend_text)
    box on
    title(subject)

    if debug_plots
        % View fissure surfaces
        figure()
        hold on
        cellfun( @(x,y,z)...
            trimesh( x(z,:), y(:,1), y(:,2), y(:,3) ), FaceArray, NodePos, f_faces )
        daspect([ 1 1 1])
        hold off
    
        % View distance contours
        figure()
        hold on
        cellfun( @(x,y,z)...
            patch( 'Faces', x, 'Vertices', y, 'FaceColor', 'interp', 'CData', z ), FaceArray, NodePos, dist_n )
        daspect([ 1 1 1])
        colorbar()
        hold off
    end
end