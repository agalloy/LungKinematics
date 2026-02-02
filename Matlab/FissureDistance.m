% WARNING! This code will not work for right lung!!!

%% Initialize Matlab
clear
clc

%% User Parameters
% Important directories
disp_dir = 'Y:\Documents\BioMOST_Research\Lung_FE\FEBio\FEBio_Runs\TLCtoFRC_PenaltyStep';
mesh_dir = 'Y:\Documents\BioMOST_Research\Lung_FE\FEBio\Meshes_v3';
mask_dir = 'Z:\Lung\lobe-by-lobe\halfsize-data\${SUBJECT}_backup';

% Important file patterns
mesh_pattern = '${SUBJECT}_LeftLung_Lobes_Mesh_v3.mat';
disp_pattern = '${SUBJECT}_LeftLung_Lobes_f0_ndata.txt';
mask_pattern = '${SUBJECT}_base${STATE}_lobemask_half.hdr';

% Subjects to analyze
subjects = "H5978"; % To analyze all subjects found in disp_dir use "all"
exclude = ""; % May be empty

% Initial and final state for FE model
t_int = [0,1];
% Corresponding states for mask files
states = ["TLC","FRC"];

% Regions to use in the FE Mesh
FE_regions = [2,3];
% Corresponding regions in the mask
mask_regions = [16,8];
% Names of regions for plots
region_names = ["Lower Lobe","Upper Lobe"];

plot_data = false;
flip_plot = true;

%% Create a subject list
if strcmp(subjects(1),"all")
    subject_pattern = replace( disp_pattern, "${SUBJECT}", "*" );
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
num_subjects = size(subject_list,1);

%% Main Loop
tic

num_regions = numel(FE_regions);
AllFissures = cell(num_subjects,1);
for i = 1:num_subjects
    subject = char(subject_list(i));
  
% Step 1: Load FE mesh data and extract lobe surfaces
    % Load mesh data
    mesh_file = fullfile( mesh_dir, mesh_pattern );
    mesh_file = replace( mesh_file, '${SUBJECT}', subject );
    load( mesh_file, "NodeArray", "ElementArray", "nID", "eID" );
    
    % Load disp data
    disp_file = fullfile( disp_dir, disp_pattern );
    disp_file = replace( disp_file, '${SUBJECT}', subject );
    t_get = t_int;
    try
        [ n_data, ~ ] = GetNodeData( disp_file, t_get );
    catch GND_exc
        fprintf('\nFailed to read disp_file. Subject: %s Model: %i \n',subject,j)
        continue
    end
        
    % Set reference  and deformed geometry to the desired time steps
    if t_get(1) > 0
        DispOffset = n_data{1}(:,2:4);
        X_Ref = NodeArray + DispOffset;
        % Offset final displacements accordingly
        DispArray = n_data{2}(:,2:4) - DispOffset;
    else
        X_Ref = NodeArray;
        DispArray = n_data{end}(:,2:4);
    end 
    X_Def = X_Ref + DispArray;
    
    % Get the surface nodes and faces for each lobe
    FE_Surfaces = struct('Faces', cell(1, num_regions), 'Nodes', cell(1, num_regions));
    for j = 1:num_regions  
        % Get the full lobar surface
        FE_Surfaces(j).Faces = FESurface( ElementArray(eID==FE_regions(j),:) );
        FE_Surfaces(j).Nodes = unique( FE_Surfaces(j).Faces );
    end
    
    % Isolate the fissures from each lobar surface
    FE_Fissures = struct( 'Faces', cell(1,num_regions), 'Nodes', cell(1,num_regions),...
                          'Dist1', cell(1,num_regions), 'Dist2', cell(1,num_regions) );
    for j = 1:num_regions
        % Determine the points on the current lobe in contact with the other lobes
        other_regions = 1:num_regions;
        other_regions(j) = [];
        cppStruct = ClosestPointTriSurfV2( FE_Surfaces(other_regions).Faces, X_Ref, X_Ref(FE_Surfaces(j).Nodes,:) );
        % Isolate the fissure nodes and faces
        f_nodes = find( cppStruct.dist < 1 | cppStruct.inside );
        FE_Fissures(j).Nodes = FE_Surfaces(j).Nodes( f_nodes );
        f_faces = find( all( ismember( FE_Surfaces(j).Faces, FE_Fissures(j).Nodes ), 2 ) );
        FE_Fissures(j).Faces = FE_Surfaces(j).Faces( f_faces, : );
    end
    
% Step 2: Load CT mask data and get distance maps for each lobe
    % Load CT masks
    mask_file = fullfile( mask_dir, mask_pattern );
    % Read first state
    mask_file1 = replace( mask_file, {'${SUBJECT}','${STATE}'}, {subject,char(states(1))} );
    mask_info1 = analyze75info( mask_file1 );
    mask1 = analyze75read( mask_info1 );
    mask1 = CleanLobeMask(mask1);
    voxel_size = mask_info1.PixelDimensions(1:3);
    % Read second state
    mask_file2 = replace( mask_file, {'${SUBJECT}','${STATE}'}, {subject,char(states(2))} );
    mask_info2 = analyze75info( mask_file2 );
    mask2 = analyze75read( mask_info2 );
    mask2 = CleanLobeMask(mask2);
    
    % Generate distance transforms of each region
    dist_map1 = cell(num_regions,1);
    dist_map2 = cell(num_regions,1);
    for j = 1:num_regions
        dist_map1{j} = logic2levelset( mask1==mask_regions(j), voxel_size );
        dist_map2{j} = logic2levelset( mask2==mask_regions(j), voxel_size );
    end

% Step 3: Evaluate the distance maps at the FE mesh surface node locations
    for j = 1:num_regions
        % Get a set of query points in the reference and deformed states
        Xq_Ref = X_Ref( FE_Fissures(j).Nodes, : );
        Xq_Def = X_Def( FE_Fissures(j).Nodes, : );
        
        % Convert spatial points to image indices
        Iq_Ref = Xq_Ref ./ voxel_size .* [1,-1,1] + [1,1,1];
        Iq_Def = Xq_Def ./ voxel_size .* [1,-1,1] + [1,1,1];
        
        % Sample the distance maps at the query points
        FE_Fissures(j).Dist1 = interp3( dist_map1{j}, Iq_Ref(:,1), Iq_Ref(:,2), Iq_Ref(:,3), 'linear' );
        FE_Fissures(j).Dist2 = interp3( dist_map2{j}, Iq_Def(:,1), Iq_Def(:,2), Iq_Def(:,3), 'linear' );
        
        if plot_data
            % Display the distance contours on the FE meshes
            v1 = nan(size(X_Ref,1),1);
            v1( FE_Fissures(j).Nodes ) = FE_Fissures(j).Dist1;
            v2 = nan(size(X_Def,1),1);
            v2( FE_Fissures(j).Nodes ) = FE_Fissures(j).Dist2;
            figure()
            subplot(1,2,1)
            patch( 'Faces', FE_Fissures(j).Faces, 'Vertices', X_Ref, 'FaceColor', 'interp', 'CData', abs(v1) )
            daspect([1 1 1])
            colorbar
            if flip_plot
                set(gca, 'Zdir', 'reverse')
                set(gca, 'Ydir', 'reverse')
            end
            title(sprintf('Segmentation distance map: Reference Configuration\n%s',region_names(j)))
            subplot(1,2,2)
            patch( 'Faces', FE_Fissures(j).Faces, 'Vertices', X_Def, 'FaceColor', 'interp', 'CData', abs(v2) )
            daspect([1 1 1])
            colorbar
            if flip_plot
                set(gca, 'Zdir', 'reverse')
                set(gca, 'Ydir', 'reverse')
            end
            title(sprintf('Segmentation distance map: Deformed Configuration\n%s',region_names(j)))

            % Display cumulative segmentation distance histograms
            bin_edges = 0 : 1 : ceil(max(max(abs(FE_Fissures(j).Dist1)),max(abs(FE_Fissures(j).Dist2))));
            figure()
            hold on
            histogram( abs(FE_Fissures(j).Dist1), bin_edges, 'Normalization', 'cdf' )
            histogram( abs(FE_Fissures(j).Dist2), bin_edges, 'Normalization', 'cdf' )
            title(sprintf('Segmentation distance normalized cumulative histogram\n%s',region_names(j)))
            legend( 'Reference Config', 'Deformed Config' )
            ylim([0,1])
        end
    end
    % Store results
    AllFissures{i} = FE_Fissures;
end

%% Analyze results
p90_ref = nan( num_subjects, num_regions );
p90_def = nan( num_subjects, num_regions );
for i = 1:num_regions
    p90_ref(:,i) = cell2mat( cellfun(@(x) prctile(abs(x(i).Dist1),90), AllFissures, 'UniformOutput', false) );
    p90_def(:,i) = cell2mat( cellfun(@(x) prctile(abs(x(i).Dist2),90), AllFissures, 'UniformOutput', false) );
end

% for i = 1:num_regions
%     p90_ref(:,i) = cell2mat( cellfun(@(x) mean(abs(x(i).Dist1)), AllFissures, 'UniformOutput', false) );
%     p90_def(:,i) = cell2mat( cellfun(@(x) mean(abs(x(i).Dist2)), AllFissures, 'UniformOutput', false) );
% end


% Display results in table
table_names = [ repmat( {'${REGION} - Ref'}, 1, num_regions ); repmat( {'${REGION} - Def'}, 1, num_regions ) ];
table_names = reshape(strrep( table_names, repmat({'${REGION}'},2,num_regions), repmat(region_names,2,1) ),1,[]);
out_table = table( 'Size', [num_subjects,1+2*num_regions], 'VariableTypes', ["string", repmat("double",1,2*num_regions)] );
out_table{:,1} = subject_list;
out_table.Properties.VariableNames = [ {'Subject'}, table_names ];
out_table(:,2*(1:num_regions)) = array2table(p90_ref);
out_table(:,2*(1:num_regions)+1) = array2table(p90_def);
disp('90th percentile segmentation distance:')
disp(out_table)

%%
toc