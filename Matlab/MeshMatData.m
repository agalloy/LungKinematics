% Create a .mat file with node and element for a mesh in .feb format

%% Initialize Matlab
clear
clc

%% User Parameters
% Input directory with mesh files
mesh_dir = 'H:\Documents\GitHub_Projects\MeshPipeline\FEBio\Meshes\TetFactorStudy';
% Directory to place .mat files
out_dir = mesh_dir;

% Input .feb file pattern
mesh_pattern = '${SUBJECT}_${SIDE}Lung_Lobes_tf${tf}_Mesh.feb';
% Output .mat file pattern
out_pattern = replace(mesh_pattern,'.feb','.mat');

% Subjects to analyze
subjects = "all"; % To analyze all subjects found in mesh_dir use "all"
exclude = ""; % May be empty

% Model parameters (1 x P string array)
model_params = ["${SIDE}","${tf}"];
% Values to set those parameters to (M x P string array)
model_values = ["Left","1.0"
                "Left","1.1"
                "Left","1.2"
                "Left","1.5"
                "Left","2.0"
                "Left","3.0"
                "Left","4.0"];
            
%% Create a subject list
if strcmp(subjects(1),"all")
    subject_pattern = replace( mesh_pattern, ["${SUBJECT}",model_params], repmat("*",1,1+numel(model_params)) );
    mesh_dir_files = dir(fullfile(mesh_dir,subject_pattern));
    subject_list = strings(length(mesh_dir_files),1);
    % Get all the subject names in folder
    for i = 1:length(mesh_dir_files)
        subject = mesh_dir_files(i).name;
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

%% Loop through each mesh and save its data in a .mat file
tic
for i = 1:num_subjects
    for j = 1:num_models
        fprintf( '\nWorking on mesh %d of %d...', [(i-1)*num_models+j, num_subjects*num_models] )
        subject = char(subject_list(i));

        % Get mesh data from .feb file
        mesh_name = replace( mesh_pattern, ["${SUBJECT}",model_params], [subject,model_values(j,:)] );
        mesh_file = fullfile( mesh_dir, mesh_name );
        [ NodeArray, ElementArray, nID, eID ] = GetFEBMesh( mesh_file );

        % Save variables into a .mat file
        out_name = replace( out_pattern, ["${SUBJECT}",model_params], [subject,model_values(j,:)] );
        out_file = fullfile( out_dir, out_name );
        save( out_file, "NodeArray", "ElementArray", "nID", "eID" )
    end
end
fprintf('\nFinished!\n')
toc