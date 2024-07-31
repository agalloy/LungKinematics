% Compile a folder full of lung FEA results into TecPlot viewable data

%% Initialize Matlab
clear
clc
addpath('Y:\Documents\BioMOST_Research\TecPlotTools')
tic

%% User Parameters
% FEBio results input directory
feb_dir = 'Y:\Documents\BioMOST_Research\Lung_FE\FEBio\FEBio_Runs\TLCtoFRC_PenaltyStep';
% FEBio mesh directory
mesh_dir = 'Y:\Documents\BioMOST_Research\Lung_FE\FEBio\Meshes_v3';
% TecPlot data output direcctory
tec_dir = 'Y:\Documents\BioMOST_Research\Lung_FE\TecPlot\FEBioRuns\TLCtoFRC_PenaltyStep';

% Important string patterns for finding files
file_pat = 'H*_RightLung_Lobes_f*.feb'; % Wildcard notation
subject_pat = 'H\d\d\d\d'; % RegExp notation
mesh_pat = '$SUBJECT_RightLung_Lobes_Mesh_v3.feb'; % $SUBJECT will be replaced

%% Loop through all .feb files in the FEBio results directory
feb_files = dir(fullfile(feb_dir,file_pat));
num_files = length(feb_files);
fprintf('\nProcessing %i files.\n',num_files);

for i = 1:num_files
    % Extract the title part of the file
    [~,title,~] = fileparts(feb_files(i).name);
    
    % Get the subject from the file name
    subject = regexp(title,subject_pat,'match');
    subject = subject{1};
    
    % Check if the appropriate mesh file for the subject exists
    mesh_name = strrep(mesh_pat,'$SUBJECT',subject);
    mesh_file = fullfile(mesh_dir,mesh_name);
    if isfile(mesh_file)
        % Load .feb file string
        fID = fopen(mesh_file);
        mesh_txt = fscanf(fID,'%c');
        fclose(fID);

        % Loop through each element domain and fill up Element connectivity array
        elem_txt = extractBetween(mesh_txt,'<Elements', ...
                                 sprintf('\r\n\t\t</Elements>') );
        ElementArray = cell( size(elem_txt,1), 1 );
        ElementID = cell( size(elem_txt,1), 1 );
        z_title = cell( size(elem_txt,1), 1 );
        for j = 1:length(elem_txt)
            % Get element type
            e_type = regexp(elem_txt{j},'type="(\w+)"','tokens');
            e_type = e_type{1}{1};
            
            % Get zone title
            z_title{j} = regexp(elem_txt{j},'name="(\w+)"','tokens');
            z_title{j} = z_title{j}{1}{1};
            
            % Correlate local element indices to global indices
            ElementID{j} = regexp(elem_txt{j},'elem id="(\d+)"','tokens');
            ElementID{j} = cellfun(@(X) str2double(X{1}), ElementID{j});
            
            % Get local connectivity data
            if strcmp(e_type,'tri3')
                ElementArray{j} = regexp(elem_txt{j},'>(\d+),(\d+),(\d+)<','tokens');
            elseif strcmp(e_type,'tet4')
                ElementArray{j} = regexp(elem_txt{j},'>(\d+),(\d+),(\d+),(\d+)<','tokens');
            end
            ElementArray{j} = cat(1,ElementArray{j}{:});
            ElementArray{j} = cellfun(@str2double,ElementArray{j});
        end
        clear elem_txt

        % Loop through each node domain and fill up Node position array
        node_txt = extractBetween(mesh_txt,'<Nodes', ...
                                 sprintf('\r\n\t\t</Nodes>') );
        NodeArray = cell( size(node_txt,1), 1 );
        NodeID = cell( size(node_txt,1), 1 );
        for j = 1:length(node_txt)
            % Create an index from local node data to global node data
            NodeID{j} = regexp(node_txt{j},'node id="(\d+)"','tokens');
            NodeID{j} = cellfun(@(X) str2double(X{1}), NodeID{j});
            
            % Get node position data
            NodeArray{j} = regexp(node_txt{j},'>([-\d\.]+),([-\d\.]+),([-\d\.]+)<','tokens');
            NodeArray{j} = cat(1,NodeArray{j}{:});
            NodeArray{j} = cellfun(@str2double,NodeArray{j});
        end
        clear node_txt
    else
        warning('Could not find a mesh file for current subject %s. Data will not be processed.',subject);
    end
    
    % Read nodal data file
    ndata_file = fullfile(feb_dir,[title,'_ndata.txt']);
    ndata = readmatrix(ndata_file);
    u = ndata(:,2:4);
    clear ndata
    NodeData = cell(size(NodeArray,1),1);
    for j = 1:length(NodeData)
        NodeData{j} = u(NodeID{j},:);
    end
    
    % Read element data file
    edata_file = fullfile(feb_dir,[title,'_edata.txt']);
    edata = readmatrix(edata_file);
    E = edata(:,2:4);
    clear edata
    
    % Post-process element data
    lambda = E(:,[3,2,1]);
    lambda = 1./sqrt(2*lambda + 1);
    J = lambda(:,1).*lambda(:,2).*lambda(:,3);
    ADI = sqrt( (lambda(:,1)./lambda(:,2)-1).^2 + (lambda(:,2)./lambda(:,3)-1).^2 );
    SRI = atan( lambda(:,3).*(lambda(:,1) - lambda(:,2)) ./ (lambda(:,2).*(lambda(:,2) - lambda(:,3))) ) *2/pi;
    
    % Assemble element data
    ElementData = cell(size(NodeArray,1),1);
    for j = 1:length(ElementData)
        ElementData{j} = [ J(ElementID{j},:), ADI(ElementID{j}), SRI(ElementID{j}) ];
    end
    
    % Renumber element connectivity on a zone-by-zone basis
    ElementArrayL = ElementArray;
    for j = 1:size(ElementArrayL,1)
        % Create an inverse index from global to local numbering
        NodeID_inv = nan(max(NodeID{j}),1);
        NodeID_inv(NodeID{j}) = (1 : length(NodeID{j}) )';
        % Renumber array to local node numbers
        ElementArrayL{j} = NodeID_inv(ElementArray{j}); 
    end
    
    % Write to TecPlot
    tecStruct.title = title;
    tecStruct.dir = tec_dir;
    tecStruct.VarNames = ["U1","U2","U3","J","ADI","SRI"];
    tecStruct.ZoneNames = z_title;
    tecStruct.NodeArray = NodeArray;
    tecStruct.ElementArray = ElementArrayL;
    tecStruct.NodeData = NodeData;
    tecStruct.ElementData = ElementData;
    WriteTecPlotData(tecStruct);
    
    % Display files remaining update
    fprintf('\n%i files remaining.\n',num_files - i);
end

%%
toc