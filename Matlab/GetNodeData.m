% Read nodal data in ascii text format from FEBio
% INPUTS:
%   data_file = Filepath to the node data file
%   t_get = Array of time points to get data for.
%       If equal to -1 gets only last step
%       If equal to -2 gets all time steps
% OUTPUTS:
%   n_data = A cell array of all the node data indexed by retrieved step number
%   t_get = Index between retrieved step number and step time (output in
%       case different from input)

function [n_data,t_get] = GetNodeData( data_file, t_get )
    % Read node data file as text
    fID = fopen(data_file);
    data_txt = fscanf(fID,'%c');
    fclose(fID);
    % Remove carriage returns from text, if present
    data_txt = erase(data_txt,sprintf('\r'));
    
    % Find the time represented by each time step in the file
    time = regexp(data_txt,'\*Time  = ([\d\.]+)','tokens');
    time = cellfun( @str2double, time )';
    if t_get == -1
        t_get = time(end);
    elseif t_get == -2
        t_get = time;
    end
    
    % Get the nodal data title
    d_title = regexp(data_txt,'\*Data  = ([\w ]+)\n','tokens');
    d_title = d_title{1}{1};
    
    % Split text into each time step saved
    step_txt = split( data_txt, d_title );
    step_txt = step_txt(2:end);
    
    % Collect data for all time steps in t_get
    n_data = cell(size(t_get));
    for i = 1:numel(t_get)
        if t_get(i) == 0
            continue
        end
        
        step_index = find(ismember(time,t_get(i)));
        if isempty(step_index)
            error('Invalid time step.')
        end
        % Remove non-data lines at the end
        step_data = step_txt{step_index};
        if step_index < length(step_txt)
            step_data = extractBefore(step_data,'*');            
        end
        % Read in nodal data
        n_data{i} = sscanf(step_data,'%d %f %f %f\n');
        n_data{i} = reshape(n_data{i},4,[])';
    end
end