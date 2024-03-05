function [timesteps, id, varargout] = readDump( filename, sortResults )
%READDUMP Reads data from the given dump file. The dump should be a file
%with atom quantities in the order 'id vargout' eg 'id vx vy vz', then read
%using [timesteps, id, vx, vy, vz] =readDump('dump.txt')

if (nargin < 2)
    sortResults = 1;
end

fHandle = fopen(filename, 'r');
cleaner = onCleanup(@() fclose(fHandle));

cf = 1; %current frame
ri = [ 0 0 0 ]; %timestep atoms
timesteps = zeros(0);

%Create atom format string for output.
of = '%d';
for i=1:nargout-2 %first two args are atomID and timestep respectively
    of = [of ' %g'];
end

    function receiveRequired(index)
        if ri(index)
            %If we've already seen this tag we must be on the next frame.
            %removed use of min for speed
            if ~ri(1) || ~ri(2) || ~ri(3)
                error('cannot advance frame: did not find all required data for current.');
            end
            
            ri = [ 0 0 0 ];
            cf = cf+1;
        end
        %whether we advance the frame or not, we must remember we have
        %seen this required item for the current frame.
        ri(index) = 1;
    end

str = fgets(fHandle);

%enumerate frames
while ~feof(fHandle)
    if length(str) >= 9
        %handle different tags
        switch str(7:9)
            
            %Timestep field
            %case 'ITEM: TIM'
            case 'TIM'
                receiveRequired(1);
                timestep = fgets(fHandle);
                [timesteps(cf),~] = sscanf(timestep, '%d');
                %case 'ITEM: ATO'
            case 'ATO'
                %check nargout is compatible with the number of spaces.
                % id_x_y_z -> number of spaces should equal nargout-2
                str = strtrim(str(13:end));
                ns = sum(str==' ');
                if ns ~= nargout-2
                    error('Dump file incompatible with specified output format.');
                end
                
                %read each atom's data.
                A = fscanf(fHandle, of, [nargout-1, natoms]);
                
                %initialise atomVars
                if ~exist('atomVars', 'var')
                    atomVars = zeros(size(A, 1), size(A, 2), 0);
                end
                
                try
                atomVars(:,:, end+1) = A;
                catch e
                   error('Error occured during readDump. Possible NaN entry for atomic properties, which may be evidence of an unphysical system.', e.message); 
                end
                
                receiveRequired(2);
                %case 'ITEM: NUM'
            case 'NUM'
                [natoms,~] = sscanf(fgets(fHandle),'%d');
                receiveRequired(3);
                
        end
    end
    
    str=fgets(fHandle);
end

if cf == 1 %not a single frame read as frame counter never advanced.
    error('Could not identify data');
end

%Prepare output for nargout
id = reshape(atomVars(1, :, :), size(atomVars, 2), size(atomVars, 3));


%reshape data to order indices correctly according to the read ids

if (sortResults)
    % sort according to the ordering of ids
    [~,j] = sort(id,1);
    for i=1:size(id, 3)
        atomVars(:,:,i) = atomVars(:, j(:,i), i);
        id(:,i) = id(j(:,i),i);
    end
end

try
    for i=1:nargout-2
        varargout{i} = reshape(atomVars(1+i, :, :), size(atomVars, 2), size(atomVars, 3));
    end
catch e
    warning('Error occured while loading simulation output. Check that ion trap parameters are stable!');
    throw e;
end

end
