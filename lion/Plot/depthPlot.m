function handle = depthPlot( x,y,z, fillColor, mSizes, markertype, lightFactors, varargin)
%DEPTHPLOT Plots a series of points but lerps the marker color and size
%between given values to give an illusion of depth

ip = inputParser;
ip.addParameter('DistanceLimits', []);
ip.KeepUnmatched = 1;
ip.parse(varargin{:});

if (size(fillColor, 2) ~= 3 && ~(size(fillColor, 1) ~= 1 || size(fillColor, 1) ~= size(x,1)))
    error('colors must be specified as Nx3');
end

if size(x,2) ~= 1 || size(y,2) ~= 1 || size(z,2) ~= 1
    error('only lists of points supported');
end

if nargin < 5
    mSizes = [100 250];
end

if nargin < 6
    markertype = 'filled';
end

if nargin < 7
    lightFactors = [ 0.5 1 ];
end

if size(fillColor) == [1 3]
   % repeat fill color for all atoms in list
   fillColor = repmat(fillColor, size(x,1), 1);
end

% Create plot varargs
plotVarargs = {};
unmatched = ip.Unmatched;
for argName=fields(unmatched)'
    plotVarargs{end+1} = argName{1};
    plotVarargs{end+1} = unmatched.(argName{1});
end


    function direction = getCameraViewDirection(axForGraph)
        %gets camera view direction for given axes
        [az,el] = view(axForGraph);
        direction = rotz(az)*rotx(-el)*[0 1 0]';
    end

    function dist = getNormalisedDepths(x,y,z,dir)
        pos = [x y z];
        cam = repmat(dir', size(x,1), 1);
        dist = sum(pos.*cam, 2);
        
        if isempty(ip.Results.DistanceLimits)
            miD = min(dist(:));
            maD = max(dist(:));
        else
            miD = ip.Results.DistanceLimits(1);
            maD = ip.Results.DistanceLimits(2);
        end
        
        %normalise distances between 0 and 1.
        dist = (dist - miD)/(maD-miD);
        dist = min(dist, 1);
        dist = max(dist, 0);
    end

    function c = getColor(dist)
        dist3 = repmat(dist, 1, 3);
        b = fillColor * lightFactors(2);
        a = fillColor .* lightFactors(1);

        c = a .* dist3 + (1 - dist3) .* b;
    end

    function s = getSize(dist)
        s = mSizes(1) .* dist + mSizes(2) .* (1-dist);
    end

    function callback()%obj,event_obj)
        %refresh distances
        %a = event_obj.Axes;
        try
        dist = getNormalisedDepths(x,y,z,getCameraViewDirection(handle.Parent));
        fillC = getColor(dist);
        markerSize = getSize(dist);
        
        %update handle with new data
        set(handle, 'CData', fillC, 'SizeData', markerSize);
        catch e
            stop(t);
            delete(t);
        end
    end

%By default camera has az=-45, el = 45. This looks along [1 1 -1]
view(-45,45);
dist = getNormalisedDepths(x,y,z,getCameraViewDirection(gca));
fillC = getColor(dist);
markerSize = getSize(dist);
handle = scatter3(x,y,z,markerSize,fillC,markertype, 'LineWidth', 1.5, plotVarargs{:});
if strcmp(markertype, 'filled')
    set(handle,'MarkerEdgeColor', 'k');
end

%Start a timer that will update the view on rotation*
t = timer('TimerFcn', @(~,~,~)callback(), 'StartDelay', 0.2, 'Period', 0.1,...
'ExecutionMode', 'fixedDelay');
start(t);

% * Note: I experimented with a callback, but it's not as reliable and also
% this way allows you to use hold and and generate multiple depth plots.

end

