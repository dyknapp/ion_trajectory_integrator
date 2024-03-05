function handle = depthPlotLines( x,y,z, fillColor, mSizes )
%DEPTHPLOT Plots a series of points but lerps the marker color inbetween
%the given values to give an illusion of depth

if any(size(fillColor) ~= [1 3])
    error('colors must be specified as 1x3');
end

if nargin < 5
    mSizes = [100 250];
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
        
        miD = min(dist(:));
        maD = max(dist(:));
        
        %normalise distances between 0 and 1.
        dist = (dist - miD)/(maD-miD);
    end

    function c = getColor(dist)
        dist3 = repmat(dist, 1, 3);
        a = repmat(fillColor .* 0.5, size(dist,1), 1);
        b = repmat(fillColor, size(dist,1), 1);

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




%handle = scatter3(x,y,z,markerSize,fillC,'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
handles = cell(size(x,1),1);
hold on;
for i=1:size(x,1)
    u = x(i,:)'; v = y(i,:)'; w = z(i,:)';
    dist = getNormalisedDepths(u,v,w,getCameraViewDirection(gca));

    markerSize = getSize(dist);
    handles{i} = patch([u' NaN], [v' NaN], [w' NaN], [dist' NaN], 'edgecolor', 'interp', 'linewidth', 2);
end
hold off;
colormap(getColor([0:0.01:1]'));


%Start a timer that will update the view on rotation
% t = timer('TimerFcn', @(~,~,~)callback(), 'StartDelay', 0.2, 'Period', 0.1,...
% 'ExecutionMode', 'fixedDelay');
% start(t);
end

