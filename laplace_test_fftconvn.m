%% Set up the problem
res = 1024;
dimensions = [res res res];
potential_array = zeros(dimensions);
boundary_conditions = zeros(dimensions);

d = 1.0;

% Boundary conditions
for x = 1:dimensions(1)
    for y = 1:dimensions(2)
        if   (x - double(dimensions(1) / 2.0))^2 ...
           + (y - double(dimensions(2) / 2.0))^2 ...
           < (double(dimensions(1)) / 5.0)^2
            boundary_conditions(x, y, round(0.4 * dimensions(3))) =  1.0;
            boundary_conditions(x, y, round(0.6 * dimensions(3))) =  1.0;
        end
    end
end

% Set up logical indexing
where_boundaries = boundary_conditions ~= 0.0;

% Apply boundary conditions
potential_array(where_boundaries) = boundary_conditions(where_boundaries);

imagesc(squeeze(potential_array(round(dimensions(1)/2), :, :)))
set(gca,'YDir','normal')

%%

filter = zeros([3 3 3]);
for x = -1:1
    for y = -1:1
        for z = -1:1
            if abs(x) == 1 && abs(y) == 0 && abs(z) == 0
                filter(x + 2, y + 2, z + 2) = 1.0;
            elseif abs(x) == 0 && abs(y) == 1 && abs(z) == 0
                filter(x + 2, y + 2, z + 2) = 1.0;
            elseif abs(x) == 0 && abs(y) == 0 && abs(z) == 1
                filter(x + 2, y + 2, z + 2) = 1.0;
            end
        end
    end
end

conv = convnfft(potential_array, filter);

%%
f = @() convnfft(potential_array, filter);
timeit(f)

%% FDM

i = 0;
while i < res * 4
    i = i + 1;
end

%%
% function res = naive_convolution(A, B)
% 
% end

%%
function A = convnfft(A, B)
    nd = 3;
    dims = 1:nd;
    dims = reshape(dims, 1, []); % row (needed for for-loop index)
    ifun = @(m,n) 1:m+n-1;
    lfftfun = @(l) l;
    % Do the FFT
    subs(1:ndims(A)) = {':'};
    for dim=dims
        m = size(A,dim);
        n = size(B,dim);
        % compute the FFT length
        l = lfftfun(m+n-1);
        A = fft(A,l,dim);
        B = fft(B,l,dim);
        subs{dim} = ifun(m,n);
    end
    % Pointwise product == convolution
    A = A.*B;
    % Back to the non-Fourier space
    for dim=dims
        A = ifft(A,[],dim);
    end
    A = real(A);
end % convnfft