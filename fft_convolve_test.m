%% signal

signal_t = linspace(0, 2*pi, 1000);
signal = sin(signal_t) + sin(signal_t * 4);

plot(signal_t, signal);

%% filter

filter_t = linspace(0, pi, 1000);
filter = sin(filter_t).^16;

plot(filter_t, filter);

%% ffts

signal_fft = fft(signal);
filter_fft = fft(filter);

%% convolve
result = convnfft(signal, filter);
plot(result);

function A = convnfft(A, B)
    nd = 2;
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