% Calculate the (l, m) spherical harmonic at (theta, phi).
%
% 12.11.2023,  dknapp: wrote the example

function Ylm = Ylm(l, m, theta, phi)
    arguments
        l int32
        m int32
        theta (1, :) double
        phi (1, :) double
    end
    theta = theta(:);
    phi = phi(:)';
    Pn = legendre(l, cos(theta), 'norm');
    % sign change required for odd positive m
    if m >= 0
      Pn = double((-1)^m) * Pn(abs(m) + 1,:);
    else
      Pn = Pn(abs(m) + 1, :);
    end
    Ylm = (1/sqrt(2*pi)) * Pn' * exp(1j * double(m) * phi);
end