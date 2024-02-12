reps = 1000000;

for idx = 1:reps
    l = randi(100) -1;
    m = randi(l + 1) - 1;
    phi = 2 * pi * rand();
    theta = 2 * pi * rand();
    [r, c] = Ylm_f(int32(l), int32(m), theta, phi)
end

%%

% l = 322;
% m = 31;
% 
% sv = min(l-m,l+m) + 1;
% ev = max(l-m,l+m);
% fact_ratio = 1.;
% if (ev > sv)
%     for factor = sv:ev
%         fact_ratio = fact_ratio * double(factor);
%     end
% end
% if ((l-m)<(l+m))
%     fact_ratio = 1./fact_ratio;
% end
% fact_ratio

function Ylm = Ylm(n,m,theta,phi)
    theta = theta(:);
    phi = phi(:)';
    Pn = legendre(n,cos(theta),'norm');
    if m >= 0
      Pn = (-1)^m*Pn(abs(m)+1,:);
    else
      Pn = Pn(abs(m)+1,:);
    end
    Ylm = (1/sqrt(2*pi))*Pn'*exp(i*m*phi);
end