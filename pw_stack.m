function [pwx] = pw_stack(x,order)
% PHASE-WEIGHTED STACKING 
% input:
% data
% phase coefficient penalty order
%
% output:
% phase-weighted stack
% phase-weight coefficient vector
S = zeros(size(x));
for ii = 1:numel(x(1,:))
S(:,ii) = hilbert(x(:,ii));
end
S = S./abs(real(S));

% get phase information
P = (angle(S));
% compute phase stack coefficients
C = 1/ii*abs(sum(exp(1i.*P),2));
% C = conv(C,hamming(150)/sum(hamming(150)),'same');
% phase-weighted cross-correlation
% stack_PW = (ch1(:,1).*C(1:numel(ch1(:,1))).^order+ch1(:,2).*C(1:numel(ch1(:,1))).^order+ch1(:,3).*C(1:numel(ch1(:,1))).^order)./3;
pwx = (x.*repmat(C,1,ii).^order);

end
