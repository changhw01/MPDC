function r = get_expansion(A,b_in, flag_silent)

% A: m-by-n matrix
% b: m-by-1 vector
%
% r: 1-by-n vector (answer)

if nargin < 3
    flag_silent = false;
end

n = size(A,2);
r = false(1,n);

b = b_in;
idx = find(b,1,'first');
n_reduction = 0;

while ~isempty(idx) && (n_reduction<n)
    idx_c = find(A(idx,:),1,'first');
    r(idx_c) = true;
    b = b & (~A(:,idx_c));        
    n_reduction = n_reduction + 1;    
    idx = find(b,1,'first');
end

b_check = any(A(:,r),2);
if any(xor(b_in,b_check))
    if ~flag_silent
        disp('Fail to compute bases expansion!')
    end
    r = NaN;
end
