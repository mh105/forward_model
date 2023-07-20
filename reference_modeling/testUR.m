function [ flag ] = testUR(R, nonUR)
% Used to test whether a reference operator is a unipolar reference
% operator by testing the three properties of unipolar reference: 
% 1) Memorylessness
% 2) Rank deficient by 1
% 3) Orthogonal projector centering

% Non-unipolar references have the last two properties as well
% 1) Rank deficient by 1
% 2) Orthogonal projector centering

if nargin < 2
    nonUR =false;
end

assert(size(R,1) == size(R,2), 'R operator is not a symmetrical matrix.')

assert(all(abs(sum(R,2)-0) < 10^-12), 'R operator does not sum to 0 for some rows.')

flag = true;

% construct a common average to multiply
avg = eye(size(R,1)) - ones(size(R,1))*1/size(R,1);

% 1) Memorylessness
if ~nonUR
    if sum(ismembertol(avg*R, avg, 10^-6), 'all') ~= length(R(:))
        flag = false;
    end
else
    if sum(ismembertol(R*avg, R, 10^-6), 'all') ~= length(R(:))
        flag = false;
    end
end

% 2) Rank deficient by 1
if rank(R) ~= size(R,1) - 1
    flag = false;
end

% 3) Orthogonal projector centering
if sum(ismembertol(pinv(R)*R, avg, 10^-6), 'all') ~= length(R(:))
    flag = false;
end

end

