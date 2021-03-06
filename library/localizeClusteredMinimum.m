function [x,y] = localizeClusteredMinimum(powerVec,signal,k,M)
% Localize a VLP positioning signal using fingerprint database
% Uses minimum difference (not squared) between signal and powerVec

aSearch = min(powerVec - signal,[],2);
[B,ind] = mink(abs(aSearch),k); % Grab k smallest (10)

% Convert back to subscript index
[v,h] = ind2sub(M,ind);

% Average the position and round to integer
y = round(mean(v));
x = round(mean(h));

% Score error with Euclidean distance from real position
% error = sqrt((real.x(i) - guess.x(i)).^2 + (real.y(i) - guess.y(i)).^2);
        
end