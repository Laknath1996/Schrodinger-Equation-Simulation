function [zc,idx] = zero_cross_count(x)
    % find the number of zero crossings
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
    idx = zci(x);
    zc = length(idx)-2;
    
    
%     N = length(x);
%     % No. of crossings
%     zc = 0;
%     for c = 2 : N
%         if x(c-1)*x(c) < 0
%             zc = zc + 1;
%         end
%     end
end

