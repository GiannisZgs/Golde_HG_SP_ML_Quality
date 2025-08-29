%% Dynamic Time Warping calculation
function dist = dtw(x, y)
    [~, dist] = dtw_c(x, y);
end

function [D, dist] = dtw_c(x, y)
    nx = length(x);
    ny = length(y);
    D = inf(nx+1, ny+1);
    D(1,1) = 0;
    x = (x - mean(x)) / std(x);
    y = (y - mean(y)) / std(y);
    for i = 2:nx+1
        for j = 2:ny+1
            cost = (x(i-1) - y(j-1))^2;
            D(i,j) = cost + min([D(i-1,j), D(i,j-1), D(i-1,j-1)]);
        end
    end
    dist = sqrt(D(nx+1, ny+1));
end
