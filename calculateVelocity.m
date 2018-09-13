function [sessionData] = calculateVelocity(sessionData)
%CALCULATEVELOCITY adds columns with distance traveled and velocity of movement

trans_x_diff = sessionData.trans_x(2:end) - sessionData.trans_x(1:end-1);
trans_y_diff = sessionData.trans_y(2:end) - sessionData.trans_y(1:end-1);
time_diff = sessionData.timestamp(2:end) - sessionData.timestamp(1:end-1);
T = [ trans_x_diff, trans_y_diff ];
dist = sqrt(sum(T .* T, 2));
sessionData.dist = [0.0; dist];
sessionData.velocity = [0.0; dist ./ double(time_diff) * 1000];
end

