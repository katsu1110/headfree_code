function e = smooth_eye(e, samplingRate)
% smooth eye positions using moving average

% transform eye-positions into velocities
len_eye = length(e);
vel = zeros(1,len_eye);
for i = 3:len_eye-2
      vel(i) = samplingRate*(e(i+2)+e(i+1)-e(i-1)-e(i-2))/6;
end

% reconstruct eye position
e = e(1) + cumsum(vel, 2)/samplingRate;