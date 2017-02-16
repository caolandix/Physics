clear; clc;

k = 0.1;
R = 0.3;
N = 100;
M1 = 1.9;
dy = R / N;
T1 = 100;
T2 = 3000;
lambda = 0.04;
for iter = 1:1:100
    yi = iter * dy;
    alpha = asin(yi / R);
    xi = R * (1 - cos(alpha));
    M1n = M1 * sin(alpha);
    M2n = 0.00137077 + 2.92163 * alpha - 4.62126 * alpha ^ 2 + 4.24972 * ...
        alpha ^ 3 - 1.86993 * alpha ^ 4 + 0.312301 * alpha ^ 5;
    gamma = atan((M2n / M1n) * (sqrt(T1 / T2)) * tan(alpha));
    V2V1 = sqrt((T2 / T1) * ((M2n / M1n) ^ 2) * exp(-(yi ^ 2) / ...
        (lambda ^ 2)) * (cos(alpha) ^ 2) + sin(alpha) ^ 2);
    Xi = (xi - k * R) * (1 - V2V1 * cos(gamma));
    Yi = yi - V2V1 * (k * R - xi) * sin(gamma);
    plot(Xi, Yi);
    hold on
    y = [Yi, Xi];
end