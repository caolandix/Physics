clear; clc;
numGraphPoints = 1000;
k = 0.1;
radius = 0.3;
N = 100;
M1 = 1.9;
dy = radius / numGraphPoints;
T1 = 300;
T2 = 3000;
lambda = 0.04;
V2 = calc_speed_of_sound(T2);

for iter = 1:1:numGraphPoints
    yi = iter * dy;
    alpha = asin(yi / radius);
    xi = radius * (1 - cos(alpha));
    M1n = M1 * sin(alpha);
    M2n = 0.00137077 + 2.92163 * alpha - 4.62126 * alpha ^ 2 + 4.24972 * alpha ^ 3 - 1.86993 * alpha ^ 4 + 0.312301 ...
        * alpha ^ 5;
    gamma = atan((M2n / M1n) * (sqrt(T1 / T2)) * tan(alpha));
    V2V1 = sqrt((T2 / T1) * ((M2n / M1n) ^ 2) * exp(-(yi ^ 2) / (lambda ^ 2)) * (cos(alpha) ^ 2) + sin(alpha) ^ 2);
    Xi = (xi - k * radius) * (1 - V2V1 * cos(gamma));
    Yi = yi - V2V1 * (k * radius - xi) * sin(gamma);
    plot(Xi, Yi);
    hold on
    y = [Yi, Xi];
end

function [retval] = trans_left_eq(M1n, M2n, k)
    val1 = 1 / (M1n * (k - 1));
    val2 = 2 * k * M1n ^ 2 - (k - 1);
    val3 = (k - 1) * M1n ^ 2 + 2;
    val4 = 2 * k * M2n ^ 2 - (k - 1);
    val5 = 2 * k * M1n ^ 2 - (k - 1);
    power = k - (1 / 2 * k);
    retval = (val1) * sqrt(abs(val2 * val3) * (val4 / val5) ^ power) - 1;
end

function [retval] = calc_speed_of_sound(time)
    retval = 643.855 * sqrt(time / 273.15);
end