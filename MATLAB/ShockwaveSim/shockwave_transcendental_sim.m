
    
%     Unused code
%     numGraphPoints = 1000;
%     radius = 0.3;
%     M1 = 1.9;
%     dy = radius / numGraphPoints;
%     lambda = 0.04;
%     V2 = calc_speed_of_sound(T2);
    
    
    % Calculate and draw the left transcendental function
    clear; clc;
    k = 0.1;
    T1 = 300;
    T2 = 3000;
    trans_eq_r = [1:50];
    
    array_M2n = [1:50];
    % trans_eq_l = [1, 20];
    k = 0;
    idx = 1;
    M1n = 1.6;
    x = 1;
    for M2n=1:0.02:10
        trans_eq_r(x) = trans_right_eq(M1n, M2n, T1, T2);
    end
    %{
    M2n = 1.0;

    % Getting the transcendental left values with a set M2n of 1.0 over a
    % range of M1n: 1.0 - 2.8
    for M1n = 1.2:0.01:2.8
        trans_eq_l(idx) = trans_left_eq(M1n, M2n, k);
        x(idx) = M1n;
        idx = idx + 1;
    end
    plot(x, trans_eq_l);
    hold on
    
    % Calculate and draw the right transcendental function
    M1n = 1.0;
    idx = 1;
    x = [1, 20];
    trans_eq_r = [1, length(x)];
    
    % Getting the transcendental right values with a set M1n of 1.5 over a
    % range of M2n: 1.0 - 3.0    
    for M2n = 1.0:0.01:3.0
        trans_eq_r(idx) = trans_right_eq(M1n, M2n, T1, T2);
        x(idx) = M2n;
        idx = idx + 1;
    end
    %}
    plot(array_M2n, trans_eq_r);
    title('Plasma Shockwaves - Left and Right Transcendental Curves');
    hold on    


function [retval] = trans_right_eq(M1n, M2n, T1, T2)
    retval = M1n * (1 - (1 / M1n ^ 2)) - M2n * (1 - (1 / M2n ^ 2)) * sqrt(T2 / T1);
end

function [retval] = trans_left_eq(M1n, M2n, k)
    power = k - (1 / 2 * k);
    retval = (1 / (M1n * (k - 1))) * sqrt(abs((2 * k * M1n ^ 2 - (k - 1)) * ((k - 1) * M1n ^ 2 + 2)) * ...
        ((2 * k * M2n ^ 2 - (k - 1)) / (2 * k * M1n ^ 2 - (k - 1))) ^ power) - 1;
end

function [retval] = calc_speed_of_sound(time)
    retval = 643.855 * sqrt(time / 273.15);
end