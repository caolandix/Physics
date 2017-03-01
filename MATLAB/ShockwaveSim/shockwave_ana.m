clear; clc;

M1 = 1.9;
R = 0.3;
T1 = 300.0;
T2 = 3000.0;
N = 50;
time = 0.5;
maxCalcs = 50;
calc_step = 1;
arr_xi = zeros(maxCalcs);
arr_yi = zeros(maxCalcs);
arr_Xi2 = zeros(maxCalcs);
arr_Yi2 = zeros(maxCalcs);
arr_Xi3 = zeros(maxCalcs);
arr_Yi3 = zeros(maxCalcs);

for i = 1:maxCalcs
    yi = i * R / N;
    arr_yi(i) = yi;

    xi = R * sqrt((R ^ 2) - (yi ^ 2));
    arr_xi(i) = xi

    alpha = acos(1 - (xi / R));
    M1n = cos(alpha) * (1 - (xi / R));
    x = M1n;
    M2n = 0.0013707729340115783 + 2.921534777588217 * M1n - ...
        4.621262346264438 * M1n ^ 2 + 4.249722375287607 * x ^ 3 - ...
        1.8699297117154714 * M1n ^ 4 + 0.3123011824283489 * M1n ^ 5;
    ratio1 = ((yi / R) ^ 2) + 10 * ((1 - xi / R) ^ 2) * ((M2n / M1n) ^ 2);
    ratio2 = sqrt(ratio1);
    xii = R - xi;
    gamma = alpha - atan(0.3162 * (M2n / M1n)) * tan(alpha);
    arr_Xi2(i) = (ratio2 * cos(gamma) - 1) * (time * R - xi);
    arr_Yi2(i) = yi - ratio2 * (time * R - xi) * sin(gamma);

    ratio13 = ((yi / R) ^ 2) + 10 * ((1 - xi / R) ^ 2);
    ratio23 = sqrt(ratio13);
    gamma1 = alpha - atan(0.3162 * tan(alpha));

    arr_Xi3(i) = (ratio23 * cos(gamma1) - 1) * (time * R - xi);
    arr_Yi3(i) = yi - ratio23 * (time * R - xi) * sin(gamma1);
end

plot(arr_xi, arr_yi, 'r');
hold on;
plot(arr_Xi2, arr_Yi2);
hold on;
plot(arr_Xi3, arr_Yi3);
hold on;
