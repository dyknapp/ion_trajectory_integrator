res = 1000;
test_matrix = normrnd(0, 1, res);
test_coords = normrnd(res / 2.0, res / 20.0, [1, 2]);

for idx = 1:2
    test_coords(idx) = max(1, test_coords(idx));
    test_coords(idx) = min(res, test_coords(idx));
end

disp(linInterpolate2D(test_matrix, test_coords(1), test_coords(2), 1.0))