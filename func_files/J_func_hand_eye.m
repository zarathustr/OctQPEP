function J = J_func_hand_eye(X, A, B)
J = 0;
if(length(size(A)) == 3)
    len = size(A, 3);
    for i = 1 : len
        res = A(:, :, i) * X - X * B(:, :, i);
        J = J + 1 / len * trace(res.' * res);
    end
else
    res = A * X - X * B;
    J = trace(res.' * res);
end
end