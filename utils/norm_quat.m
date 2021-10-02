% 
% LibQPEP: A Library for Globally Optimal Solving Quadratic Pose Estimation Problems (QPEPs),
%          It also gives highly accurate uncertainty description of the solutions.
%
%
% Article: 
%      Wu, J., Zheng, Y., Gao, Z., Jiang, Y., Hu, X., Zhu, Y., Jiao, J., Liu, M. (2020)
%           Quadratic Pose Estimation Problems: Unified Solutions, 
%           Solvability/Observability Analysis and Uncertainty Description 
%           in A Globally Optimal Framework.
%
%
% Authors:      Jin Wu and Ming Liu
% Affiliation:  Hong Kong University of Science and Technology (HKUST)
% Emails:       jin_wu_uestc@hotmail.com; eelium@ust.hk
% Websites:     https://zarathustr.github.io
%               https://ram-lab.com


function q_ = norm_quat(q)
N = 1 ./ sqrt(q(:, 1) .* q(:, 1) + q(:, 2) .* q(:, 2) + q(:, 3) .* q(:, 3) + q(:, 4) .* q(:, 4));
q_ = [q(:, 1) .* N, q(:, 2) .* N, q(:, 3) .* N, q(:, 4) .* N];
end