function [y] = problem3(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = length(x);
%y = zeros(n, 1);
m = (1:n)'./n;

y = m.*exp(x) -1;
return
% for i = 1:n
%     y(i) = (i/n)*exp(x(i)) - 1;
% end