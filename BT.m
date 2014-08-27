function [ B ] = BT()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
B = @(beta,T)(1-exp(-beta*T))/beta;
end