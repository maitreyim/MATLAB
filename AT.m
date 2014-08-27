function [ A ] = AT()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
A = @(alpha,beta,sigma,T) exp(((sigma^2)/(2*beta^2) - alpha/beta)*T + ...
            (alpha/(beta^2) - (sigma^2)/(beta^3))*(1 - exp(-1*beta*T)) + ...
            ((sigma^2)/(4*(beta^3)))*(1 - exp(-2*beta*T)));
end
