%%% extractDistribution.m
%%% 31/07/2019
%%% Written by Soma Olasz
%%% 
%%% This script is created to extract the distribution function
%%% from the last time step of a NORSE calculation to Python
%%% languge.
%%% 
%%% 
%%% 
%%% 
%%% 
%%% 
%%% 

function double = extractDistribution(object)

    % take the last column of the NORSE distribuion field o.f
    double = object.f(:,end);
    

end