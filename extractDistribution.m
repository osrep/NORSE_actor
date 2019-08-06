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

function double = extractDistribution(NORSEobject)

    % take the last column of the NORSE distribuion field o.f
    double = NORSEobject.f(:,end);
    

end