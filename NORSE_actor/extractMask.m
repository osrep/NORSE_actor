%%% extractMask.m
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

function double = extractMask(NORSEobject)

    % take the mask 
    double = NORSEobject.runawayRegion.mask(:,1);
    

end