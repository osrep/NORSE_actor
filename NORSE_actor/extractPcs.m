%%% extractPcs.m
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

function double = extractPcs(NORSEobject)

    % take the mask 
    double = NORSEobject.runawayRegion.pcs(:,1);
    

end