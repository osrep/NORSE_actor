%%% extractMask.m
%%% 31/07/2019
%%% Written by Soma Olasz
%%% 
%%% This script is created to extract the mask from a 
%%% NORSE calulcation used to determine the runaway
%%% region in momentum space.
%%% 
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