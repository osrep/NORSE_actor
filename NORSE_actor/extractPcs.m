%%% extractPcs.m
%%% 31/07/2019
%%% Written by Soma Olasz
%%% 
%%% This script is created to extract the critical momentum
%%% values from a NORSE calulcation used to determine the 
%%% runaway region in momentum space.
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