%%% createStructure.m
%%% 08/01/2019
%%% Written by Soma Olasz
%%% 
%%% This script is written in order to create a Matlab structure from
%%% arbitrary input variables. The purpose of this file is to be able to
%%% create Matlab structures from Python.
%%% It requires an even number of input arguments. From every two input, 
%%% the first variable is the data, the second is the name of the field of 
%%% the structure. The names must be in string format.
%%%
%%% Usage
%%%
%%%     struct = createStructure(data, fieldName, ...)
%%%


function struct = createStructure(varargin)

    if mod(nargin, 2) ~= 0        
        error('Need an even number of input arguments.');
    end

    for i = 1 : nargin/2
        
        struct.(varargin{2*i}) = varargin{2*i-1}';
        
    end

end