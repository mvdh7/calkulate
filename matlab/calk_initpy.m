function calk_initpy(python_exe)
% Check if Python is loaded in MATLAB, and load it if not
% Matthew P. Humphreys
% Written by Matthew P. Humphreys, last updated 2018-12-20

% === INPUT ============
% python_exe: string, path to calkenv Python executable

[~,~,pyloaded] = pyversion;

if ~pyloaded
    pyversion(python_exe);
end %if
    
end %function calk_initpy
