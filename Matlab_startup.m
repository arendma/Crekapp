%Startup file to add dependencies to environment

%%%% EDIT FOR YOUR INSTALLATION %%%%%%
COBRA_PATH='~/Software/cobratoolbox'; % Path of cobratoolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%add paths to matlab path
ckadpath(COBRA_PATH)
ckadpath('Program/', true)

initCobraToolbox(false)

function ckadpath(p, recursive)
if nargin<2
    recursive=false;
end
cur_path=path();
        if isfolder(p)
            if ~contains(cur_path, p)
                if recursive
                    addpath(genpath(p))
                else
                    addpath(p)
                end
            end
        else
            error(append(p, " is not a directory"))
        end
end
