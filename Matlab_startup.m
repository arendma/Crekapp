%Startup file to add dependencies to environment
%%%% EDIT FOR YOUR INSTALLATION %%%%%%
COBRA_PATH='~/Software/cobratoolbox'; % Path of cobratoolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add paths to matlab path
ckadpath(COBRA_PATH)
ckadpath('Program/', true)

initCobraToolbox(false)

function ckadpath(p, recursive)
if nargs<2
    recursive=false;
end
cur_path=path();
        if isfolder(path)
            if ~contains(cur_path, path)
                if recursive
                    addpath(genpath(path))
                else
                    addpath(p)
                end
            end
        else
            error(append(p, " is not a directory"))
        end
end
