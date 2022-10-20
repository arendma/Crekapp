%Startup file to add dependencies to environment
%%%% EDIT FOR YOUR INSTALLATION %%%%%%
COBRA_PATH='~/Software/cobratoolbox'; % Path of cobratoolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%add paths to matlab path
ckadpath({COBRA_PATH, "Program/"})

initCobraToolbox(false)

function ckadpath(paths)
cur_path=path();
    for i=1:length(paths)
        if isfolder(paths{i})
            if ~contains(cur_path, paths{i})
                addpath(paths{i})
            end
        else
            error(append(paths{i}, " is not a directory"))
        end
    end
end
