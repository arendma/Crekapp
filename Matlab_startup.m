%Startup file to add dependencies to environment


if isunix
    %%%% EDIT FOR YOUR UNIX INSTALLATION %%%%%%
    COBRA_PATH='~/Software/cobratoolbox'; % Path of cobratoolbox
    RAVEN_PATH='~/Software/RAVEN'; %Path of Raven toolbox
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %add paths to matlab path
    for p={COBRA_PATH,'Program/','Program/deps/NIDLE', 'Program/deps/utilities'}
        ckadpath(p{:})
    end
    initCobraToolbox(false)
    run(fullfile(RAVEN_PATH, 'installation/checkInstallation.m'))
    
elseif ispc
    for p={'Program/','Program/deps/NIDLE', 'Program/deps/utilities'}
        ckadpath(p{:})
    end
else
    error("OS not supported")
end

function ckadpath(p, recursive)
%Function to check if a certain directory is already included in the matlab
%path variable and if not add it to the path
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
