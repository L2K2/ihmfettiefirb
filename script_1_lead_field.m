%% This example demonstrates computing TMS-induced E-field with a BEM
%
% Copyright (c) 2020, Lari M. Koponen
%
% Contact:  lari.koponen@aalto.fi
% Version: 2020-09-28
%

% Empty the workspace.
clear;

% The first step is to install HBFTMS library and add it to your path.
%
% It is available at https://github.com/MattiStenroos/hbf_tms_realtime
%
% The exact steps required to add it to your path depend on where you
% install it. For the most straightforward case, where you cloned the
% repository to the same folder with this script, you can run command like:
%
% addpath(genpath(fullfile(pwd(),'hbf_tms_realtime-master')));

% The second step is to install HBF library, and then add it to your path.
%
% The library is available at https://github.com/MattiStenroos/hbf_lc_p
%
% Download it. After that, locate and run 'hbf_SetPaths()' to set paths.

%% Load the model geometry, coil model and a set of 35 target coordinates.
%
% For this purpose, you must first download the supplementary dataset from
% the Mendeley Data repository, or more specifically, the named three files
% from that dataset. I have placed them in folder 'input', but edit the
% script as you see fit to cover your use case.
%
% The data is at https://doi.org/10.17632/bnhb3mgkfb.1
%
%   COILS
%
% The coil models are made by uniform sampling of the coil, similarly to
% the reference coil model in (Stenroos & Koponen, 2019). Such coil model
% is very inefficient in terms of the required number of discretisation
% points. As the computation is reasonably fast even with these models, we
% did not update this work to use our new optimised computational coil
% models from the refered work. (A two layer version of the 42 dipole coil
% model should be virtually indistinguishable from the 2,568 dipole
% high-resolution coil model in the file loaded below, whilst having less
% dipoles than even the low-resolution model with 136 dipoles in the file
% loaded below.)
%
%   Stenroos, M., & Koponen, L.M. (2019). Real-time computation of the
%   TMS-induced electric field in a realistic head model.
%   NeuroImage, 203, 116159.
%   https://doi.org/10.1016/j.neuroimage.2019.116159
%
%   COORDINATES
%
% The targets cover a 15 mm by 10 mm area of the scalp, and define a local
% coordinate system for each of those points. For each target, we compute
% the E-field at 10 degree increments in angle.
%
%   MESHES
%
% Contains the surface meshes for all models in (Koponen et al., 2019).
%
%   Koponen, L.M., Stenroos, M., Nieminen, J.O., Jokivarsi, K., Grohn, O.,
%   and Ilmoniemi, R.J. (2019). Individual head models for estimating the
%   TMS-induced electric field in rat brain. bioRxiv.
%   https://doi.org/10.1101/2019.12.23.886861
%

load('input/coils.mat');
load('input/coordinates.mat');
load('input/meshes.mat');

N=length(coordinates);
angles=-90:10:80;

%% Compute all three comparison models
%
% This will take a few minutes on a laptop. Go make yourself coffee or tea.
%
% For hardware, 8 GiB RAM should be enough, but in that case, I recommend
% closing any other applications prior to running this script. With 16 GiB,
% no such limitations should exist. (Although, during the matrix inversions
% your computer might not be very responsive as MATLAB will use 100% of
% your CPU, even for a multi-core system.)
%

mkdir('output');

for model=1:3

    switch model
        case 1 % spineless, 2C*
            filename='output/lead_field_2C_compare_1.mat';
            bmeshes={compare_1.skull,compare_1.torso};
            ci=0.33*[1/50 1]; 
            co=[ci(2:end) 0];
        case 2 % with spine, 2C
            filename='output/lead_field_2C_compare_2.mat';
            bmeshes={compare_2.skull,compare_2.torso};
            ci=0.33*[1/50 1]; 
            co=[ci(2:end) 0];
        case 3 % just the torso, 1C
            filename='output/lead_field_1C_compare.mat';
            bmeshes={compare_1.torso};
            ci=0.33; 
            co=0;
    end

    % Compute boundary potentials for the brain mesh vertices
    % ... double-layer operators for BEM
    D=hbf_BEMOperatorsPhi_LC(bmeshes);
    % ... BEM transfer matrix
    Tphi_full=hbf_TM_Phi_LC(D,ci,co,length(bmeshes));
    % ... potential matrix
    Phi=single(hbf_LFM_Phi_LC(bmeshes,Tphi_full,brain.p));
    clear D Tphi_full

    % Compute weighted potentials of (Stenroos & Koponen, 2019)
    % Eqs. 4-5
    Phiw=hbftms_WeightedPhi(bmeshes,Phi,ci,co);
    clear Phi

    % Compute dipoles for fast beta integrals of (Stenroos & Koponen, 2019)
    % Eq. 6-7
    [betaQpos,betaQmom]=hbftms_BetaQDipoles(bmeshes);

    % First axis = long axis of figure-of-eight coil
    % Second axis = short axis of figure-of-eight coil
    % Third axis = normal axis
    coiltemplate=coils_HR;
    coil=coiltemplate;

    Ns=3*size(brain.p,1);
    Np=N*length(angles);

    lead_field=zeros(Np,Ns);
    index=0;
    for i=1:N
        position=coordinates(i).origin;
        rotation=coordinates(i).axes;
        for angle_d=angles
            index=index+1;

            % Compute coil orientation, rotate the coil around the target
            % coordinate axis by the desired angle in 'angle_d'
            ROTATION=[cosd(angle_d) sind(angle_d) 0;-sind(angle_d) cosd(angle_d) 0;0 0 1]*rotation;
            % Apply rotation and translation
            coil.QP=coiltemplate.QP*ROTATION+position;
            coil.QN=coiltemplate.QN*ROTATION;

            % Compute the field for given coil
            Bs=Phiw*hbftms_BetaQ(coil.QP,coil.QN,coil.QW,betaQpos,betaQmom);
            Bp=hbftms_BpFlux_xyz(coil.QP,coil.QN,coil.QW,brain.p);
            lead_field(index,:)=Bp+Bs;
        end
    end

    save(filename,'lead_field');

end
