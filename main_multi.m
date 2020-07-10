%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Test case using WCAWE                             %
%                                                                         %
%                             March 2020                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------------------------------------------------------------
% Init main program
%--------------------------------------------------------------------------

warning('off', 'MATLAB:nearlySingularMatrix');
format short eng;


%--------------------------------------------------------------------------
% Folders
%--------------------------------------------------------------------------

% Add folders for Mumps, WCAWE, and Mesh functions/files
global Meshfold WCAWEfold Mumps DataMap Derivatives
Meshfold = 'Matrices';
WCAWEfold = 'WCAWE';
Mumps = 'Mumps';
DataMap = 'DataMap';
Derivatives = 'Derivatives';

addpath(genpath(strcat(pwd,'/',Meshfold)));
addpath(genpath(strcat(pwd,'/',WCAWEfold)));
addpath(genpath(strcat(pwd,'/',Mumps)));
addpath(genpath(strcat(pwd,'/',DataMap)));
addpath(genpath(strcat(pwd,'/',Derivatives)));

%--------------------------------------------------------------------------
% Input data for the problem
%--------------------------------------------------------------------------

% Input parameters for Matlab calculation
flag.rerun = 0; % to recalculate FreeFem++ matrices
flag.recalculated = 1; % allow WCAWE and/or FE recalculation
flag.calculateFE = 1;  % calculate FE solution
flag.calculateWCAWE = 1; % calculate WCAWE solution

flag.plot_solution2D = 0; % plot solution (interpolation) for 2D problem
flag.plot_solution3D = 0; % plot solution (interpolation) for 2D problem

flag.plotcomparison = 0; % plot comparison between FE and WCAWE
flag.comparisonMULTI = 0;

flag.convert2VTK = 0; % convert Pout.mat into a .ktf file
flag.converge_sizemesh = 0;
flag.compare_FE_WCAWE = 1;

flag.get_matrices = 1;
if flag.convert2VTK || flag.converge_sizemesh || flag.compare_FE_WCAWE
    flag.get_matrices = 0;
    flag.rerun = 0;
    flag.recalculated = 0;
end

% define timing struct to access to time calculation of each method                                                    
timing.freefem = 0;
timing.WCAWE = 0;
timing.computeFE = 0;

mesh.file = 'Truck';
sizemesh_file = load('sizemesh.txt');
sizemesh = sizemesh_file(end);

% Source
source = [0.8,1.15,0.5];  % coordinates of the source 
Qc = 0.01;

% Material parameters
rho = 1.2;
c0 = 340;
P0 = 2e-5;

% Frequency range
param.fmin = 130; % 300
param.fmax = 150; % 6000
param.f_range = [param.fmin param.fmax];
param.freqincr = 0.5; % 20
param.freq = param.fmin:param.freqincr:param.fmax; % frequency range
param.nfreq = length(param.freq);

param.freqref = [140];    % those frequencies are the frequencies point for
                          % Padé expension
param.nfreqref = length(param.freqref); 

% Input data for the loop over expansion orders. Note that for each
% frequency sweep the number of vectors from the WCAWE basis will depend on
% the number of point for Padé expension. For instance, if we have 2 points
% for expansion, and nvecfreq=5 (order of expansion), we will have 15
% vectors in the basis, 5 by intervals.
param.nvecfreqmin = 10;
param.nvecfreqmax = 10;
param.incrvec = 5;
param.vecfreqrange = param.nvecfreqmin:param.incrvec:param.nvecfreqmax;


%--------------------------------------------------------------------------
% Build intervals for the Parametric sweep
%--------------------------------------------------------------------------
param = build_interval(param);

% generate folders needed for computation
genfolders(mesh,param)

%--------------------------------------------------------------------------
% Matrices calculated with Freefem++
%--------------------------------------------------------------------------
if flag.get_matrices
    matrix_names = ["Hr.txt","Hi.txt",...
                    "Qr.txt","Qi.txt"];
                
    [FEmatrices,ndof,timing,flag] = get_matrices(timing,flag,mesh,matrix_names,param);
    Nodes = FEmatrices.Nodes;
    LHS = FEmatrices.LHS;
    nLHS = length(LHS);
    [~,id_source] = min((FEmatrices.Nodes(:,1)-source(1)).^2 + ...
                        (FEmatrices.Nodes(:,2)-source(2)).^2 + ...
                        (FEmatrices.Nodes(:,3)-source(3)).^2);
    % Save data only for FE solution
    save(['Matrices/',mesh.file,'/',folder.path1,'/','DATA_sizemesh_',num2str(sizemesh),'.mat'],'FEmatrices','param','timing');
end

       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag.recalculated
    
    %--------------------------------------------------------------------------
    % Calculate reference frequency response if problem updated, otherwise load
    %--------------------------------------------------------------------------
    t_0 = cputime;
    if flag.calculateFE == 1
    disp('#######################');
    disp('Calculate FE solution');
    disp('#######################');
    SOLFE = zeros(ndof,param.nfreq); %size ndof,nfreq
    id = initmumps();
    id = zmumps(id);
    id.JOB = 1;
    id = zmumps(id,LHS{1}+LHS{2});
    % Frequency loop calculation
    for ii=1:param.nfreq
       tic;
       disp(['[FE] Frequency = ',num2str(param.freq(ii))]);
       Aglob = sparse(size(LHS{1},1),size(LHS{1},2));
       for kk = 1:nLHS
          Aglob = Aglob + LHScoeffderiv{kk,1}(param.freq(ii))*LHS{kk};
       end %kk
       id.JOB = 5;
       id.RHS = RHScoeffderiv{1}(param.freq(ii))*RHS_amp;
       id = zmumps(id, Aglob);
       resp_P = id.SOL;
       SOLFE(:,ii) = resp_P;
       toc;
    end % ii
    id.JOB = -2; id = zmumps(id);
    timing.computeFE = cputime-t_0;
    end
    
    
    for nvecfreq=param.vecfreqrange

        % Initialize array of RHSderiv, Cell array linear comb coefficients derivative functions
        deriv_deg = [param.nvecfreqmax];

        if exist('Derivatives/derivative_orders.mat','file') ~= 2
            disp('#################################');
            disp('Recalculate all cross derivatives');
            disp('#################################');
            coeff_LHS = {@(f) 1, @(f) -(2*pi*f/c0)^2};
            coeff_RHS = @(f) 1i*2*pi*f;
            create_cross_derivatives(FEmatrices,coeff_LHS,...
                                     coeff_RHS,deriv_deg,'f');
        end
        load('Derivatives/derivative_orders.mat');
        if ~isempty(find(derivative_orders-deriv_deg<0))
            disp('#################################');
            disp('Recalculate all cross derivatives');
            disp('#################################');
            coeff_LHS = {@(f) 1, @(f) -(2*pi*f/c0)^2};
            coeff_RHS = @(f) 1i*2*pi*f;
            create_cross_derivatives(FEmatrices,coeff_LHS,...
                                     coeff_RHS,deriv_deg,'f');
        end 

        LHScoeffderiv = cell(nLHS,nvecfreq+1);
        RHScoeffderiv = cell(1,nvecfreq+1);
        [LHScoeffderiv,RHScoeffderiv] = get_coeff_deriv_matrices(...
                                LHScoeffderiv,RHScoeffderiv,nvecfreq,nLHS);

        coeff_derivgen_fun = @(freq) cellfun(@(cellfunc) cellfunc(freq),LHScoeffderiv);

        % RHS
        RHS_amp = zeros(ndof,1);
        RHS_amp(id_source) = Qc; % RHS = i*omega*RHS
        RHSderivmulti = cell(1,param.nfreqref);
        for ii=1:param.nfreqref
            RHSderiv = cell(1,nvecfreq+1);
            for kk=1:nvecfreq+1
                RHSderiv{kk} = RHScoeffderiv{kk}(param.freqref(ii))*RHS_amp;    %derivatives of RHS at freq0
            end
            RHSderivmulti{ii}=RHSderiv;
        end

        % Fill Cell array of linear comb coefficients derivatives at
        % freqref(ii)
        coeff_deriv_multi = cell(1,param.nfreqref);
        for ii=1:param.nfreqref
            coeff_deriv_multi{ii}=coeff_derivgen_fun(param.freqref(ii));
            % Fix coeff_deriv by replacing all NaN values by 0
            tmpidxnan = find(isnan(coeff_deriv_multi{ii}));
            coeff_deriv_multi{ii}(tmpidxnan) = 0;
            tmpidxinf = find(isinf(coeff_deriv_multi{ii}));
            coeff_deriv_multi{ii}(tmpidxinf) = 0;
        end

        %-----------------------------------------------------------------------
        %Recalculation of WCAWE basis
        %-----------------------------------------------------------------------

        if flag.calculateWCAWE
               Wtrans = [];
               for ii=1:param.nfreqref
                  [Wtranstmp,Ucoeff,timing] = WCAWE_basis(LHS,coeff_deriv_multi{ii},RHSderivmulti{ii},nvecfreq,timing);
                   Wtrans = [Wtrans Wtranstmp];
               end
               [uu,vv,ww] = svd(Wtrans,0);
               iiselect = find(diag(vv)>vv(1,1)*1e-15);
               Wtranssvd = uu(:,iiselect);
               nsvd = size(Wtranssvd,2);
               output = sprintf("[SVD:Info] Number of selected vector %d/%d",nsvd,size(Wtrans,2));
               disp(output);
               [SOLWCAWE] = Solve_WCAWE(LHS,LHScoeffderiv,RHS_amp,Wtranssvd,param.freq);
        end

        %--------------------------------------------------------------------------
        % Saves
        %--------------------------------------------------------------------------
        if flag.calculateFE
            path1 = ['[' num2str(param.f_range(1)) '_' num2str(param.f_range(2)) ']'];
            save(['Matrices/',mesh.file,'/',path1,'/SOLFE','_sizemesh_',num2str(sizemesh),'.mat'],'SOLFE');
        end
        if flag.calculateWCAWE
            path2 = ['[' num2str(param.f_range(1)) '_' num2str(param.f_range(2)) ']',...
                    '/[' num2str(nvecfreq) '][',replace(num2str(param.freqref),' ','_') ']'];
            save(['Matrices/',mesh.file,'/',path2,'/SOLWCAWE','_sizemesh_',num2str(sizemesh),'.mat'],'SOLWCAWE');
        end

    end

end


if flag.plot_solution2D == 1
    Nodes = load(strcat('Matrices/',edp.file,'/','Nodes_PML_',edp.file,'.txt'));
    nvec_basis = num2str(input('Choose the number of vector of the basis'));
    SOLWCAWE = struct2cell(load(strcat('Matrices/',edp.file,'/','PoutWCAWE_',edp.file,'_nvec',...
                                 nvec_basis,'_[',num2str(f_range),']','.mat')));
    plotsurf2D(mesh,Nodes,real(SOLWCAWE{1}(:,20)));
end

if flag.plot_solution3D == 1
    Nodes = load(strcat('Matrices/',edp.file,'/','Nodes_PML_',edp.file,'.txt'));
    nvec_basis = num2str(input('Choose the number of vector of the basis : '));
    SOLWCAWE = struct2cell(load(strcat('Matrices/',edp.file,'/','PoutWCAWE_',edp.file,'_nvec',...
                                 nvec_basis,'_[',num2str(f_range),']','.mat')));
    plotsurf3D(mesh,Nodes,real(SOLWCAWE{1}(:,1)),3,0);
    disp('[main_simulation] Warning | if the plot is blank, change to plotsurf2D');
end


if flag.comparisonMULTI
    SOLFE = cell(1,4);
    nvecbasis=40;
    PoutFE = struct2cell(load(strcat('Matrices/',edp.file,'/','[',num2str(f_range),']/','Pout_',edp.file,'.mat')));
    SOLFE{1} = PoutFE{1};
    Pout1 = struct2cell(load(strcat('Matrices/',edp.file,'/','[',num2str(f_range),']/','PoutWCAWE_',edp.file,'_nvec',num2str(nvecbasis),'_[',num2str([300]),']','.mat')));
    SOLFE{2} = Pout1{1};
    Pout2 = struct2cell(load(strcat('Matrices/',edp.file,'/','[',num2str(f_range),']/','PoutWCAWE_',edp.file,'_nvec',num2str(nvecbasis),'_[',num2str([250 350]),']','.mat')));
    SOLFE{3} = Pout2{1};
    Pout4 = struct2cell(load(strcat('Matrices/',edp.file,'/','[',num2str(f_range),']/','PoutWCAWE_',edp.file,'_nvec',num2str(nvecbasis),'_[',num2str([225 275 325 375]),']','.mat')));
    SOLFE{4} = Pout4{1};
    hold on
    plot(freq,real(SOLFE{1}(obs_point,:)),'+','DisplayName',['FE']);
    for ii=2:4
        plot(freq,real(SOLFE{ii}(obs_point,:)),'DisplayName',['NbPoint(Padé): ' num2str(2^(ii-2))]);
    end
    legend;
    hold off
    
    
end

if flag.converge_sizemesh
    clear FEmatrices param SOLFE;
    try
        meanFE = cell(length(sizemesh_file),1);
        title_VALUES_1 = cell(length(sizemesh_file),1);
        for ii=1:length(sizemesh_file)
            DATA = struct2cell(load(['Matrices/',mesh.file,'/',folder.path1,'/','DATA_sizemesh_',num2str(sizemesh_file(ii)),'.mat']));
            FEmatrices = DATA{1};
            title_VALUES_1{ii} = [num2str(size(FEmatrices.Nodes,1)) ' ndofs'];
            SOLFE = struct2cell(load(['Matrices/',mesh.file,'/',folder.path1,'/SOLFE','_sizemesh_',num2str(sizemesh_file(ii)),'.mat']));
            SOLFE = SOLFE{1};
            meanFE{ii} = mean(real(SOLFE(FEmatrices.surf_nodes,:)),1);
        end
    catch
         disp("[COMPARISON] : error in the file");
         meanFE = {zeros(length(param.freq),1)};
         title_VALUES_1 = {''};
    end
    arg.config = 'converge';
    arg.xlabel = 'freq';
    arg.ylabel = 'mean pressure (Pa)';
    arg.title = '';
    arg.save_name = ['Convergence FE ' replace(num2str(sizemesh_file'),' ','_')];
    show_graph(arg,meanFE,title_VALUES_1,mesh,param);
end


if flag.compare_FE_WCAWE
    path1 = ['[' num2str(param.f_range(1)) '_' num2str(param.f_range(2)) ']'];
    DATA = struct2cell(load(['Matrices/',mesh.file,'/',path1,'/','DATA_sizemesh_',num2str(sizemesh),'.mat']));
    FEmatrices = DATA{1};
    
    title_VALUES_1 = cell(length(param.vecfreqrange) + 1,1); % for FE solution + all WCAWE solution
    title_VALUES_2 = cell(length(param.vecfreqrange),1); % for FE solution + all WCAWE solution
    VALUES_1 = cell(length(param.vecfreqrange) + 1,1);
    VALUES_2 = cell(length(param.vecfreqrange),1);
    
    %read FE solution
    SOLFE = struct2cell(load(['Matrices/',mesh.file,'/',path1,'/SOLFE','_sizemesh_',num2str(sizemesh),'.mat']));
    SOLFE = SOLFE{1};
    VALUES_1{1} = mean(real(SOLFE(FEmatrices.surf_nodes,:)),1);
    title_VALUES_1{1} = 'FE';

    %read WCAWE solution
    for ii=1:length(param.vecfreqrange)
        path2 = ['[' num2str(param.f_range(1)) '_' num2str(param.f_range(2)) ']',...
                '/[' num2str(param.vecfreqrange(ii)) '][',replace(num2str(param.freqref),' ','_') ']'];
        SOLWCAWE = struct2cell(load(['Matrices/',mesh.file,'/',path2,'/SOLWCAWE','_sizemesh_',num2str(sizemesh),'.mat']));
        SOLWCAWE = SOLWCAWE{1};
        VALUES_1{ii+1} = mean(real(SOLWCAWE(FEmatrices.surf_nodes,:)),1);
        title_VALUES_1{ii+1} = ['WCAWE [' num2str(param.vecfreqrange(ii)) ' vec]'];
    end
    % calculate relative error
    for ii=1:length(param.vecfreqrange)
        VALUES_2{ii} = abs((mean(real(SOLWCAWE(FEmatrices.surf_nodes,:)),1) - mean(real(SOLFE(FEmatrices.surf_nodes,:)),1))./mean(real(SOLFE(FEmatrices.surf_nodes,:)),1));
        title_VALUES_2{ii} = ['WCAWE [' num2str(param.vecfreqrange(ii)) ' vec]'];
    end
    arg1.config = 'converge';
    arg1.xlabel = 'freq';
    arg1.ylabel = 'mean pressure (Pa)';
    arg1.title = '';
    arg1.save_name = ['Convergence FE ' replace(num2str(sizemesh_file'),' ','_')];
    arg2.config = 'converge';
    arg2.xlabel = 'freq';
    arg2.ylabel = 'relative error';
    arg2.title = '';
    arg2.save_name = ['Convergence FE ' replace(num2str(sizemesh_file'),' ','_')];
    show_graph(arg1,VALUES_1,title_VALUES_1,mesh,param);
    show_graph(arg2,VALUES_2,title_VALUES_2,mesh,param);
end

%--------------------------------------------------------------------------
% convert
%--------------------------------------------------------------------------
if flag.convert2VTK
    clear FEmatrices;
    sizemesh_ARRAY = load('sizemesh.txt');
    sizemesh = sizemesh_ARRAY(end);
    DATA = struct2cell(load(['Matrices/',mesh.file,'/',folder.path1,'/','DATA_sizemesh_',num2str(sizemesh),'.mat']));
    FEmatrices = DATA{1};
    param = DATA{2};
    %%%
    PARTITION = cell(1);
    
    PARTITION{1} = {'scalar',...
                    1:1:FEmatrices.size_system,...
                    1:1:FEmatrices.size_system,...
                    'PRESSURE'};
    
    range = {1:1:1};

    if flag.calculateFE
        SOLFE = struct2cell(load(['Matrices/',mesh.file,'/',folder.path1,'/SOLFE','_sizemesh_',num2str(sizemesh),'.mat']));
        SOLFE = SOLFE{1};
        convertGEO2VTK(FEmatrices,mesh,sizemesh,SOLFE,PARTITION,param,range);
    end
    if flag.calculateWCAWE
        SOLWCAWE = struct2cell(load(['Matrices/',mesh.file,'/',folder.path2,'/SOLWCAWE','_sizemesh_',num2str(sizemesh),'.mat']));
        SOLWCAWE = SOLWCAWE{1};
        convertGEO2VTK(FEmatrices,mesh,sizemesh,SOLWCAWE,PARTITION,param,range);
    end
end



