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
flag.recalculated = 0; % allow WCAWE and/or FE recalculation
flag.calculateFE = 1;  % calculate FE solution
flag.calculateWCAWE = 1; % calculate WCAWE solution
flag.cat_basis = 0; % concatenate the sub basis from WCAWE to form one by SVD


flag.plot_solution2D = 0; % plot solution (interpolation) for 2D problem
flag.plot_solution3D = 0; % plot solution (interpolation) for 2D problem

flag.plotcomparison = 0; % plot comparison between FE and WCAWE
flag.comparisonMULTI = 0;

flag.convert2VTK = 0; % convert Pout.mat into a .vtk file
flag.converge_sizemesh = 1;
flag.compare_FE_WCAWE = 0;

flag.get_matrices = 1;
if flag.convert2VTK || flag.converge_sizemesh || flag.compare_FE_WCAWE
    flag.get_matrices = 0;
    flag.rerun = 0;
    flag.recalculated = 0;
end

% define timing struct to access to time calculation of each method                                                    
timing.freefem = 0;
timing.WCAWE = 0;
timing.FE = 0;

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
param.fmin = 100;
param.fmax = 200;
param.f_range = [param.fmin param.fmax];
param.freqincr = 0.5; % 20
param.freq = param.fmin:param.freqincr:param.fmax; % frequency range
param.nfreq = length(param.freq);

% those frequencies are the frequencies point for Padé expension
param.freqref = [120];    
param.nfreqref = length(param.freqref); 

% interval_construct enables us to build sub basis for WCAWE by
% by using the ref frequencies as we want. We can choose which ref freq to
% add for each sub basis. The number of sub basis for WCAWE before SVD is
% equal to length(interval_construct)
param.interval_construct = {[1]};

% Input data for the loop over expansion orders. Note that for each
% frequency sweep the number of vectors from the WCAWE basis will depend on
% the number of point for Padé expension. For instance, if we have 2 points
% for expansion, and nvecfreq=5 (order of expansion), we will have 15
% vectors in the basis, 5 by intervals.
param.nvecfreqmin = 20;
param.nvecfreqmax = 20;
param.incrvec = 5;
param.vecfreqrange = param.nvecfreqmin:param.incrvec:param.nvecfreqmax;

% path to store data
path1 = ['[' num2str(param.f_range(1)) '_' num2str(param.f_range(2)) ']'];

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
end

       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag.recalculated
    
    coeff_LHS = {@(f) 1, @(f) -(2*pi*f/c0)^2};
    coeff_RHS = @(f) 1i*2*pi*f;
    
    % RHS
    RHS_amp = zeros(ndof,1);
    RHS_amp(id_source) = Qc; % RHS = i*omega*RHS
    
    %--------------------------------------------------------------------------
    % Calculate reference frequency response if problem updated, otherwise load
    %--------------------------------------------------------------------------
    if flag.calculateFE == 1
        tic;
        disp('***********************');
        disp('*Calculate FE solution*');
        disp('***********************');
        SOLFE = zeros(ndof,param.nfreq); %size ndof,nfreq
        id = initmumps();
        id = zmumps(id);
        id.JOB = 1;
        id = zmumps(id,LHS{1}+LHS{2});
        % Frequency loop calculation
        for ii=1:param.nfreq
           disp(['[FE] Frequency = ',num2str(param.freq(ii))]);
           tic;
           Aglob = sparse(size(LHS{1},1),size(LHS{1},2));
           for kk = 1:nLHS
              Aglob = Aglob + coeff_LHS{kk}(param.freq(ii))*LHS{kk};
           end %kk
           id.JOB = 5;
           id.RHS = coeff_RHS(param.freq(ii))*RHS_amp;
           id = zmumps(id, Aglob);
           resp_P = id.SOL;
           SOLFE(:,ii) = resp_P;
           toc;
        end % ii
        id.JOB = -2; id = zmumps(id);
        timing.FE = toc;
        ouput = sprintf('[FE] CPUtime for FE frequency sweep: %.4f s',timing.FE);
        disp(ouput);
    end
    
    
    for nvecfreq=param.vecfreqrange

        % Initialize array of RHSderiv, Cell array linear comb coefficients derivative functions
        deriv_deg = [param.nvecfreqmax];

        if exist('Derivatives/derivative_orders.mat','file') ~= 2
            disp('***********************************');
            disp('*Recalculate all cross derivatives*');
            disp('***********************************');
            create_cross_derivatives(FEmatrices,coeff_LHS,...
                                     coeff_RHS,deriv_deg,'f');
        end
        load('Derivatives/derivative_orders.mat');
        if ~isempty(find(derivative_orders-deriv_deg<0))
            disp('***********************************');
            disp('*Recalculate all cross derivatives*');
            disp('***********************************');
            create_cross_derivatives(FEmatrices,coeff_LHS,...
                                     coeff_RHS,deriv_deg,'f');
        end 

        LHScoeffderiv = cell(nLHS,nvecfreq+1);
        RHScoeffderiv = cell(1,nvecfreq+1);
        [LHScoeffderiv,RHScoeffderiv] = get_coeff_deriv_matrices(...
                                LHScoeffderiv,RHScoeffderiv,nvecfreq,nLHS);

        coeff_derivgen_fun = @(freq) cellfun(@(cellfunc) cellfunc(freq),LHScoeffderiv);

        
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
            timing.WCAWE = cputime;
            tic;
            BASIS = cell(1,length(param.interval_construct));
            for ii=1:length(BASIS)
                BASIS{ii} = build_basis(LHS,coeff_deriv_multi,RHSderivmulti,nvecfreq,param.interval_construct{ii});
            end
            timing.WCAWE = toc;
            outputdisplay = sprintf('[WCAWE] CPUtime for building of WCAWE basis : %.4f s',timing.WCAWE);
            disp(outputdisplay);
            
            SOLWCAWE = cell(1,length(param.n_sub_range));
            if flag.cat_basis
                disp("[WCAWE] Compute SVD for concatenation of all sub basis");
                BASIS = [BASIS{:}];
                [uu,vv,~] = svd(BASIS,0);
                iiselect = find(diag(vv)>vv(1,1)*1e-15);
                BASIS_cat = uu(:,iiselect);
                nsvd = size(BASIS_cat,2);
                output = sprintf("[SVD:Info] Number of selected vector %d/%d",nsvd,size(BASIS,2));
                disp(output);
                disp("[WCAWE] Solution phase");
                SOLWCAWE = Solve_WCAWE(LHS,LHScoeffderiv,RHS_amp,BASIS_cat,param.freq);
            else
                disp("[WCAWE] Solution phase");
                for ii=1:param.n_sub_range
                    SOLWCAWE{ii} = Solve_WCAWE(LHS,LHScoeffderiv,RHS_amp,BASIS{ii},param.sub_range{ii});
                end
                SOLWCAWE = [SOLWCAWE{:}];
            end
            
        end

        %--------------------------------------------------------------------------
        % Saves
        %--------------------------------------------------------------------------
        if flag.calculateFE
            save(['Matrices/',mesh.file,'/',path1,'/SOLFE','_sizemesh_',num2str(sizemesh),'.mat'],'SOLFE');
        end
        if flag.calculateWCAWE
            path2 = ['[' num2str(param.f_range(1)) '_' num2str(param.f_range(2)) ']',...
                    '/[' num2str(nvecfreq) '][',replace(num2str(param.freqref),' ','_') ']'];
            save(['Matrices/',mesh.file,'/',path2,'/SOLWCAWE','_sizemesh_',num2str(sizemesh),'.mat'],'SOLWCAWE');
        end
        % Save data only for FE solution
        save(['Matrices/',mesh.file,'/',path1,'/','DATA_sizemesh_',num2str(sizemesh),'.mat'],'FEmatrices','param','timing');

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


if flag.converge_sizemesh
    clear FEmatrices param SOLFE;
    
    meanFE = cell(length(sizemesh_file),1);
    title_VALUES_1 = cell(length(sizemesh_file),1);
    for ii=1:length(sizemesh_file)
        DATA = struct2cell(load(['Matrices/',mesh.file,'/',path1,'/','DATA_sizemesh_',num2str(sizemesh_file(ii)),'.mat']));
        FEmatrices = DATA{1};
        param = DATA{2};
        title_VALUES_1{ii} = [num2str(size(FEmatrices.Nodes,1)) ' ndofs'];
        SOLFE = struct2cell(load(['Matrices/',mesh.file,'/',path1,'/SOLFE','_sizemesh_',num2str(sizemesh_file(ii)),'.mat']));
        SOLFE = SOLFE{1};
        meanFE{ii} = mean(real(SOLFE(FEmatrices.surf_nodes,:)),1);
    end

    arg.config = 'converge';
    arg.xlabel = 'freq';
    arg.ylabel = 'mean pressure (Pa)';
    arg.title = '';
    arg.save_name = ['Convergence FE ' replace(num2str(sizemesh_file'),' ','_')];
    arg.save_path = path1;
    show_graph(arg,meanFE,title_VALUES_1,param);
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
        VALUES_2{ii} = abs((mean(real(SOLWCAWE(FEmatrices.surf_nodes,:)),1) - mean(real(SOLFE(FEmatrices.surf_nodes,:)),1))./mean(1+real(SOLFE(FEmatrices.surf_nodes,:)),1));
        title_VALUES_2{ii} = ['WCAWE [' num2str(param.vecfreqrange(ii)) ' vec]'];
    end
    arg1.config = 'converge';
    arg1.xlabel = 'freq';
    arg1.ylabel = 'mean pressure (Pa)';
    arg1.title =  '';
    arg1.save_name = ['Comparison FE WCAWE ' replace(num2str(param.vecfreqrange),' ','_')];
    arg1.save_path = ['/Matrices/' mesh.file '/' path1];
    arg2.config = 'converge';
    arg2.xlabel = 'freq';
    arg2.ylabel = 'relative error';
    arg2.title =  '';
    arg2.save_name = ['Relative error comparison' replace(num2str(param.vecfreqrange),' ','_')];
    arg2.save_path = ['/Matrices/' mesh.file '/' path1];
    show_graph(arg1,VALUES_1,title_VALUES_1,param);
    show_graph(arg2,VALUES_2,title_VALUES_2,param);
end

%--------------------------------------------------------------------------
% convert
%--------------------------------------------------------------------------
if flag.convert2VTK
    clear FEmatrices;
    sizemesh_ARRAY = load('sizemesh.txt');
    sizemesh = sizemesh_ARRAY(end);
    DATA = struct2cell(load(['Matrices/',mesh.file,'/',path1,'/','DATA_sizemesh_',num2str(sizemesh),'.mat']));
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
        SOLFE = struct2cell(load(['Matrices/',mesh.file,'/',path1,'/SOLFE','_sizemesh_',num2str(sizemesh),'.mat']));
        SOLFE = SOLFE{1};
        convertGEO2VTK(FEmatrices,mesh,sizemesh,SOLFE,PARTITION,param,range);
    end
    if flag.calculateWCAWE
        SOLWCAWE = struct2cell(load(['Matrices/',mesh.file,'/',path2,'/SOLWCAWE','_sizemesh_',num2str(sizemesh),'.mat']));
        SOLWCAWE = SOLWCAWE{1};
        convertGEO2VTK(FEmatrices,mesh,sizemesh,SOLWCAWE,PARTITION,param,range);
    end
end



