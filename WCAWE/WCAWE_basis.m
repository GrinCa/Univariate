function Wtrans = WCAWE_basis(listLHS,coeff_deriv,RHSderiv,nvecfreq)

disp('************************');
disp('Calculating WCAWE basis');
disp('************************');


ndof = size(listLHS{1},1);
nmatglob = length(listLHS);

%warning off all

% ndof = 100;
% nvecfreq = 4;

Vtilde = zeros(ndof,nvecfreq);
V = zeros(ndof,nvecfreq);
U = [];



%--------------------------------------------------------------------------
%Initialisation of random cell matrices
%--------------------------------------------------------------------------


% RHSderiv = cell(1,4);
% Kglob = rand(ndof);
% Mglob = rand(ndof);
% listLHS = {Kglob,Mglob};
% Aglob = rand(ndof,ndof);
% nmatglob = length(listLHS);
% coeff_deriv = rand(2,nvecfreq);

% for ii=1:nvecfreq
%     Vtilde(:,ii) = zeros(ndof,1);
%     V(:,ii) = zeros(ndof,1);
%     %RHSderiv{ii} = rand(ndof,1);
% end

Aglob = sparse(ndof,ndof);
for kk=1:nmatglob
    Aglob = Aglob + coeff_deriv(kk,1)*listLHS{kk};
end

  

% RHSderiv{1,1} = ones(ndof,1);


%--------------------------------------------------------------------------
% MDWCAWE Algorithm
%--------------------------------------------------------------------------

id = initmumps();
id = zmumps(id);
id.JOB = 4;
id = zmumps(id,Aglob);

id.JOB = 3;
id = zmumps(id,Aglob);
Vtilde(:,1) = id.SOL;
U(1,1) = norm(Vtilde(:,1));
V(:,1) = Vtilde(:,1)/U(1,1);    % normalization

for nn=2:nvecfreq
    term1 = 0;
    term2 = 0;
    for mm=1:nn-1
        PU1 = Pu1(nn,mm,U);
        term1 = term1 + RHSderiv{mm+1}*PU1(1,nn-mm);
    end% mm
    for mm=2:nn-1
        PU2 = Pu2(nn,mm,U);
        A_m = sparse(ndof,ndof);
        for kk=1:nmatglob
            A_m = A_m + coeff_deriv(kk,mm+1)*listLHS{kk};
        end
        term2 = term2 + A_m*V(:,1:nn-mm)*PU2(:,nn-mm);
    end% mm
    A_1 = sparse(ndof,ndof);
    for kk=1:nmatglob
        A_1 = A_1 + coeff_deriv(kk,2)*listLHS{kk};
    end
    id.RHS = term1-term2-A_1*V(:,nn-1);
    id = zmumps(id,Aglob);
    Vtilde(:,nn) = id.SOL;
    for alpha=1:nn-1
        U(alpha,nn) = V(:,alpha)'*Vtilde(:,nn);
        Vtilde(:,nn) = Vtilde(:,nn) - U(alpha,nn)*V(:,alpha);
    end %alpha
    U(nn,nn) = norm(Vtilde(:,nn));
    V(:,nn) = Vtilde(:,nn)/U(nn,nn);
end %nn

% delete Mumps instance
id.JOB = -2;
id = zmumps(id);

Wtrans = V;
Ucoeff = U;


function Pu = Pu1(nn,mm,U)
    Pu = eye(nn-mm);
    for t=1:mm
        Pu = Pu/U(t:nn-mm+t-1,t:nn-mm+t-1);
    end
end

function Pu = Pu2(nn,mm,U)
    Pu = eye(nn-mm);
    for t=2:mm
        Pu = Pu/U(t:nn-mm+t-1,t:nn-mm+t-1);
    end
end


end





