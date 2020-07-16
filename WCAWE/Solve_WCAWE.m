function [SOLWCAWE] = Solve_WCAWE(listLHS,coeff_deriv_fun,RHS,Wtrans,freq_interval)


nmat_glob = length(listLHS);
ndof = size(listLHS{1},1);

%--------------------------------------------------------------------------
% Project matrices on WCAWE basis and setup RHS
%--------------------------------------------------------------------------
for ii = 1:nmat_glob
   listLHS{ii} = sparse(Wtrans'*listLHS{ii}*Wtrans);
end %ii

%--------------------------------------------------------------------------
% Setup projected RHS
%--------------------------------------------------------------------------
RHS = sparse(Wtrans'*RHS);

%--------------------------------------------------------------------------
% Frequency/Fi loops calculation
%--------------------------------------------------------------------------

nfreq = length(freq_interval);

SOLWCAWE = zeros(ndof,nfreq);


for ii=1:nfreq
      Aglob_red = sparse(size(listLHS{1},1),size(listLHS{1},2));
      for kk = 1:nmat_glob
         Aglob_red = Aglob_red + coeff_deriv_fun{kk,1}(freq_interval(ii))*listLHS{kk};
      end %kk
      solp = Aglob_red\(1i*2*pi*freq_interval(ii)*RHS);
      resp_P = Wtrans*solp;
      SOLWCAWE(:,ii) = resp_P;
end % ii


end

