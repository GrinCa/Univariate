function Wtranssvd = build_basis(LHS,coeff_deriv_multi,RHSderivmulti,nvecfreq,pade_points)

Wtrans = [];
for ii=1:length(pade_points)
    Wtranstmp = WCAWE_basis(LHS,coeff_deriv_multi{pade_points(ii)},RHSderivmulti{pade_points(ii)},nvecfreq);
    Wtrans = [Wtrans Wtranstmp];
end

[uu,vv,~] = svd(Wtrans,0);
iiselect = find(diag(vv)>vv(1,1)*1e-15);
Wtranssvd = uu(:,iiselect);
nsvd = size(Wtranssvd,2);
output = sprintf("[SVD:Info] Number of selected vector %d/%d",nsvd,size(Wtrans,2));
disp(output);

end