function [D1,D2,A1,A2] = SCDL(X1,X2,D1,D2,K,Eps,Numiters,fixedDCatom)
% Simaltaneous and coupled dictionary learning

[m,na] = size(D1);

if fixedDCatom ==1
    D1(:,1) = 1/sqrt(m)*ones(m,1);
    D2(:,1) = D1(:,1);
end

j=1;
if fixedDCatom ==1
    j =2;
end

Iind = [];

for i = 1: Numiters
    
    [A1,A2] = SOMP_C(X1,X2,D1,D2,K,Eps);   % simaltaneous sparse coding (with coupled dictionaries)
    
    % atom update (K-SVD)
    for ii = j:na
        inds = find(A1(ii,:));
        if ~isempty(inds)
            Dt1 = D1;
            Dt2 = D2;
            Dt1(:,ii) = 0;
            Dt2(:,ii) = 0;
            
            R1 = X1(:,inds) - Dt1*A1(:,inds);
            R2 = X2(:,inds) - Dt2*A2(:,inds);
            [u,s,v] = svds(R1,1);
            D1(:,ii) = u;
            A1(ii,inds) = s*v';
            [u,s,v] = svds(R2,1);
            D2(:,ii) = u;
            A2(ii,inds) = s*v';
            
        else
            R = vecnorm(X1-D1*A1).*vecnorm(X2-D2*A2);
            R(Iind) = 0;
            [~,iind] = max(R);
            Iind = [Iind iind];
            D1(:,ii) = X1(:,iind)/norm(X1(:,iind));
            D2(:,ii) = X2(:,iind)/norm(X2(:,iind));            
        end
    end
    
    
    
    
end





%     RMS1  = sqrt(sum((X1-D1*A1).^2,1:2))/numel(X1);
%     RMS2  = sqrt(sum((X2-D2*A2).^2,1:2))/numel(X2);
%     S1 = nnz(A1);
%     S2 = nnz(A2);

%     fprintf('&&&&&&& %3g \t %2.3e \t\t %3.3e \t\t %3.3e \t\t %3.3e \n', i, RMS1,RMS2,S1,S2);

end
% fprintf('&&&&&&& CDL Done! \n');
% fprintf('________________________________________________\n');



