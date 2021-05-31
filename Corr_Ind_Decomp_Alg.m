function [Z1,Z2,E1,E2,D1,D2,A1,A2] = Corr_Ind_Decomp_Alg(X1,X2,D1,D2,opts)
% correlated-independent components decomposition algorithm
%
% Problem:
%  Columnwise Decomposition of a pair of matrices (X1 and X2) into their
%  inter-modality dependent components (Z1 = D1*A1 and Z2 = D2*A2 where
%  Supp(A1)=Supp(A2)) and their inter-modality independent components
%  (E1 and E2) where Corr(E1,E2) is minimized in a columnwise manner.
%
% Inputs:
%   X1:               input 1
%   X2:               input 2
%   D1:               dictionary 1
%   D2:               dictionary 2
%   (optionals:)
%   opts.DL           perform dictionary learning (true-false, default true)
%   opts.Numiters     number of decomosition iterations (default 5)
%   opts.k            number of nonzero entries in saprse vectors (default 5)
%   opts.DLiters      number of Dictionary leaning iterations (default 1)
%   opts.rho          optimization penalty term (default 10)
%   opts.print        print the results (true-false, default false)
%
% Outputs
%   Z1,Z2             the correlated components
%   E1,E2             the independent components
%   D1,D2             the dictionaries (updated if opts.DL is true)
%   A1,A2             sparse representation matrices with identical support
%_______________________________________________________________________________________________________
%
% Optimization problem:
%
%   minimize_(A1,A2,D1,D2,E1,E2)   ||((E1 - mean(E1)).*(E2 - mean(E2)))./(std(E1).*std(E2))||_F^2
%                                + rho/2*||D1*A1 + E1 - X1||_F^2
%                                + rho/2*||D2*A2 + E2 - X2||_F^2
%
%                      s.t       ||[A1]_i||_0 <=T , forall i
%                                supp(||A1||_0) = supp(||A2||_0)
%                                ||[D1]_t|| = 1, ||[D2]_t|| = 1, forall t
%
%________________________________________________________________________________________________________
% Alternating optimization steps:
%
% Step 1:
%   {A1,A2,D1,D2} = argmin_{A1,A2,D1,D2}    ||D1*A1 + E1 -  X1||_F^2
%                                       + ||D2*A2 + E2 -  X2||_F^2,
%
%                      s.t       ||[A1]_i||_0 <=T , forall i
%                                supp(||A1||_0) = supp(||A2||_0)
%                                ||[D1]_t|| = 1, ||[D2]_t|| = 1, forall t
%
% Step 2:
%   {E1,E2} = argmin_{E1,E2}    ||((E1 - mean(E1)).*(E2 - mean(E2)))./(std(E1).*std(E2))||_F^2
%                             + rho/2*||D1*A1 + E1 - X1 + U1||_F^2
%                             + rho/2*||D2*A2 + E2 - X2 + U2||_F^2
%

%__________________________________________________________________________
%                          Initialization

[m,np,opts] = Parse_input(opts,X1,D1);
% m: size of columns
% np: number of columns

E1 = zeros(m,np); E2 = E1; % independent components
delta = 1e-6; % threshold for minimum considerable variance (for stabilization)
Eps = 5e-4; % approxiation accuracy threshold
rho = opts.rho; % penalty parameter

k_seq = round(linspace(1,opts.k,opts.Numiters)); % gradual increasing of sparsity (warm start)

%__________________________________________________________________________
%
if opts.print==true
    fprintf('Decomposition ... \n');
    fprintf('iter \t ||DA1+E1-X1|| \t ||DA2+E2-X2|| \t ||E[E1,E2]|| \n');
end

%__________________________________________________________________________
%                               Decomposition iterations
%__________________________________________________________________________

for t = 1:opts.Numiters
    %----------------------------------------------------------------------
    %                            {A,D} update
    if opts.DL == true
        [D1,D2,A1,A2] = SCDL(X1-E1,X2-E2,D1,D2,k_seq(t),Eps,opts.DLiters,1);
        
    else                        % A update only
        [A1,A2] = SOMP_C(X1-E1,X2-E2,D1,D2,k_seq(t),Eps);
    end
    
    Z1 = D1*A1;
    Z2 = D2*A2;
    
    %----------------------------------------------------------------------
    %                              E update
    

    
    X1_ = X1-Z1;
    X2_ = X2-Z2;
    E1 = X1_;
    E2 = X2_;
    
    for Eit = 1:t
        
        M1 = mean(E1,1);
        M2 = mean(E2,1);
        V1 = var(E1,1);
        V2 = var(E2,1);
        
        V12 = V1.*V2;
        
        inds = V12>delta;
        C1 = ((E2(:,inds)-M2(:,inds)).^2)./(V12(:,inds));
        E1(:,inds) = (C1.*M1(:,inds) + rho*(X1_(:,inds)))./(C1 + rho);
        
        C2 = ((E1(:,inds)-M1(:,inds)).^2)./(V12(:,inds));
        E2(:,inds) = (C2.*M2(:,inds) + rho*(X2_(:,inds)))./(C2 + rho);
        
    end
    
    %----------------------------------------------------------------------
    %                            Progress
    
    if opts.print ==true
        [f1val,f2val,f3val] = ObjFuncVal(X1,Z1,E1,X2,Z2,E2,delta);
        fprintf('%3g \t %10.3e \t %10.3e \t %10.3e \n', t, f1val,f2val,f3val);
    end
    
end

end


function [m,np,opts] = Parse_input(opts,X,D)

if ~isfield(opts,'print')
    opts.print = false;
end
if ~isfield(opts,'DL')
    opts.DL = true;
end
if ~isfield(opts,'DLiters')
    opts.DLiters = 1;
end
if ~isfield(opts,'Numiters')
    opts.Numiters = 5;
end
if ~isfield(opts,'k')
    opts.k = 5;
end
if ~isfield(opts,'rho')
    opts.rho = 10;
end

[m,~] = size(D);
np = size(X,2);

end

function [f1val,f2val,f3val] = ObjFuncVal(X1,Z1,E1,X2,Z2,E2,delta)
N = numel(X1);
f1val = norm(X1-Z1-E1,'fro')^2/N; % first approximation error
f2val = norm(X2-Z2-E2,'fro')^2/N; % second approximation norm
V1 = var(E1,1);
V2 = var(E2,1);
V12 = V1.*V2;
V12(V12<delta)=delta;
f3val = norm(((E1- mean(E1)).*(E2-mean(E2)))./sqrt(V12),'fro')^2/N;
end
