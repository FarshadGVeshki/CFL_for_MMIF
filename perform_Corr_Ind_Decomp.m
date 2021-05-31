function [Iz1,Iz2,Ie1,Ie2,D1,D2,A1,A2] = perform_Corr_Ind_Decomp(I1,I2,D1,D2,opts)
% Decomposition of a pair of images (I1 and I2) into their inter-modality 
% dependent components (Iz1 and Iz2) and their inter-modality independent
% components (Ie1 and Ie2), via coupled feature learning and independent
% components approximation.
%
% Inputs:
%   I1:               input 1
%   I2:               input 2
%   D1:               dictionary 1
%   D2:               dictionary 2
%   (optionals:)
%   opts.DL           perform dictionary learning (true-false, default true)
%   opts.Numiters     number of decomosition iterations (default 5)
%   opts.k            number of nonzero entries in saprse vectors (default 10)
%   opts.DLiters      number of Dictionary leaning iterations (default 5) 
%   opts.rho          optimization penalty term (default 10)
%   opts.ss           sliding step for patch extaction (default 1)
%   opts.print        print the results (true-false, default false) 
%   opts.plot         plot the decomposed components (true-false, default false) 
%
% Outputs
%   Iz1,Iz2           the correlated components
%   Ie1,Ie2           the independent components
%   D1,D2             the dictionaries (updated if opts.DL is true)
%   A1,A2             sparse representation matrices with identical support

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isfield(opts,'ss')
    opts.ss = 1;
end

if ~isfield(opts,'plot')
    opts.plot = false;
end

if ~isfield(opts,'delta')
    opts.delta = 1e-6;
end

if ~isfield(opts,'print_res')
    opts.print_res = false;
end

d = opts.delta;
bb = sqrt(size(D1,1));
[nr,nc] = size(I1);


X1 = mexExtractPatches(I1,bb,opts.ss);
X2 = mexExtractPatches(I2,bb,opts.ss);

[Z1,Z2,E1,E2,D1,D2,A1,A2] = Corr_Ind_Decomp_Alg(X1,X2,D1,D2,opts);

Ie1 = mexCombinePatches(E1,zeros(nr, nc),bb,0,opts.ss);
Ie2 = mexCombinePatches(E2,zeros(nr, nc),bb,0,opts.ss);

Iz1 = [];
Iz2 = [];

if opts.print_res == true
V1 = std(E1,1);
V2 = std(E2,1);
V12 = V1.*V2;
inds = (V12.^2>d);

C = abs(((E1(:,inds)- mean(E1(:,inds))).*(E2(:,inds)-mean(E2(:,inds))))./sqrt(V12(:,inds)));
Fvals(2) = mean(C(:));
Iz1 = mexCombinePatches(Z1,zeros(nr, nc),bb,0,opts.ss);
Iz2 = mexCombinePatches(Z2,zeros(nr, nc),bb,0,opts.ss);
Fvals(1) = immse([I1 I2], [Iz1+Ie1 Iz2+Ie2]);
fprintf('Decomposition results: \n');
fprintf(' avg. MSE \t avg. abs corr \n');
fprintf('%10.3e \t %10.3e \n', Fvals(1),Fvals(2));
end



if opts.plot == true
    Iz1 = mexCombinePatches(Z1,zeros(nr, nc),bb,0,opts.ss);
    Iz2 = mexCombinePatches(Z2,zeros(nr, nc),bb,0,opts.ss);
    figure(111), colormap(gray)
    
    subplot(231);
    imagesc(I1),colorbar
    ax = gca;
    axis(ax,'off')
    xlabel(ax,'I_1')
    ax.XLabel.Visible = 'on';
    
    subplot(232);
    imagesc(Iz1),colorbar
    ax = gca;
    axis(ax,'off')
    xlabel(ax,'I^z_1')
    ax.XLabel.Visible = 'on';
    
    subplot(233);
    imagesc(Ie1),colorbar
    ax = gca;
    axis(ax,'off')
    xlabel(ax,'I^e_1')
    ax.XLabel.Visible = 'on';
    
    
    subplot(234);
    imagesc(I2),colorbar
    ax = gca;
    axis(ax,'off')
    xlabel(ax,'I_2')
    ax.XLabel.Visible = 'on';    
    
    subplot(235);
    imagesc(Iz2),colorbar
    ax = gca;
    axis(ax,'off')
    xlabel(ax,'I^z_2')
    ax.XLabel.Visible = 'on';
        
    subplot(236);
    imagesc(Ie2),colorbar
    ax = gca;
    axis(ax,'off')
    xlabel(ax,'I^e_2')
    ax.XLabel.Visible = 'on';
    
    
end

end
