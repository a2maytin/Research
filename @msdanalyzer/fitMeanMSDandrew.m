function  varargout = fitMeanMSDandrew(obj)
%%Fit the weighted averaged MSD by a linear function. 
% Modified version of fitMeanMSD.
%
% obj.fitMeanMSD computes and fits the weighted mean MSD by a
% straight line y = a * x. The fit is therefore valid only for
% purely diffusive behavior. Fit results are displayed in the
% command window.
%
%
% [fo, gof] = obj.fitMeanMSD(...) returns the fit object and the
% goodness of fit.

if ~obj.msd_valid
    obj = obj.computeMSD;
end

ft = fittype('poly1');
mmsd = obj.getMeanMSD;

t = mmsd(:,1);
y = mmsd(:,2);
w = 1./mmsd(:,3);

w(1) = 1;

[fo, gof] = fit(t, y, ft, 'Weights', w);

ci = confint(fo);
str = sprintf([
    'Estimating D through linear weighted fit of the mean MSD curve.\n', ...
    'D = %.3e with 95%% confidence interval [ %.3e - %.3e ].\n', ...
    'Goodness of fit: R² = %.3f.' ], ...
    fo.p1/2/obj.n_dim, ci(1)/2/obj.n_dim, ci(2)/2/obj.n_dim, gof.adjrsquare);
disp(str)

if nargout > 0
    varargout{1} = fo;
    if nargout > 1
        varargout{2} = gof;
    end
end

end