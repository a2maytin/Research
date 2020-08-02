function  varargout = fitMeanMSDanom(obj, clip_factor,x1,x3)
%%FITMEANMSDANOM Fit the weighted averaged MSD by a nonlinear function.
%
% obj.fitMeanMSDanom computes and fits the weighted mean MSD by a
% nonlinear function. Fit results are displayed in the
% command window.
%
% obj.fitMeanMSD(clip_factor) does the fit, taking into account
% only the first potion of the average MSD curve specified by
% 'clip_factor' (a double between 0 and 1). If the value
% exceeds 1, then the clip factor is understood to be the
% maximal number of point to take into account in the fit. By
% default, it is set to 0.25.
%
% [fo, gof] = obj.fitMeanMSD(...) returns the fit object and the
% goodness of fit.

if nargin < 2
    clip_factor = 0.25;
end

if ~obj.msd_valid
    obj = obj.computeMSD;
end

ft = fittype('poly1');
mmsd = obj.getMeanMSD;

t = mmsd(:,1);
y = mmsd(:,2);
w = 1./mmsd(:,3);

% Clip data, never take the first one dt = 0
if clip_factor < 1
    t_limit = 2 : round(numel(t) * clip_factor);
else
    t_limit = 2 : min(1+round(clip_factor), numel(t)); %Note: I changed 2 to 1 back to 2
end
t = t(t_limit);
y = y(t_limit);
w = w(t_limit);
w = ones([1,length(w)]);  %Note: I got rid of weights

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[fo, gof] = fit(t, y, ft, 'Weights', w); 
%[fo, gof] = fit(t, y, ft);

fun = @(r)(r(1)*t.^r(2)+r(3))-y;
x0 = [x1,1,x3];
x = lsqnonlin(fun,x0);
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

str = sprintf([
    'Estimating D through nonlinear fit of the weighted EATA MSD curve (first two delays).\n', ...
    'D = %.3e'], x(1)/2/obj.n_dim);
disp(str)

if nargout > 0
    varargout{1} = fo;
    if nargout > 1
        varargout{2} = gof;
        if nargout > 2
            varargout{3} = x;
        end
    end
end

end