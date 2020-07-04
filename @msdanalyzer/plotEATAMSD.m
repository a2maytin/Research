function varargout = plotEATAMSD(obj, ha, indices)
%%PLOTEATAMSD Plot the ensemble average of the MSD curves.
%
% obj,plotEATAMSD computes and plots the ensemble average of all MSD
% curves.
%
% obj,plotEATAMSD(ha) plots the curve in the axes with the
% specified handle.
%
%
% h = obj,plotEATAMSD(...) returns the handle to the line plotted.
%
% [h, ha] = obj,plotEATAMSD(...) also returns the handle of the
% axes in which the curve was plotted.

if nargin < 4
    indices = [];
end

msmsd = obj.getEATAMSD(indices);

if nargin < 3
    errorbar = false;
    if nargin < 2
        ha = gca;
    end
end


h = plot(ha, msmsd(:,1), msmsd(:,2), 'k', ...
        'LineWidth', 2);


obj.labelPlotMSD(ha);

if nargout > 0
    varargout{1} = h;
    if nargout > 1
        varargout{2} = ha;
    end
end

end