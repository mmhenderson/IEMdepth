function h = prepFigForExport(h, pubVsPoster)
% A tool to set several figure properties to standard preferences (e.g.
% axis line width = 1, axis color is black, etc.). They're slightly
% different for publications vs. posters.
%
% VAV 5/15/18

if nargin == 0
    h = gcf;
    pubVsPoster = 1; % assumes you're making publication figures, not poster
end
if nargin == 1
    pubVsPoster = 1;
end

if pubVsPoster
    % You're making publication figures
    fontSize = 10;
else
    % You're making poster figures
    fontSize = 18;    
    h.PaperOrientation = 'landscape';
end
h.Color = 'w';
% h.FontSize = fontSize;

allAx = h.Children;

%% Properties common to both

for i = 1:numel(allAx)
    if isfield(allAx(i),'XColor')
        allAx(i).XColor = [0 0 0];
        allAx(i).YColor = [0 0 0];
        allAx(i).ZColor = [0 0 0];
    end
    allAx(i).LineWidth = 1;
    allAx(i).FontSize = fontSize;
    allAx(i).FontName = 'Arial';
    if isfield(allAx(i),'XLabel')
        allAx(i).XLabel.Color = [0 0 0];
    end
    if isfield(allAx(i),'YLabel')
        allAx(i).YLabel.Color = [0 0 0];
    end
    if isfield(allAx(i),'Box')
        allAx(i).Box = 'off';
    end
end