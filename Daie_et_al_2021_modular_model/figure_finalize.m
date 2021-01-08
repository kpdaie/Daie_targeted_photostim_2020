function figure_finalize(fs);
if nargin == 0
    fs = 8;
end
ch = get(gcf,'children');
for i = 1:length(ch);
    try
        ch(i).FontSize = fs;
        ch(i).FontName = 'Arial';
        ch(i).Box = 'off';
        ch(i).TickDir = 'out';
    end
end