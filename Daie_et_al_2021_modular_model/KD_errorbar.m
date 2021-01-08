function KD_errorbar(x,y,e,color,marker,markersize,Xer);
if ~iscell(color)
    color = {color};
end
if nargin < 4
    color = {'k'}
end
if nargin < 5
    marker = 'o-';
end
if nargin < 6
    markersize = 5;
end
if nargin < 7
    Xer = 0;
end
if length(color) == 1;
    color{2} = color{1};
end
X = reshape(x,1,[]);
Y = reshape(y,1,[]);
E = reshape(e,1,[]);
Y = [Y+E;Y-E];
if Xer == 0
    plot([X;X],Y,'color',color{1});hold on;
else
    X = reshape(x,1,[]);
    Y = reshape(y,1,[]);
    E = reshape(e,1,[]);
    X = [X+E;X-E];
    plot(X,[Y;Y],'color',color{1});hold on;
end
plot(x,y,marker,'color',color{1},'markersize',markersize,'markerfacecolor',color{2});hold on;
hold off