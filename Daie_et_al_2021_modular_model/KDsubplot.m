function ax=KDsubplot(row,col,num,marg)


if length(num) == 1
    nums(1) = ceil(num/col);
    nums(2) = round(mod(num,col+.1));
    num = nums;
end
if numel(marg)==1;
    marg(2)=marg;
end
if num(1)>row | num(2)>col
    'no'
else
    set(gcf,'units','inches');
    % num=fliplr(num);
    num(1)=row-num(1)+1;
    P=get(gcf,'position');
    A(3)=(P(3)-marg(1)*(col+1))/col;
    A(4)=(P(4)-marg(2)*(row+1))/row;
    
    A(1)=(num(2)-1)*A(3)+marg(1)*num(2);
    A(2)=(num(1)-1)*A(4)+marg(2)*num(1);
    axes
    set(gca,'units',get(gcf,'units'));
    set(gca,'position',A);
end
ax=gca;