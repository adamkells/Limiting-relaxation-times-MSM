% function by Adam Kells to fit to relaxation times calculated at different
% lag times

function [x1] =optimize_fit(xdata,mm_all)
factor=1/1000;
mm_all=mm_all*factor;

F=@(x,xdata)xdata.*x(2)./(xdata+x(1)*x(2));
x0(1)=(mm_all(2,1)-mm_all(1,1))/(xdata(2)-xdata(1));
x0(2)=max(mm_all(:,1));

[x1,~,~,~,~] = lsqcurvefit(F,x0,xdata',mm_all);
limit_mm=x1(2);
mm_all=mm_all/factor;
x1(1)=x1(1)*factor;
x1(2)=x1(2)/factor;
