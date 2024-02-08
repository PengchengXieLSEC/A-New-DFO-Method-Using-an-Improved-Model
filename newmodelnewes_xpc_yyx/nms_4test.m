function [x, f, exitflag, output] = nms_4test(fun, x0, options,~)

if (~exist('nms'))
    warning('nms seems unavailable on your computer. Forget about testing it.');
    x = x0;
    f = feval(x0);
else
    
    [x,f,nf]=nms(fun,x0,10^(-6));
    exitflag=0;
    output.funcCount=nf;
end
return;
