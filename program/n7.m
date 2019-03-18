clear
stem= [0.451471758	0.232299125	0.069809069	0.056682578	0.08194113	0.10779634];
    n=7;
    fid=fopen('prf_test.txt','w');

    vmin=1;
for r = -1:0.01:-0.02
    if (r==0) 
        x = 1:1:n-1;    
        a = sum (1./x);
        y= (1./x)./a;
        y= vpa(y,6)
        v1 = mean((y-stem).^2);
        data=[r v1];
        fprintf(fid,'%d\t',data);
        fprintf(fid,'\n');
     else
        x = 1 : 1: n-1;
      % r
        f = @(x)quad(@(t)(1-exp(-2.*r.*(1-t)))./(1-exp(-2.*r)).*(1./(t.*(1-t)))*(prod(1:n)/(prod(1:x).*prod(1:n-x))).*(t.^x).*((1-t).^(n-x)),0,1);
        b= arrayfun(f,x);
        y1 = sum(b);
        y= b./y1;
        v1 = mean((y-stem).^2);
    
        data=[r v1];
     %data
        if v1 < vmin
            vmin=v1;
            vmin = vpa(vmin,6);
            min = [r vmin];
        end
    fprintf(fid,'%d\t',data);
    fprintf(fid,'\n');
 
   end
end
   fclose(fid);
          min