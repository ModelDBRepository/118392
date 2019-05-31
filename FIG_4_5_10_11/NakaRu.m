%this is normalised, M=1


function frate=NakaRu(inp,semisat)

%if inp<=0
%    frate=0
%else
    %frate=inp/semisat;
    frate=((inp+abs(inp))/2).^2./(semisat.^2+inp.^2);
frate(find(isnan(frate)))=0;
%end
end