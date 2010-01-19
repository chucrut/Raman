% [J,selec] = sfs(f,d,m,fobj,p)
%
% Toolbox: Balu
%    Sequential Forward Selection for fatures f according to ideal classification d.
%    m features will be selected. fobj = 1 uses Fisher objetctive function, fobj = 0
%    uses Sp @Sn=100%. p is the a priori probability of each class.
%    J is the value of the objective function and selec are the number of the
%    selected features.
%
% D.Mery, PUC-DCC, Apr. 2008
% http://dmery.ing.puc.cl
%

function [J,selec] = sfs(f,d,m,fobj,p)

selec  = []; %selected features
J    = [];
Jmax = 0;
k=1;

dn = max(d)-min(d)+1; % number of classes

if (not(exist('p')))
    p = ones(dn,1)/dn;
end

while (k<=m)
    nuevo = 0;
    for i=1:size(f,2)
        if (k==1) | (sum(selec==i)==0)
            s = [selec i];
            fs = f(:,s);
            if (fobj)
                Js = jfisher(fs,d,p);
            else
                Js = sp100(fs,d);
            end
            if (Js>Jmax)
                ks = i;
                Jmax = Js;
                nuevo = 1;
            end
        end
    end
    if (nuevo)
        selec = [selec ks];
        J = [J Jmax];
        clf
        bar(J)
        hold on
        for i=1:length(selec)
            text(i-0.4,J(i)*1.05,sprintf('%d',selec(i)));
        end
        pause(0)
        k = k + 1;
    else
        disp('no more improvement, the sequantial search is interrupted.');
        k = 1e10;
    end
end
