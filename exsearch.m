% [J,selec] = exsearch(f,d,m,fobj,p)
%
% Toolbox: Balu
%    Exhaustive Search Selection for fatures f according to ideal classification d.
%    m features will be selected. fobj = 1 uses Fisher objetctive function, fobj = 0
%    uses Sp @Sn=100%. p is the a priori probability of each class.
%    J is the value of the objective function and selec are the number of the
%    selected features.
%
% D.Mery, PUC-DCC, Apr. 2008
% http://dmery.ing.puc.cl
%

function [J,selec] = exsearch(f,d,m,fobj,p)



if (not(exist('p')))
    dn = max(d)-min(d)+1; % number of classes
    p = ones(dn,1)/dn;
end

M = size(f,2);

N = nchoosek(M,m);
pause(1)

% if (N>10000)
%     ok = input(sprintf('Exhaustive Search needs %d evaluations... continue [y/n]?',N))
%     if (s=='n')
%         error('Exhaustive search for feature selection interrupted.')
%     end
% end

T = nchoosek(1:M,m);

Jmax = 0;
for i=1:N
    fs = f(:,T(i,:));
    if (fobj)
        Js = jfisher(fs,d,p);
    else
        Js = sp100(fs,d);
    end
    if (Js>Jmax)
        ks = i;
        Jmax = Js;
    end
end
selec = T(ks,:);
J = Jmax;