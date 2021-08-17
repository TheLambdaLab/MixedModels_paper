function  [Y]=powerlaw(b,X)
% se log-log dedomena
Y=b(1)+b(2)*X;
% se mh log-log
% Y=b(1)*X.^b(2);