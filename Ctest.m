function result = CtestS()
%%testing compound option
S=[0.8,0.9,1,1.1,1.2];
% level [ 2 2 2 ];
level1est=[];
for i=1:5
    [est,~] = MainFuncC(1,S(i));
    level1est=[level1est est];
end
level1est = level1est*100;
save('level1sparse.mat','level2est');

% level [ 3 3 3 ];
level2est=[];
for i=1:5
    [est,~] = MainFuncC(2,S(i));
    level2est=[level2est est];
end
level2est = level2est*100;
save('level2sparse.mat','level3est');

true_est = [0.1032 0.6079 1.5228 2.4339 3.0889];
result=[level1est; level2est; true_est];
save('compoundTestingSparse.mat','resultSparse');

