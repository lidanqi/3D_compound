% sparse grid combination step
% American Compound Option
function [estimation, details] = MainFuncAC(requiredlevel,S)
% set benchmark for comparison: MC + POSR
benchmark = [0.1072, 0.6119, 1.5618, 2.5233, 3.1928];
benchmark_idx = floor((S-0.7)/0.1);
benchmark = benchmark(benchmark_idx);
% start 
timer = cputime;
sums=zeros(3,1);
dimension=3;
details=[];
for level=requiredlevel:requiredlevel+2

level_s = level - requiredlevel + 1;
list=zeros(nchoosek(level+dimension-1,dimension-1),4);

answers=[];

index=1;

for i=0:level
    for j=0:(level-i)
        k=level-i-j;
        % number of points on this grid of combination (i,j,k)
        points = (2^i+1)*(2^j+1)*(2^k+1);
        [~,~,est,~] = mainAC(S,128,'level',[i j k]);
        temp = est(2);
    %  temp = 0;
        details=[details;i j k temp];
        if isnan(temp)
            temp=0;
        end
        answers = [answers; temp];
        fprintf('%5d %5d %5d %10f\n',i,j,k,temp);
        list(index,:)=[i j k points];
        index=index+1;
    end
end
sums(level_s)=sum(answers);

end

sums';
estimation = sums(1) - 2*sums(2) + sums(3); 
fprintf('\n==============================================================\n');
fprintf(['\nThe estimated American Mother Option price: %6f'  ...
         '\n                              Benchmark   : %6f ' ...
         '\n                          Total time spent: %4f s\n'],  ...
        estimation*100,benchmark,cputime-timer);

