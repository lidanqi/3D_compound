% sparse grid combination step : American option
function [estimation, details] = MainFuncA(requiredlevel,)
% set benchmark for comparison:  POSR
S = [0.8, 0.9, 1, 1.1, 1.2];
benchmark = [19.9987, 10.9820, 5.4899, 2.6295, 1.2388];
benchmark_idx = floor((S-0.7)/0.1);
benchmark = benchmark(benchmark_idx);
% start
timer = cputime;

sums=zeros(3,5);
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
        fprintf('levels: %d %d %d \n', i,j,k);
        [~,temp] = mainA(S,200,'level',[i j k]);
    %  temp = 0;
        details=[details;i j k temp'];
        if isnan(temp)
            temp=0;
        end
        answers = [answers; temp'];
      %  fprintf('%5d %5d %5d %10f\n',i,j,k,temp);
        list(index,:)=[i j k points];
        index=index+1;
    end
end
sums(level_s)=sum(answers);


end

sums';
estimation = sums(1) - 2*sums(2) + sums(3); 
estimation = estimation*100;

fprintf('The estimated Option price: '); fprintf('%6g ',estimation); 
fprintf('\n              Benchmark   : '); fprintf('%6g ',benchmark);
fprintf('\n          Total time spent: %4d s \n',cputime-timer);

