%[estimation, details] = MainFunc(requiredlevel,S)
function [estimation, details] = MainFunc(requiredlevel,S)
tic
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
        [~,temp] = main(S,500,'level',[i j k]);
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

% for printidx=1:size(list,1)
%  fprintf('\n%5d %5d %5d %10d',list(printidx,1),list(printidx,2),list(printidx,3),list(printidx,4));
% end
end

sums';
estimation = sums(1) - 2*sums(2) + sums(3); 
toc
