% sparse grid combination step
% American Compound Option
function [estimation, details, timespent] = MainFuncAC(requiredlevel)
% set benchmark for comparison: MC + POSR
S = [0.8, 0.9, 1, 1.1, 1.2];
benchmark_A = [19.9987, 10.9820, 5.4899, 2.6295, 1.2388];
benchmark_AC = [0.1072, 0.6119, 1.5618, 2.5233, 3.1928];

% start 
timer = cputime;
sums=zeros(3,10);
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
        % generate a file containing estimation 
        filename=[num2str(i) num2str(j) num2str(k) '.mat'];
        IScomputed = exist(filename);
        if (IScomputed==0)
            [~,~,est,~] = mainAC(S,256,'level',[i j k]);
            save(filename,'est');
        else
            loadedEst = load(filename);
            est=loadedEst.est;
        end
        tempA = est(:,1);
        tempAC = est(:,2);
    %  temp = 0;
        details=[details;i j k tempA' tempAC'];
        if isnan(tempA)
            tempA=0;
        end
        answers = [answers; tempA' tempAC'];
        list(index,:)=[i j k points];
        index=index+1;
    end
end
sums(level_s,:)=sum(answers);

end

sums';
estimation = (sums(1,:) - 2*sums(2,:) + sums(3,:))*100; 
estimation_A = estimation(1:5);
estimation_AC = estimation(6:10);

fprintf('=====================================================================');
fprintf('\nEstid Daughter Option price: '); fprintf('%6g ',estimation_A); 
fprintf('\n               Benchmark   : '); fprintf('%6g ',benchmark_A);
fprintf('\n---------------------------------------------------------------------');
fprintf('\nEstid Daughter Option price: '); fprintf('%6g ',estimation_AC); 
fprintf('\n               Benchmark   : '); fprintf('%6g ',benchmark_AC);
timespent = cputime-timer;
fprintf('\nTotal time spent: %4d s \n',timespent);

