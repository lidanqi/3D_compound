% sparse grid combination step
% Compound Option
function [estimation, details,timespent] = MainFuncC(requiredlevel)
% set benchmark for comparison: MC + POSR
S = [0.8, 0.9, 1, 1.1, 1.2];
benchmark = [16.9512, 9.6856, 5.0227, 2.4618, 1.1774];
benchmark_C = [0.1032, 0.6079, 1.5228, 2.4339, 3.0889];

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
            [~,~,est,~] = mainC(S,128,'level',[i j k]);
            save(filename,'est');
        else
            loadedEst = load(filename);
            est=loadedEst.est;
        end
        temp = est(:,1);
        tempC = est(:,2);
    %  temp = 0;
        details=[details;i j k temp' tempC'];
        if isnan(temp)
            temp=0;
        end
        answers = [answers; temp' tempC'];
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
fprintf('\n               Benchmark   : '); fprintf('%6g ',benchmark);
fprintf('\n---------------------------------------------------------------------');
fprintf('\nEstid Mother Option price: '); fprintf('%6g ',estimation_AC); 
fprintf('\n               Benchmark   : '); fprintf('%6g ',benchmark_C);
timespent = cputime-timer;
fprintf('\nTotal time spent: %4d s \n',timespent);

