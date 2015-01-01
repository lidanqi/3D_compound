function [ list ] = Points( requiredlevel )
% input a level value, list out all possible combinations and number of
% points 
dimension=3;
% three directions i,j,k satisfying i+j+k=level and >=0, <=level

count_point=[0 0 0];
count_grid=count_point;

for level=requiredlevel:requiredlevel+2
    
list=zeros(nchoosek(level+dimension-1,dimension-1),4);

index=1;
for i=0:level
    for j=0:(level-i)
        k=level-i-j;
        % number of points on this grid of combination (i,j,k)
        points = (2^i+1)*(2^j+1)*(2^k+1);
        list(index,:)=[i j k points];
        index=index+1;
    end
end

% for printidx=1:size(list,1)
%  fprintf('\n%5d %5d %5d %10d',list(printidx,1),list(printidx,2),list(printidx,3),list(printidx,4));
% end
count_point(level-requiredlevel+1)=sum(list(:,4));
count_grid(level-requiredlevel+1)=index-1;
end

fullpoints=(2^requiredlevel+1)^dimension;
sparsepoints_max = max(count_point);
sparsepoints=sum(count_point);
sparsegrids=sum(count_grid);


diff_percent=(fullpoints-sparsepoints)/fullpoints;
fprintf('\nAt level %d, \n  A Full grid has %d points',requiredlevel,fullpoints);
fprintf('\n  Sparse grid has in total %d points, max: %d points',sparsepoints,sparsepoints_max);
fprintf(', by %d + %d + %d ',count_point(1),count_point(2),count_point(3));
fprintf('\n  Sparse grid has in total %d grids',sparsegrids);
fprintf(', by %d + %d + %d ',count_grid(1),count_grid(2),count_grid(3));
fprintf('\nReduction of %8f percecnt realised\n\n',diff_percent*100);

end

