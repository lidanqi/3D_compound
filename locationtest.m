function result=locationtest(l1,l2,l3)   
    N1=2^l1;
    N2=2^l2;
    N3=2^l3;
    
    
    % testing
    mat_size = (N1+1)*(N2+1)*(N3+1);
    location_index = zeros(mat_size,3);
    D=zeros(mat_size,4);
    for i=1:N1+1
        for j=1:N2+1
            for k=1:N3+1
              loc = locate_D(i,j,k);
              [ii ,jj,kk]=unlocate(loc);
              if (ii-i)*(jj-j)*(kk-k)~=0
               fprintf('Not the same for %4d, expected: %2d %2d %2d, output: %2d %2d %2d\n',loc,i,j,k,ii,jj,kk);
              end
            end
        end
    end
    
    
    function location = locate_D(i,j,k)
     % mapping matric D(i,j,k) to D((N1+1)*(N2+1)*(N3+1))
     location = (i-1)*(N2+1)*(N3+1)...
              + (j-1)*(N3+1)...
              + (k);
     location_index(location,:)=[i j k];     
    end

    function [ii,jj,kk] = unlocate(loc)
      ii=location_index(loc,1);
      jj=location_index(loc,2);
      kk=location_index(loc,3);
    end
  result=1;
end