% this function is used for multilinear interpolation
function estimation=interpolation(Svec,vvec,rvec,resultVect,S,v,r)
sizeS = size(Svec,1);
sizev = size(vvec,1);
sizer = size(rvec,1); 
[Y,X,Z] =meshgrid(vvec,Svec,rvec);

% reshape result vector into 3-D
[res] = mymesh(resultVect);

% estimation = interp3( X,Y,Z,resultVect ,S,v,r);
 F = griddedInterpolant(X,Y,Z,res,'cubic');
% F = griddedInterpolant(X,Y,Z,res,'linear');
estimation = F(S,v,r);
% estimation = interp3(x,y,z, res, S,v,r); 


function [resultVect] = mymesh(V)
    N1 = sizeS-1;
    N2 = sizev-1;
    N3 = sizer-1;
    
    resultVect=zeros(sizeS,sizev,sizer);
%     x=resultVec;y=x; z=x; 
        for ii=1:N1+1
            for jj=1:N2+1
                for kk=1:N3+1
%                     x(ii,jj,kk)=s(ii);
%                     y(ii,jj,kk)=v(jj);
%                     z(ii,jj,kk)=r(kk);
                    resultVect(ii,jj,kk) = V(locate_D(ii,jj,kk));
                end
            end
        end
    
%% mapping
%   location mapping: 3D to 1D
    function location = locate_D(i,j,k)
     % mapping matric D(i,j,k) to D((N1+1)*(N2+1)*(N3+1))
     location = (i-1)*(N2+1)*(N3+1)...
              + (j-1)*(N3+1)...
              + (k);
     location_index(location,:)=[i j k];     
    end

%   location mapping: 1D to 3D
    function [ii,jj,kk] = unlocate(loc)
      ii=location_index(loc,1);
      jj=location_index(loc,2);
      kk=location_index(loc,3);
    end
end

end

