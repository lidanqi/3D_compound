    function final_est =interpolation()
        i=floor(N1/2)+1;
        j=floor(N2/2)+1;
        k=floor(N3/2)+1;
        
     if (l1==0 && l2~=0 && l3~=0)
         final_est=mean([D_new(locate_D(1,j,k)) D_new(locate_D(2,j,k))]);
         return;
     end
     if (l1~=0 && l2==0 && l3~=0)
         final_est=mean([D_new(locate_D(i,1,k)) D_new(locate_D(i,2,k))]);
         return;
     end
     if (l1~=0 && l2~=0 && l3==0)
         final_est=mean([D_new(locate_D(i,j,1)) D_new(locate_D(i,j,2))]);
         return;
     end
     
     if l1~=0
         final_est = interpolation2(vmin,vmax,rmin,rmax,D_new(locate_D(i,1,1)),D_new(locate_D(i,1,2)),D_new(locate_D(i,2,1)),D_new(locate_D(i,2,2)));
         return;
     end
     if l2~=0
         final_est = interpolation2(Smin,Smax,rmin,rmax,D_new(locate_D(1,j,1)),D_new(locate_D(1,j,2)),D_new(locate_D(2,j,1)),D_new(locate_D(2,j,2)));
         return;
     end
     if l3~=0
         final_est = interpolation2(Smin,Smax,vmin,vmax,D_new(locate_D(1,1,k)),D_new(locate_D(1,2,k)),D_new(locate_D(2,1,k)),D_new(locate_D(2,2,k)));
         return;
     end
     
    end