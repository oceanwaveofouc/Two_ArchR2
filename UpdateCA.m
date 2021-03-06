function [CA ,flagg,flag_offspring]= UpdateCA(CA,New,MaxSize,W)
% Update CA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    n_old=length(CA); 
    n_offspring=length(New); %记录new里面是否存活的flag
    flagg=true(1,n_old);
    flag_offspring=false(1,n_offspring);
    
    CA = [CA,New];
    N  = length(CA);
 
    if N <= MaxSize  %没有超过最大，直接退出
        return;
    end
    
    %% Calculate the fitness of each solution

    
    zmin = min(CA.objs,[],1);
    zmax = max(CA.objs,[],1);
     [Rank,Norm] = R2Ranking(CA.objs,W,zmin,zmax);
     [~,rank]    = sortrows([Rank,Norm]);
     

     

      CA  = CA(rank(1:MaxSize));
      Rank        = Rank(rank(1:MaxSize));
      Norm        = Norm(rank(1:MaxSize));
      
     index=rank(1:MaxSize)  ;  
     index=index-n_old;

     flagg(find(index>0))=false;
     flag_offspring(find(index>0))=true;
      
 
end