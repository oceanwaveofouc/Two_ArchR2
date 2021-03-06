function DA = UpdateDA(CA,DA,New,flag,MaxSize,p)
% Update DA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Find the non-dominated solutions
    DA_oldsize=length(DA);
    DA = [DA,New];
    Next   = 1:length(DA);
    Choose = Association(CA.objs ,DA.objs,MaxSize,DA_oldsize,flag);
    Next   = Next(Choose);
    % Population for next generation
    DA = DA(Next);
    
end


function Choose = Association(PopObj1,PopObj2,N,DA_oldsize,flag)
% Association operation in the algorithm

    [N1,~] = size(PopObj1);
    
    [N2,M] = size(PopObj2);
    
%     Choosess = false(1,N2);
%     for i=1:N2
%         for j=1:N1
%             if sum( PopObj2(i,:)==PopObj1(j,:)   )==M
%             Choosess(i)=true;
%             end
%         end
%     end
%     
%     aa=sum(Choosess)
    PopObj = [PopObj2];
    N1=0;
    %% Normalization
    Zmin   = min(PopObj,[],1);
    Zmax   = max(PopObj,[],1);
    PopObj = (PopObj-repmat(Zmin,size(PopObj,1),1))./repmat(Zmax-Zmin,size(PopObj,1),1);
    
    %% Calculate the fitness value of each solution
    fit = sum(PopObj,2);
    
    %% Angle between each two solutions
    angle = acos(1-pdist2(PopObj,PopObj,'cosine'));
    
    %% Niching
%     Choose = false(1,N2);
    
     Choose = [false(1,DA_oldsize) flag]  ;
    
%      if ~any(Choose)
        % Select the extreme solutions first
        [~,extreme]        = min(pdist2(PopObj2,eye(M),'cosine'),[],1);
        Choose(N1+extreme) = true;
        % Select the first M best converged solutions
        [~,rank] = sort(fit(N1+1:end));
        Choose(N1+rank(1:min(M,length(rank)))) = true;
%       end
    
    
    while sum(Choose) < N
        % Maximum vector angle first
        Select  = find(Choose);
        Remain  = find(~Choose);
        [~,rho] = max(min(angle(Remain,Select),[],2));
        Choose(Remain(rho)) = true;
        % Worse elimination
        if ~all(Choose)
            Select      = [Select,Remain(rho)];
            Remain(rho) = [];
            [~,mu]      = min(min(angle(Remain,Select),[],2));
            [theta,r]   = min(angle(Remain(mu),Select));
            if theta < pi/2/(N+1) && fit(Select(r)) > fit(Remain(mu))
                Choose(Select(r))  = false;
                Choose(Remain(mu)) = true;
            end
        end
    end
end