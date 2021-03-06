function Two_ArchR2(Global)
% <algorithm> <T>
% Two-archive algorithm 2
% CAsize --- --- Convergence archive size
% p      --- --- The parameter of fractional distance

%------------------------------- Reference --------------------------------
% H. Wang, L. Jiao, and X. Yao, Two_Arch2: An improved two-archive
% algorithm for many-objective optimization, IEEE Transactions on
% Evolutionary Computation, 2015, 19(4): 524-541.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [CAsize,p] = Global.ParameterSet( (Global.N) ,1/Global.M);
    
    %% Generate random population
    [W,Global.N] = UniformPoint(Global.N,Global.M);

    Population = Global.Initialization();
   [ CA,~]= UpdateCA([],Population,CAsize);
%      Population = Global.Initialization();
%     DA = UpdateDA(CA,Population,[],Global.N,p);
     DA=CA;
     Zmin = min(Population(all(Population.cons<=0,2)).objs,[],1);
%     DA=[];
    %% Optimization
    t=0;
    while Global.NotTermination([CA DA])
        t=t+1;
        [ParentC,ParentM] = MatingSelection(CA,DA,Global.N); % ParentC CA+DA  , ParentM,DA 交配池默认都为N大小
         Offspring         = [GA(ParentC) ,GA(ParentM)];  
         Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
        
         [CA,Choosess,flag_offspring] = UpdateCA(CA,Offspring,CAsize,W);
         aa=sum(Choosess);
         flag_offspring=false;
         chooseflagoff=[false(1,aa)    flag_offspring ];
         DA = UpdateDA(CA,DA,[CA(Choosess) Offspring],chooseflagoff,Global.N,p);

    end
end