clear all; 
clc; 


Pmax        = [455,455,130,130,162,80,85,55,55,55];                                             % Maximal power output (MW)
Pmin        = [150,150,20,20,25,20,25,10,10,10];                                                % Minimal power output (MW)
a           = [1000,970,700,680,450,370,480,660,665,670];                                       % $/h
b           = [16.19,17.26,16.6,16.5,19.7,22.26,27.74,25.92,27.27,27.79];                       % $/MWh
c           = [0.00048,0.00031,0.002,0.00211,0.00398,0.00712,0.00079,0.00413,0.00222,0.00173];  % $/MWh^2
SUH         = [4500,5000,550,560,900,170,260,30,30,30];                                         % Start-up HOT cost ($)
SUC         = [9000,10000,1100,1120,1800,340,520,60,60,60];                                     % Start-up COLD cost ($)
Tcold       = [5,5,4,4,4,2,2,0.5,0.5,0.5];                                                      % Cooling time unit
InitialTON  = [8,8,0,0,0,0,0,0,0,0];                                                            % Initial ON unit
InitialTOFF = [0,0,-5,-5,-6,-3,-3,-1,-1,-1];                                                    % Initial OFF unit
MDT         = [8,8,5,5,6,3,3,1,1,1];                                                            % Min. Down Time (MIN OFF)
MUT         = [8,8,5,5,6,3,3,1,1,1];                                                            % Min. UP Time (MIN ON)
x0          = [1,1,0,0,0,0,0,0,0,0];                                                            % Initial condition


DEMAND = [700; 750; 850; 950; 1000; 1100; 1150; 1200; 1300; 1400; 1450; 1500; 1400; 1300; 1200; 1050; 1000; 1100; 1200; 1400; 1300; 1100; 900; 800;];

Reserves    = 0.1*DEMAND;
[Y,h]       = sort(Pmax,'descend'); 

SA          = 10; 
G           = 10; 
T           = 24; 

MaxIt       = 50;                     % Maximal iteration
c1 = 0.18;                            % range of c1 is [0, 0.18]
c2 = 0.82;                            % range of c2 is [0.19, 1]
[Y,h]       = sort(Pmax,'descend'); 

redfox       = zeros(SA,G*T); 
ab          = 1/sqrt(2)*ones(2,G*T); 


for n = 1:SA 
    redfoxl = initial(G,T); 
    [redfox(n,:),redfoxTOFF(:,:,n)] = Verify(redfoxl,G,T,x0,MUT,MDT,InitialTON,InitialTOFF,Pmax,DEMAND,Reserves,h); 
    [fitness(n)] = ObjectFitness(redfox(n,:),G,T,x0,MDT,Pmax,Pmin,DEMAND,a,b,c,SUH,SUC,Tcold,redfoxTOFF(:,:,n)); 
end 
Gmbest = redfox(1,:); 
Gfitness = fitness(1); 

for n = 1:SA 
    Pmbest(n,:) = redfox(n,:); 
    Pfitness(n) = fitness(n); 
    if Gfitness > Pfitness(n) 
        Gmbest = Pmbest(n,:); 
        Gfitness = Pfitness(n); 
    end 
end 
 
for n = 1:SA 
    ab1(:,:,n) = ab; 
end

Distance_Fox_Rat = zeros(SA, G*T);
Gfitness_values = zeros(1, MaxIt);
MinT = inf;
for k = 1:MaxIt 

    for n = 1:SA 
       [redfox(n,:),ab1(:,:,n)] = Update(redfox(n,:),G,T,Pmbest(n,:),Gmbest,fitness(n),Pfitness(n),Gfitness,c2,c1,k,MaxIt,ab1(:,:,n),MinT,n); 
       [redfox(n,:),redfoxTOFF(:,:,n)] = Verify(redfox(n,:),G,T,x0,MUT,MDT,InitialTON,InitialTOFF,Pmax,DEMAND,Reserves,h); 
       [fitness(n)] = ObjectFitness(redfox(n,:),G,T,x0,MDT,Pmax,Pmin,DEMAND,a,b,c,SUH,SUC,Tcold,redfoxTOFF(:,:,n)); 
  
        if Pfitness(n) >= fitness(n) 
            Pmbest(n,:) = redfox(n,:); 
            Pfitness(n) = fitness(n); 
        end 
        if Gfitness > Pfitness(n) 
        Gmbest = Pmbest(n,:); 
        Gfitness = Pfitness(n); 
        end         
    end
    Gfitness_values(k) = Gfitness;
end 



options=optimset('LargeScale','off','MaxPCGIter',50,'PrecondBandWidth',inf); 
for t = 1:T 
    H = 2*diag(c.* Gmbest(G*(t-1)+1:G*t)); 
    C = b.* Gmbest(G*(t-1)+1:G*t); 
    Aeq = ones(1,10).* Gmbest(G*(t-1)+1:G*t); 
    VLB = Pmin; 
    VUB = Pmax; 
    beq = DEMAND(t); 
    p0 = 0.5*(VLB+VUB); 
    P = quadprog(H,C',[],[],Aeq,beq,VLB,VUB,p0,options); 
    Gptbest(:,t) = P'.*Gmbest(G*(t-1)+1:G*t);
end 
 
Gptbest'
disp('-----------------------------------------------------------------------');
resu = Gfitness;
fprintf('Total Cost is $ %.6f\n', resu);
disp('-----------------------------------------------------------------------');
figure;
plot(1:MaxIt, Gfitness_values, 'LineWidth', 1.5);
xlabel('Generation');
ylabel('Total Generation Cost $');
title('Graph of Total Cost vs Generation');