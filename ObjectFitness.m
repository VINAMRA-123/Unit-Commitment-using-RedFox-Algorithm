function [Fobject]=ObjectFitness(x,G,T,x0,MDT,Pmax,Pmin,PD,a,b,c,SUH,SUC,Tcold,TOFF) 

 options=optimset('LargeScale','off','MaxPCGIter',50,'PrecondBandWidth',inf); 
 for t = 1:T 
    H = [];  
    C = [];  
    Aeq = []; 
    VLB = []; 
    VUB = []; 
    i = 0; 
    for g = 1:G 
        if x(g+G*(t-1)) == 1 
            i = i+1; 
            H(i,i) = 2*c(g); 
            C(i) = b(g); 
            Aeq(i) = 1; 
            VLB(i) = Pmin(g); 
            VUB(i) = Pmax(g); 
            p0 = 0.5*(VLB+VUB); 
       end 
    end 
 
    beq = PD(t);
    H;
    C';
    Aeq;
    beq;
    VLB;
    VUB;
    p0;
    [Pt,F] = quadprog(H,C',[],[],Aeq,beq,VLB,VUB,p0,options);
   
    T
    F
    a
    x(G*(t-1)+1:G*t)'
    FC(t) = F+a*x(G*(t-1)+1:G*t)';
    FC(t)

%% Start-up cost 

    if t == 1 
        for g = 1:G 
            if x(g+G*(t-1)) == 1&&x0(g) == 0 
                 
                    SU = SUH(g); 
                
            end 
        end 
    else 
        for g = 1:G 
            if x(g+G*(t-1)) == 1&&x(g+G*(t-2))==0 
                 if TOFF(g,t-1) >= MDT(g)&& TOFF(g,t-1) <= MDT(g)+Tcold(g) 
                       SU = SUH(g); 
                 elseif TOFF(g,t-1) > MDT(g)+Tcold(g) 
                         SU = SUC(g); 
                 end                             
            end    
        end 
    end 
end  

  Fobject = SU + sum(FC); % Total cost