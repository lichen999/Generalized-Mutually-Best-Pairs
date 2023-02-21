clc;
clear;
%% Initializing parameters
A = 0:0.25:1; % alpha
L = 0:0.25:1;  % lambda 
beta = 0:0.05:1; % senstivity of distance in preferences
delta = 0:0.05:1;   % sensitivity of distance in priorities
m = 50; % # schools
q= 20;  % capacity per school
n = 1000; % # repetitions

%% Simulations 
a = 1; % choose the number 
alpha = A(a); 

l= 5; % choose the number  
lambda=L(l);
Q=transpose(q*repmat((1:1)',1,m));
MatMBP = cell(n,1);  % sequential mutually best pairs  
MatGMB = cell(n,1); % generalized mutually best pairs
MatTD = cell(n,1); % TTC = DA
MatE = cell(n,1);  % DA is efficient

for rep=1:n
        disp(rep);
        matmbp = zeros (length(delta),length(beta));
        matgmb = zeros (length(delta),length(beta));
        mattd = zeros (length(delta),length(beta));
        mate = zeros (length(delta),length(beta));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % Components for Priority and Preferences 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        % Distance Dist: each row is a student, and each column is a
        % school, each entry is a randomly generated number -- distance --
        % associated to student i in school j, on [-1,0]
        
        Dist = -1 + rand (m*q,m); 

        % G is the common component for priorities, each entry is a
        % student
        G = rand (m*q,1);
        
        % eta is the random coponent for priorities
        eta = rand(m,m*q); 

        % Common component V for preferences 
        V = ones (m*q,m);
        common = rand (m,1);

        for i = 1:m
            V(:,i)=common(i,1)*V(:,i);
        end 

        % Common component X for preferences: different values
        X = ones (m*q,m);
        common2 = rand(m,1);
        for i = 1:m
            X(:,i)=common2(i,1)*X(:,i);
        end
        ai= rand(m*q,1);

     % Random component epsilon 
        e = rand(m*q,m);

        % Combining all partsI: priorities  
        for k = 1: length(beta)

            Priorities = alpha*(beta(k)*Dist'+(1-beta(k))*G')+ (1-alpha)*eta;
            schoolP = zeros(m,q*m);

            for w=1:m
                T=Priorities(w,:);
                [foo,idx] = sort(T(:), 'descend');
                schoolP(w,:)=idx; 
            end
          
        % Combining all parts II: preferences 
             for j =1: length(delta)
     
                Pref = lambda*(delta(j)*Dist + (1-delta(j))*(ai.*X+V)) + (1-lambda)*e;
                studentP = zeros (m*q,m);
    
                for w = 1:q*m   
                    T=(Pref(w,:));
                    [foo,idx] = sort(T(:), 'descend');
                    studentP(w,:)=idx;
                end 

                mbppref = mbp(studentP,schoolP,m,q);
                matmbp(k,j)=mbppref;

                gmbppref = gmbp(studentP,schoolP,m,q);
                matgmb(k,j) = gmbppref;
    
                tt=TTC(studentP,schoolP,Q);
                dd=dastudent(studentP,schoolP,Q);
                mattd(k,j)=all(tt==dd);

                ee=pareto(dd,studentP,schoolP);
                mate(k,j)=ee;
             end
        end 

        MatMBP{rep} = matmbp;
        MatGMB{rep} = matgmb;
        MatTD{rep} = mattd;
        MatE{rep} = mate;

end

save([ 'MatMBP_a' num2str(alpha) '_l' num2str(lambda) '.mat' ],'MatMBP');
save([ 'MatTD_a' num2str(alpha) '_l' num2str(lambda) '.mat' ],'MatTD');
save([ 'MatE_a' num2str(alpha) '_l' num2str(lambda) '.mat' ],'MatE');
save([ 'MatGMB_a' num2str(alpha) '_l' num2str(lambda) '.mat' ],'MatGMB');

%% 3D graphs
n=1000;
perc_td=sum(cat(3,MatTD{1:n}),3)/n; % first concatenate on the 3rd dimention, then sum
perc_e = sum(cat(3,MatE{1:n}),3)/n;
perc_mbp = sum(cat(3,MatMBP{1:n}),3)/n;
perc_gmbp = sum(cat(3,MatGMB{1:n}),3)/n;

figure 
surf (beta,delta,perc_mbp,'FaceColor','b','FaceAlpha',0.5); hold on
surf (beta,delta, perc_gmbp,'FaceColor','r','FaceAlpha',0.5); hold on
surf (beta,delta,perc_e,'FaceColor','g','FaceAlpha',0.5);
legend('Sequential MBP', 'GMBP','DA is efficient','FontSize',12, 'Location', 'best');
xlabel('\beta','FontSize',14); 
ylabel('\delta','FontSize',14);
zticks(0)