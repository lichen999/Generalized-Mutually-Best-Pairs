
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CODE TO RUN TOP TRADING CYCLE (TTC) MECHANISM 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS:
% studentP: (nbstudents x maxchoice) matrix where each row corresponds to a
% student ROL, with schools labels 1 to nbschool (if ROL is shorter than
% maxchoice then zeros should appear - they will be replaced by the student id)
% schoolP: (nbschools x nbstudents) recording students' priorities in each
% school.
% Q is a nbschool x 1 vector of capacities.
% OUTPUT
% Matching is a nbstudents X 1 vector with zeros for unassigned students
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Matching=TTC(studentP,schoolP,Q)


% INITIALISATION / DEFINITION OF VARIABLES / ORGANISATION OF DATA

[nbstudents,maxchoice]=size(studentP);
% In TTC all agents are coded from 1 to nbstudents + nbschools
% Recode schoolids in students' preference so that it starts at nbstudents+1
studentP=studentP+(studentP>0)*nbstudents;
studentP=[studentP [1:nbstudents]']; % will be useful when we loop to update studentP
maxchoice=maxchoice+1;
[nbschools, maxranking]=size(schoolP);
schoolP=[schoolP [nbstudents+1:nbstudents+nbschools]']; % will be useful when we loop to update schoolP
maxranking=maxranking+1;
M=zeros(nbstudents,1);
assigned=zeros(1, nbstudents+nbschools);

% CHECK INTEGRITY OF INPUTS
for i=1:nbstudents
if any(diff(sort(studentP(i,:)))==0)
    print('Some student(s) rank the same school twice');
end

end
for i=1:nbschools
if any(diff(sort(schoolP(i,:)))==0)
    print('Some school(s) rank the same student twice');
end
end

%Complete student preferences and school priorities so that students point
%to themselve when their ROL is exhausted and same thing for schools
studentP=studentP+(studentP==0).*repmat([1:nbstudents]',1,maxchoice);
schoolP=schoolP+(schoolP==0).*repmat([1:nbschools]',1,maxranking);

% make sure that schools with no capacity point to themselves
schoolP=(1-repmat((Q==0), 1, maxranking)).*schoolP + repmat((Q==0).*[1:nbschools]',1,maxranking);

% make sure students point to a school with positive capacity or update
% preferences accordingly otherwise
for i=1:nbstudents
    if Q(studentP(i:1)-nbstudents)==0
        % remember nbstudents was added to all school labels
        studentP(i,:)=[studentP(i,2:maxchoice) i];
        % only accounting for first choice but ok given how the loop is coded
    end
end


% TTC algo:

% Pointing vector 1 x (nbstudents+nbschools) record first choice of students
% and first priority student of schools.
Pointing=horzcat(studentP(:,1)', schoolP(:,1)');

S=0; % initialisation, S=nbstudents+nbschools when TTC ends because we've
% reached a stage where all students and schools point to themselves (! S =
% nbstudents+nbschools is a necessary condition not a sufficient condition
% for convergence. Will be sufficient if preferences and priorities are
% fully updated at each stage.

while S < nbstudents+nbschools

% create the sparse matrix corresponding to the (nbstudents + nbschools) x
% (nbstudents + nbschools) graph matrix where students and schools point
% to the preferred partner.
% G = sparse(i,j,v, m, n) generates a sparse m x n matrix G from the
% triplets i, j, and v such that G(i(k),j(k)) = v(k). The max(i)-by-max(j)
% output matrix has space allotted for length(v) nonzero elements. sparse
% adds together elements in v that have duplicate subscripts in i and j.

G=sparse([1:(nbstudents+nbschools)],Pointing,true,(nbstudents+nbschools),(nbstudents+nbschools));

[S,C] = graphconncomp(G);

% S = the number of connected components in the graph
% S = nbstudents+nbschools when all students assigned and all schools have
% exhausted their list, but it can be equal to that number earlier
% C = 1 x (nbstudents+nbschools) vectors that identifies the cycle to which
% the node belongs. If the node does not belong to any cycle, then it forms
% a component by itself.
% See http://nl.mathworks.com/help/bioinfo/ref/graphconncomp.html

% For each cycle, we identify the students and schools that participate and
% update the matching matrix M accordingly and update school capacities

for i=1:S
    if sum(C==i)>1 %Only consider real cycles (with length more than 1)
       assigned=(C==i); % records the students and schools in the cycle
       assignedstudent=find(assigned(1:nbstudents));
       M=M+(assigned(1:nbstudents).*Pointing(1:nbstudents))'; % updates matching of students
       Q=Q-assigned(nbstudents+1:nbstudents+nbschools)'; % updates capacity
    end
end

% Update pointing vector:
% Schools:
% we check if remaining first priority student is not assigned and
% remove him otherwise (we loop to make sure the school does not point to
% an already assigned student

for i=1:nbschools
    if schoolP(i,1)<nbstudents+1 % that school points to one student
        for k=1:find(schoolP(i,:)>nbstudents,1)-1 % k is well defined given that we've added the school pointing to itself at the end
            if M(schoolP(i,1))>0  % that student is assigned
                schoolP(i,:)=[schoolP(i,2:maxranking) i+nbstudents];
            end
        end
    end
end
% schools with zero capacity point to themselves
schoolP=repmat((1-(Q==0)),1, maxranking).*schoolP+repmat((Q==0).*[nbstudents+1:nbstudents+nbschools]',1,maxranking);

Pointing(nbstudents+1:nbstudents+nbschools)=schoolP(:,1);

% students:
% assigned students point to themselves
Pointing(1:nbstudents)=(1-(M==0))'.*[1:nbstudents]+(M==0)'.*Pointing(1:nbstudents);

% student preferences and pointing vector are updated by removing the schools with no capacity
for i=1:nbstudents
    if M(i)==0
        for k=1:find(studentP(i,:)==i,1)-1
            if Q(studentP(i,1)-nbstudents)==0
                studentP(i,:)=[studentP(i,2:maxchoice) i];
                Pointing(i)=studentP(i,1);
            end
        end
    end
end

end

% if we reach this point then, all students are pointing to themselves
% (i.e. they're assigned or have no more acceptable schools) and schools
% are pointing to themselves

%Final Matching with original school codes

Matching=(M>0).*(M-nbstudents);

end
