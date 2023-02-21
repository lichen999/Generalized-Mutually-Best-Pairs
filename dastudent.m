%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         CODE TO RUN  STUDENT PROPOSING-DA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUTS
% Q is a nbschool x 1 vector of capacities
% INPUTS:
% studentP: (nbstudents x maxchoice) matrix where each row corresponds to a
% student ROL, with schools labels 1 to nbschool (if ROL is shorter than
% maxchoice then zeros should appear - they will be replaced by the student id)
% schoolP: (nbschools x nbstudents) recording students' priorities in each
% school.
% Q is a nbschool x 1 vector of capacities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function match=dastudent2(studentP,schoolP,Q)


% INITIALISATION / DEFINITION OF VARIABLES / ORGANISATION OF DATA

[nbstudents,~]=size(studentP);
[nbschools, maxranking]=size(schoolP);

% notProposed records elements of studentP that corresponds to schools to
% which students have not proposed yet.
% Initialisation: all values = 1 except for 0 entries of studentP
notProposed=true(size(studentP));%.*(studentP>0); this is not logical operator

% keep track of which student is unassigned
unassigned=true(nbstudents,1);
% keep track of tentative assignment
propTent=false(nbschools, maxranking);
matcht=zeros(nbstudents,1);

iter=0;
%============ loop begins ===============
while any(any(notProposed(find(unassigned~=0),:),2)) && ~all(matcht)
    %any(unassign) && any(any(notProposed,2))
     iter=iter+1;

% students proposing...
proposer=find(any(unassigned,2));
canprop=any(notProposed(proposer,:),2);
proposers=proposer(canprop);
favSchool=zeros(nbstudents,1);

for i=proposers(1:end)'
    schooltopropose=studentP(i,notProposed(i,:));
    favSchool(i)=schooltopropose(1);
    [~,order]=find(studentP(i,:)==schooltopropose(1));
    notProposed(i,order)=false;
end

% school deciding...
school=nonzeros(unique(favSchool(proposers))); % records which schools are listed first
for j=school(1:end)'
    propReceived=(propTent(j,:)==true);
    comp=find(favSchool==j);
    propReceived(comp)=true;
    allcomp=find(propReceived==true);
    matcht(propTent(j,:))=0;
    propTent(j,:)=false;
    if sum(propReceived)<=Q(j)
        best=allcomp;
    else
        candidates=intersect(schoolP(j,:),allcomp,'stable');
        best=candidates(1:Q(j));
    end
    propTent(j,best)=true;
    matcht(best)=j;
end

unassigned=(~all(matcht,2));
%====== loop finshed ======================

end
match=matcht;
end
