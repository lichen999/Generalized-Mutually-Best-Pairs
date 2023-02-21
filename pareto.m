% ================================================
% CODE TO CHECK IF DA IS EFFICIENT
% May, 2022
% ==================================================
% studentP= [1 2 3;1 3 2; 2 1 3];
% schoolP = [3 2 1 ;1 3 2; 2 3 1];
% Q=[1;1;1];

function E =pareto(DA,studentP,schoolP)

% parameters
nstudents = length(studentP);
nschools = length(schoolP);
studentid = 1:nstudents;

% students who haven't received their first choices would like to exchange
POS = zeros(nstudents,1); %position assigned
for i=1:nstudents
    if DA(i,1) ~= 0
        POS(i,1) = studentid(studentP(i,:) ==DA(i,1));
    end
end

totrade = POS~=1; % totrade is an index vector =true if a student is not assigned to his first
% choice
tradeid = studentid(totrade); % ids of these students

% Construct digraph for trading
% Part I: Students who are not assigned to first choices, and their more preferred choices
id1 = cell(1,length(tradeid));  % ids of students who would like to trade
better = cell(1,length(tradeid)); % schools that the students would like to trade

for k = 1:length(better)
    better{k} = studentP(tradeid(k),1:(POS(tradeid(k))-1));
    id1{k}=repmat(tradeid(k),1,length(better{k}));
end

better = cell2mat(better); % transform better to matrix,
better = better + size(studentP,1);  % start to count school id as of
% nstudents+1

id1 = cell2mat(id1);  % transform id1 to matrix, and make x copies equal to the number of more preferred choices

% Part II: Students who are currently assigned to a shool

current = cell(1,nschools); %for each school, the students assigned to it
id2 = cell(1,nschools);    % the ids of these schools

for k = 1:length(schoolP)
     current{k} = studentid(DA(:)==k);
     id2{k} =repmat(k,1,length(current{k})) ;
end

current = cell2mat(current);
id2 = cell2mat(id2);
id2 = id2 + size(studentP,1); % relabel schools

% combine the nodes
node1 =[id1 id2];  % node1 contains the ids for students who would like to trade, and all the schools
node2 = [better current]; % node 2 contains the schools that students would like to exchange (part I), and the students that are currently held by schools (part II)

% diagraph where the 1st element of node1 is connected to 1st element
% of node2, 2nd element of node2 is connected to 2nd element of node2, and
% so on
G = digraph(node1,node2);

[~,bsize]=conncomp(G);

% if one component involves more than and equal to 4 nodes, then there is a
% cycle: student points to a more preferred school, which is held by a
% student (let the school points to the student), who points to another
% school, which is held by a student (let the school points to the
% student), so for a cycle, we need minimum 4 connections.
if any(bsize>=4)
E = false;
else
    E = true;

end
end
