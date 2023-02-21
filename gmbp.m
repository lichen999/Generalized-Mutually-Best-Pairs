%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                        %
%   	CODE TO CHECK Mutally Best Pairs on truncated market                                       %
%                                                                         %
%  			Aug 2022                                                   %                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% student preferences - studentP
% school priorities - schoolP
% nbschools:m
% nbseats:q
% Output: a logical number = true if satisfying HTTqt, = false otherwise 

function bmp = gmbp(studentP,schoolP,m,q) 

% parameteres 
nbstudents = length(studentP);
Q=transpose(q*repmat((1:1)',1,m));
idstudent = (1:nbstudents)'; % student and school ID to keep track

% I. remove unattractive schools for students 
NTstudent = true(nbstudents,m); % Not truncated matrix for students, true = not removed
NTschool = true (m,nbstudents); % Not truncated matrix for schools (update according to student truncation)

% SEARCHschool = false(m,nbstudents); % Index matrix that stores currently searched top candidates 
% TopschoolNew = false(m,nbstudents);

schoolNewTop = true;
schoolsTopStudent = zeros(m,q);
schoolsTopUpdate = zeros(m,q);

while schoolNewTop == true % problem here 

    for ss = 1:m
    
        schoolPNT = schoolP(ss,NTschool(ss,:)); % for this school, students that haven't been removed
        topstudents = schoolPNT(1:q);  % top q students
        schoolsTopStudent(ss,:) = topstudents; 
        
        for ii=topstudents(1:end)
    
    %         ps=schoolP(ss,:)==ii;
    %         SEARCHschool(ss,ps) = true;
            
            studentPNT = studentP(ii,NTstudent(ii,:)); % for this student, schools that haven't been removed
            
            pos = find(studentPNT(:)==ss); % position of the school where the student has priority
    
            if pos<length(studentPNT) % if not the last choice (there is an unattractive school)
    
                NTstudent(ii,pos+1:end) = false;
               
                schooltodelete = studentPNT(pos+1:end);
    
                for jj = schooltodelete(1:end)
    
                    poss = schoolP(jj,:)==ii;
    
                    NTschool(jj,poss) = false;
                end 
    
            end
        end 
    end

    for ss=1:m
        schoolPNT = schoolP(ss,NTschool(ss,:));
        topstudents = schoolPNT(1:q);
        schoolsTopUpdate(ss,:) = topstudents;
    end

    if schoolsTopStudent==schoolsTopUpdate
        schoolNewTop = false;
    end 

end

% % Update truncated preferences and priorities 
% studentPt = zeros(size(studentP));
% schoolPt = zeros(size(schoolP));
% 
% for ii=1:nbstudents
%     for ss = 1:m
%         if NTstudent(ii,ss)==true
%         studentPt(ii,ss) = studentP(ii,ss);
%         end
%     end 
% end 
% 
% for ss=1:m
%     for ii = 1:nbstudents
%         if NTschool(ss,ii) ==true 
%             schoolPt (ss,ii) = schoolP(ss,ii);
%         end 
%     end 
% end 


% II. Checking Mutally Best Pairs 
% null hypothesis: htt = true, and whenever we can't find a top top pair, htt is changed to false
bmp=true;

% student and school are logical matrix, true = haven't been deleted (after
% finding top top pair)
student = NTstudent;
school  = NTschool; 

% candidate is a vector indicating whether the student is to be considered
% for top top pair
candidate = true(nbstudents,1);

while any(bmp) && any(candidate)

    % tt is logical vector to indicate whether a student can find a school
    % that is each other's top choice
    
    tt = false(size(studentP,1),1); 
    students_cand = idstudent(candidate);
    
    for ii = students_cand(1:end)'
        
        % consider only preferences that haven't be deleted due to TT 
        pref_cand = studentP(ii,:);
        pref_cand = pref_cand(student(ii,:));
        
        % consider of priorities that haven't be deleted due to TT
        priority_cand = schoolP(pref_cand(1),:);
        
        %Delete students with a top top match before 
        priority_cand = priority_cand(school(pref_cand(1),:)); 
       
        % top q available capacity 
        q_available = Q(pref_cand(1));
        tt (ii) = ismember(ii,priority_cand(1:q_available));
    end 
    
    if all (tt ==0) % if all elements of tt == 0, no top top to be found
       
        bmp=false;
        
    else
        
        %%We first implement all the top-top pairs that we found
        while any(tt)
        
        %if there are multiple top top pairs, just pick one
        topstudent = idstudent(tt); 
        topstudent = topstudent(1);       
            
        %%HERE BE CAREFUL! Previously was later, but I need to compute
        %%this before candidate is changed to zeros
        pref_cand2 = studentP(topstudent,:);
        pref_cand2 = pref_cand2(student(topstudent,:));
        topschool = pref_cand2(1);
        
        Q(topschool)= Q(topschool)-1;
        
        % for this top top pair, change this student' preferences indicator to
        % false
        student(topstudent,:)=false;
        
        % also change candidate indicator
        candidate(topstudent) = false;

        %And also clear at tt vector
        
        tt(topstudent)=0;
        
        % do the same for the schools, and change to false when capacity is reduced
        % to zero
        %topschool = studentP(topstudent);
        %%HERE BE CAREFULL! changed place to above
        %pref_cand2 = studentP(topstudent,:);
        %pref_cand2 = pref_cand2(student(topstudent,:));
        %topschool = pref_cand2(1);
        
        %Q(topschool)= Q(topschool)-1

        if Q(topschool) == 0
           school(topschool,:) = false;
            % schoolfor other students, whenever their preferences contain this top school,
            % change to false also

            otherstudent= idstudent(idstudent~=topstudent);

            for ii = otherstudent(1:end)'
                student(ii,studentP(ii,:)==topschool) = false;
            end 
        end 

          % for schools other than this top school, also remove this student from
    % their priority lists
          for kk = 1:m   
              school(kk,schoolP(kk,:)==topstudent) =false;
          end 
          
        end
    end     

end 

end 
