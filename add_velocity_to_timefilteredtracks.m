% TimeFilteredTracks{tr}(:,35)= frame to frame velocity
% TimeFilteredTracks{tr}(:,36)= frame to frame relative velocity
% TimeFilteredTracks{tr}(:,37)= 1h interval velocity
% TimeFilteredTracks{tr}(:,38)= 2h interval velocity
% TimeFilteredTracks{tr}(:,39)= 1h interval relative velocity
% TimeFilteredTracks{tr}(:,40)= 2h interval relative velocity
% TimeFilteredTracks{tr}(:,41)= avg frame to frame relative velocity
% TimeFilteredTracks{tr}(:,42)= avg 1h relative velocity
% TimeFilteredTracks{tr}(:,43)= avg 2h relative velocity
% TimeFilteredTracks{tr}(:,44)= avg aspect ratio over 3 frames
% TimeFilteredTracks{tr}(:,45)= avg aspect ratio over 12 frames (1h)
% TimeFilteredTracks{tr}(:,46)= avg cell shape over 3 frames
% TimeFilteredTracks{tr}(:,47)= avg cell shape over 12 frames (1h)
% TimeFilteredTracks{tr}(:,48)= cell shape change frame to frame
% TimeFilteredTracks{tr}(:,49)= frame to frame volume difference



for tr=1:length(time_filtered_Tracks)
time_filtered_Tracks{tr}(1:end,35:49)=NaN;
end


%% get rid of all zeros and replacing them with NaNs
%if you can replace this with a cleaner easier way, that would be great

delete_zeros = cell(size(time_filtered_Tracks));
for tr=1:length(time_filtered_Tracks) %%setting up the section where zeros will be deleted
    if ~isempty(time_filtered_Tracks{tr})
        if size(time_filtered_Tracks{tr},2) > 19
        delete_zeros{tr}=NaN(size(time_filtered_Tracks{tr},1),31);
        delete_zeros{tr} = time_filtered_Tracks{tr}(:,5:34);
        else 
        delete_zeros{tr}=NaN(size(time_filtered_Tracks{tr},1),15);
        delete_zeros{tr} = time_filtered_Tracks{tr}(:,5:19);
        end
    else
        delete_zeros{tr} = [];
    end
end

delete_zeros = cellfun(@(M) subsasgn(M, substruct('()', {(M()==0)}), NaN), delete_zeros, 'uniform', 0); %% replacing zeros with NaNs

for tr=1:length(time_filtered_Tracks) %% plugging into time_flitered_tracks
    if ~isempty(time_filtered_Tracks{tr})
        if size(delete_zeros{tr},2) > 19
        time_filtered_Tracks{tr}(:,5:34) = delete_zeros{tr};
        else 
        time_filtered_Tracks{tr}(:,5:19) = delete_zeros{tr};
        end
    end
end

clear delete_zeros;

%%

for tr=1:length(time_filtered_Tracks)
    %for t=1:length(time_filtered_Tracks{tr})-1
    if ~isempty(time_filtered_Tracks{tr}) 
        time_filtered_Tracks{tr}(1:end-1,35)= sum((time_filtered_Tracks{tr}(2:end,9:11)-time_filtered_Tracks{tr}(1:end-1,9:11)).^2,2).^.5;
        
        time_filtered_Tracks{tr}(time_filtered_Tracks{tr}(2:end,1)-time_filtered_Tracks{tr}(1:end-1,1)~=1,35)=NaN;
    end
end

        
%% add relative velocity
% assuming you have '"all" Tracks' 'TrackNeighborTracks' which used all Tracks and 'time_filtered_Tracks' with velocities 


for tr=1:length(time_filtered_Tracks)
    %for t=1:length(time_filtered_Tracks{tr})-1
    if ~isempty(time_filtered_Tracks{tr})%dont look at deleted tracks
        
        idx_match=arrayfun( @(tar) find(Tracks{tr}(:,1) == time_filtered_Tracks{tr}(tar,1)),1:size(time_filtered_Tracks{tr},1)); %since the indices are different between the time filtered tracks (which has deleted tracks) and tracks (also TrackNeighborsTracks) we have to match them so the indicies are the same and we can link TrackNeighborsTracks time_filtered_tracks a %tar= time in this array, just a variable
        
        for t=1:length(idx_match)-1
            if isnan(time_filtered_Tracks{tr}(t,35))
                time_filtered_Tracks{tr}(t,36)=NaN;
            else
                NNT=TrackNeighborTracks{tr}{idx_match(t)}; %Nearest Neighbor Track, basically a list of numbers that refrence the index of TrackNeighborTracks within a a certain track(aka tr)
                NNT(arrayfun(@(trn) isempty(time_filtered_Tracks{trn}),NNT))=[]; %marks the neighbours that got deleted in the filtering of the neighnours list. the ending part does the initializing of NNT for trn 
                NNT(NNT==tr)=[]; %get rid of the track itself as a neighbor
                NNT(~arrayfun( @(trn) min(ismember([time_filtered_Tracks{tr}(t,1),time_filtered_Tracks{tr}(t,1)+1],time_filtered_Tracks{NNT(trn)}(:,1))),1:length(NNT)))=[]; %delete the marked neighbours by looking for the min aka zeros and then inverting it so that we can get rid of the 1s instead of the zeros that act as a mark

             %find(time_filtered_Tracks{NNT(1)}(:,1)==time_filtered_Tracks{tr}(t,1))  %fine the time where the track and its neighbor are equalequalize the time between the main track and the neighbor, since each track has its own time and time isnt universal
                
                if isempty(NNT)
                    time_filtered_Tracks{tr}(t,36)=NaN;
                else
                    avg_vel_NNT=mean(cell2mat(arrayfun(@(trn) time_filtered_Tracks{NNT(trn)}(find(time_filtered_Tracks{NNT(trn)}(:,1)==time_filtered_Tracks{tr}(t,1)+1),9:11)-time_filtered_Tracks{NNT(trn)}(find(time_filtered_Tracks{NNT(trn)}(:,1)==time_filtered_Tracks{tr}(t,1)),9:11),1:length(NNT),'UniformOutput',false)')); %calculate the velocity vector diff of each neighbor and then avgs them over all  
                    time_filtered_Tracks{tr}(t,36)=sum((time_filtered_Tracks{tr}(t+1,9:11)-time_filtered_Tracks{tr}(t,9:11)-avg_vel_NNT).^2,2).^.5; %get the relative velocity 
     
                end
            end
        end
    end
end


%% Relative velocity nuclei 

for tr=1:length(time_filtered_Tracks)
    %for t=1:length(time_filtered_Tracks{tr})-1
    if ~isempty(time_filtered_Tracks{tr})%dont look at deleted tracks
        
        idx_match=arrayfun( @(tar) find(Tracks{tr}(:,1) == time_filtered_Tracks{tr}(tar,1)),1:size(time_filtered_Tracks{tr},1)); %since the indices are different between the time filtered tracks (which has deleted tracks) and tracks (also TrackNeighborsTracks) we have to match them so the indicies are the same and we can link TrackNeighborsTracks time_filtered_tracks a %tar= time in this array, just a variable
        
        for t=1:length(idx_match)-1
            if isnan(time_filtered_Tracks{tr}(t,35))
                time_filtered_Tracks{tr}(t,50)=NaN;
            else
                NNT=TrackNeighborTracks{tr}{idx_match(t)}; %Nearest Neighbor Track, basically a list of numbers that refrence the index of TrackNeighborTracks within a a certain track(aka tr)
                NNT(arrayfun(@(trn) isempty(time_filtered_Tracks{trn}),NNT))=[]; %marks the neighbours that got deleted in the filtering of the neighnours list. the ending part does the initializing of NNT for trn 
                NNT(NNT==tr)=[]; %get rid of the track itself as a neighbor
                NNT(~arrayfun( @(trn) min(ismember([time_filtered_Tracks{tr}(t,1),time_filtered_Tracks{tr}(t,1)+1],time_filtered_Tracks{NNT(trn)}(:,1))),1:length(NNT)))=[]; %delete the marked neighbours by looking for the min aka zeros and then inverting it so that we can get rid of the 1s instead of the zeros that act as a mark

             %find(time_filtered_Tracks{NNT(1)}(:,1)==time_filtered_Tracks{tr}(t,1))  %fine the time where the track and its neighbor are equalequalize the time between the main track and the neighbor, since each track has its own time and time isnt universal
                
                if isempty(NNT)
                    time_filtered_Tracks{tr}(t,50)=NaN;
                else
                    avg_vel_NNT=mean(cell2mat(arrayfun(@(trn) time_filtered_Tracks{NNT(trn)}(find(time_filtered_Tracks{NNT(trn)}(:,1)==time_filtered_Tracks{tr}(t,1)+1),2:4)-time_filtered_Tracks{NNT(trn)}(find(time_filtered_Tracks{NNT(trn)}(:,1)==time_filtered_Tracks{tr}(t,1)),2:4),1:length(NNT),'UniformOutput',false)')); %calculate the velocity vector diff of each neighbor and then avgs them over all  
                    time_filtered_Tracks{tr}(t,50)=sum((time_filtered_Tracks{tr}(t+1,2:4)-time_filtered_Tracks{tr}(t,2:4)-avg_vel_NNT).^2,2).^.5; %get the relative velocity 
     
                end
            end
        end
    end
end



%% Velocity every 1h and 2h 


for tr=1:length(time_filtered_Tracks)
    %for t=1:length(time_filtered_Tracks{tr})-1
    if ~isempty(time_filtered_Tracks{tr})
        for tt=1:size(time_filtered_Tracks{tr},1)
            if ismember(time_filtered_Tracks{tr}(tt,1)+12, time_filtered_Tracks{tr}(:,1))
                ttt= find(time_filtered_Tracks{tr}(:,1) == time_filtered_Tracks{tr}(tt,1)+12);
                time_filtered_Tracks{tr}(tt,37)= sum((time_filtered_Tracks{tr}(ttt,9:11)-time_filtered_Tracks{tr}(tt,9:11)).^2,2).^.5;
            else
                time_filtered_Tracks{tr}(tt,37)=NaN;
            end    
            %time_filtered_Tracks{tr}(time_filtered_Tracks{tr}(12:end,1)-time_filtered_Tracks{tr}(1:end-11,1)~=1,37)=NaN;
        end
    end
end


for tr=1:length(time_filtered_Tracks)
    %for t=1:length(time_filtered_Tracks{tr})-1
    if ~isempty(time_filtered_Tracks{tr})
        for tt=1:size(time_filtered_Tracks{tr},1)
            if ismember(time_filtered_Tracks{tr}(tt,1)+24, time_filtered_Tracks{tr}(:,1))
                ttt= find(time_filtered_Tracks{tr}(:,1) == time_filtered_Tracks{tr}(tt,1)+24);
                time_filtered_Tracks{tr}(tt,38)= sum((time_filtered_Tracks{tr}(ttt,9:11)-time_filtered_Tracks{tr}(tt,9:11)).^2,2).^.5;
            else
                time_filtered_Tracks{tr}(tt,38)=NaN;
            end    
            %time_filtered_Tracks{tr}(time_filtered_Tracks{tr}(12:end,1)-time_filtered_Tracks{tr}(1:end-11,1)~=1,37)=NaN;
        end
    end
end



        
%% add relative velocity 1h
% assuming you have '"all" Tracks' 'TrackNeighborTracks' which used all Tracks and 'time_filtered_Tracks' with velocities 


for tr=1:length(time_filtered_Tracks)
    %for t=1:length(time_filtered_Tracks{tr})-1
    if ~isempty(time_filtered_Tracks{tr})%dont look at deleted tracks
        
        idx_match=arrayfun( @(tar) find(Tracks{tr}(:,1) == time_filtered_Tracks{tr}(tar,1)),1:size(time_filtered_Tracks{tr},1)); %since the indices are different between the time filtered tracks (which has deleted tracks) and tracks (also TrackNeighborsTracks) we have to match them so the indicies are the same and we can link TrackNeighborsTracks time_filtered_tracks a %tar= time in this array, just a variable
        
        for t=1:length(idx_match)-1
            if isnan(time_filtered_Tracks{tr}(t,37))
                time_filtered_Tracks{tr}(t,39)=NaN;
            else
                NNT=TrackNeighborTracks{tr}{idx_match(t)}; %Nearest Neighbor Track, basically a list of numbers that refrence the index of TrackNeighborTracks within a a certain track(aka tr)
                NNT(arrayfun(@(trn) isempty(time_filtered_Tracks{trn}),NNT))=[]; %marks the neighbours that got deleted in the filtering of the neighnours list. the ending part does the initializing of NNT for trn 
                NNT(NNT==tr)=[]; %get rid of the track itself as a neighbor
                NNT(~arrayfun( @(trn) min(ismember([time_filtered_Tracks{tr}(t,1),time_filtered_Tracks{tr}(t,1)+12],time_filtered_Tracks{NNT(trn)}(:,1))),1:length(NNT)))=[]; %delete the marked neighbours by looking for the min aka zeros and then inverting it so that we can get rid of the 1s instead of the zeros that act as a mark

             %find(time_filtered_Tracks{NNT(1)}(:,1)==time_filtered_Tracks{tr}(t,1))  %find the time where the track and its neighbor are equal and qualize the time between the main track and the neighbor, since each track has its own time and time isnt universal
                
                if isempty(NNT)
                    time_filtered_Tracks{tr}(t,39)=NaN;
                else
                    avg_vel_NNT=mean(cell2mat(arrayfun(@(trn) time_filtered_Tracks{NNT(trn)}(find(time_filtered_Tracks{NNT(trn)}(:,1)==time_filtered_Tracks{tr}(t,1)+12),9:11)-time_filtered_Tracks{NNT(trn)}(find(time_filtered_Tracks{NNT(trn)}(:,1)==time_filtered_Tracks{tr}(t,1)),9:11),1:length(NNT),'UniformOutput',false)')); %calculate the velocity vector diff of each neighbor and then avgs them over all  
                    
                    time_filtered_Tracks{tr}(t,39)=sum((time_filtered_Tracks{tr}(find(time_filtered_Tracks{tr}(:,1)==time_filtered_Tracks{tr}(t,1)+12),9:11)-time_filtered_Tracks{tr}(find(time_filtered_Tracks{tr}(:,1)==time_filtered_Tracks{tr}(t,1)),9:11)-avg_vel_NNT).^2,2).^.5; %get the relative velocity 
     
                end
            end
        end
    end
end

%% add relative velocity 2h
% assuming you have '"all" Tracks' 'TrackNeighborTracks' which used all Tracks and 'time_filtered_Tracks' with velocities 


for tr=1:length(time_filtered_Tracks)
    %for t=1:length(time_filtered_Tracks{tr})-1
    if ~isempty(time_filtered_Tracks{tr})%dont look at deleted tracks
        
        idx_match=arrayfun( @(tar) find(Tracks{tr}(:,1) == time_filtered_Tracks{tr}(tar,1)),1:size(time_filtered_Tracks{tr},1)); %since the indices are different between the time filtered tracks (which has deleted tracks) and tracks (also TrackNeighborsTracks) we have to match them so the indicies are the same and we can link TrackNeighborsTracks time_filtered_tracks a %tar= time in this array, just a variable
        
        for t=1:length(idx_match)-1
            if isnan(time_filtered_Tracks{tr}(t,38))
                time_filtered_Tracks{tr}(t,40)=NaN;
            else
                NNT=TrackNeighborTracks{tr}{idx_match(t)}; %Nearest Neighbor Track, basically a list of numbers that refrence the index of TrackNeighborTracks within a a certain track(aka tr)
                NNT(arrayfun(@(trn) isempty(time_filtered_Tracks{trn}),NNT))=[]; %marks the neighbours that got deleted in the filtering of the neighnours list. the ending part does the initializing of NNT for trn 
                NNT(NNT==tr)=[]; %get rid of the track itself as a neighbor
                NNT(~arrayfun( @(trn) min(ismember([time_filtered_Tracks{tr}(t,1),time_filtered_Tracks{tr}(t,1)+24],time_filtered_Tracks{NNT(trn)}(:,1))),1:length(NNT)))=[]; %delete the marked neighbours by looking for the min aka zeros and then inverting it so that we can get rid of the 1s instead of the zeros that act as a mark

             %find(time_filtered_Tracks{NNT(1)}(:,1)==time_filtered_Tracks{tr}(t,1))  %fine the time where the track and its neighbor are equalequalize the time between the main track and the neighbor, since each track has its own time and time isnt universal
                
                if isempty(NNT)
                    time_filtered_Tracks{tr}(t,40)=NaN;
                else
                    avg_vel_NNT=mean(cell2mat(arrayfun(@(trn) time_filtered_Tracks{NNT(trn)}(find(time_filtered_Tracks{NNT(trn)}(:,1)==time_filtered_Tracks{tr}(t,1)+24),9:11)-time_filtered_Tracks{NNT(trn)}(find(time_filtered_Tracks{NNT(trn)}(:,1)==time_filtered_Tracks{tr}(t,1)),9:11),1:length(NNT),'UniformOutput',false)')); %calculate the velocity vector diff of each neighbor and then avgs them over all  
                    
                    time_filtered_Tracks{tr}(t,40)=sum((time_filtered_Tracks{tr}(find(time_filtered_Tracks{tr}(:,1)==time_filtered_Tracks{tr}(t,1)+24),9:11)-time_filtered_Tracks{tr}(find(time_filtered_Tracks{tr}(:,1)==time_filtered_Tracks{tr}(t,1)),9:11)-avg_vel_NNT).^2,2).^.5; %get the relative velocity 
     
                end
            end
        end
    end
end

%% add running averages

for tr=1:length(time_filtered_Tracks)
    if ~isempty(time_filtered_Tracks{tr}) && size(time_filtered_Tracks{tr},2)>38
        intermed=[];
        intermed(:,1)=time_filtered_Tracks{tr}(1,1):time_filtered_Tracks{tr}(end,1);
        for t=1:size(intermed,1)
            if ismember(intermed(t,1),time_filtered_Tracks{tr}(:,1)) % is t filtered out
                intermed(t,2:6)=time_filtered_Tracks{tr}(time_filtered_Tracks{tr}(:,1)==intermed(t,1),[36,39,40,6,7]);
            else%% put NaNs in
                intermed(t,2:6)=[NaN,NaN,NaN,NaN,NaN];    
            end    
        end
         intermed(:,7)=movmean(intermed(:,2),12,'omitnan'); %avg frane to frame rel vel
         intermed(:,8)=movmean(intermed(:,3),12,'omitnan'); %avg 1h rel vel
         intermed(:,9)=movmean(intermed(:,4),12,'omitnan'); %avg 2h rel vel
         intermed(:,10)=movmean(intermed(:,5),12,'omitnan');%avg ar 1h
         intermed(:,11)=movmean(intermed(:,5),3,'omitnan'); %avg ar 3 frames
         intermed(:,12)=movmean(intermed(:,6),12,'omitnan');%avg shape 1h
         intermed(:,13)=movmean(intermed(:,6),3,'omitnan'); %avg shape 3 frames
         intermed(1:end-1,14)=intermed(2:end,13)-intermed(1:end-1,13); %change of shape
         if intermed(end,14)==0
            intermed(end,14)=NaN;
         end
         
        if size(time_filtered_Tracks{tr}(:,41:48))~=size(intermed(:,7:14))
            intermed(isnan(intermed(1:end-1,5)),:)=[];
        end
        
        time_filtered_Tracks{tr}(:,41:48)= intermed(:,7:14);
    end
end

%% add frame to frame volume difference

for tr=1:length(time_filtered_Tracks)
    %for t=1:length(time_filtered_Tracks{tr})-1
    if ~isempty(time_filtered_Tracks{tr}) 
        time_filtered_Tracks{tr}(2:end,49)= time_filtered_Tracks{tr}(2:end,8)-time_filtered_Tracks{tr}(1:end-1,8);
        
        time_filtered_Tracks{tr}(time_filtered_Tracks{tr}(2:end,1)-time_filtered_Tracks{tr}(1:end-1,1)~=1,49)=NaN;
    end
end

%% setting last input to nan in most columns
for tr=1:length(time_filtered_Tracks)
if ~isempty(time_filtered_Tracks{tr}) 
time_filtered_Tracks{tr}(end,35:49)=NaN;
end
end
        