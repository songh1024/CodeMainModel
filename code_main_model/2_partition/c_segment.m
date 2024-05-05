clear all;

h=150;

load('final_choice.mat');
load('input.mat','distance','Pi');
pop_den=Pi;

nmkt=length(location_index_final);
seg=unique(location_index_final);
nseg=length(seg);

segment(1:nmkt,1:3)=nan; % market index, segment, group
segment(:,1)=1:1:nmkt;
segment(:,2)=location_index_final';
segment(:,3)=2; % the default is group 2

% pop den > h
segment(pop_den>150,3)=1;
for i_seg=1:nseg
    seg_mkt=find(segment(:,2)==seg(i_seg));

    % d > 50
    dist=distance(seg_mkt,seg_mkt);
    dist(dist==0)=1e7; % get rid of diagonal
    dist_closest=min(dist);
    segment(seg_mkt(dist_closest>50),3)=1; % pop den > h
end

load('comm_loc.mat','comm_loc');
bank86=find(comm_loc(:,1)==1); % call markets with bank in 86 group 0
segment(bank86,3)=0;

save('segment.mat','segment')
