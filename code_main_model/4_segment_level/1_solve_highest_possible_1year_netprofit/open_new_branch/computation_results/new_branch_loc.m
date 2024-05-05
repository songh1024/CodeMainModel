clear all;

load('data.mat','bank_profit','nmkt')

profit=sum(bank_profit,2);

% A
cd ..; cd ..; load('no_new_branch\computation_results\data.mat','bank_profit_total');
load('data.mat','d')
profit(d(1,:)==0)=0;
profit(profit>0)=profit(profit>0)-bank_profit_total;

% B
% profit(profit>0)=profit(profit>0)-profit(1); % profit(1) is the benchmark where no new branch is opened

cd ..; cd ..; cd ..; cd ..; 
load('comm_loc.mat','comm_loc');
load('2_partition\segment.mat','segment');
cd 4_segment_level\1_solve_highest_possible_1year_netprofit\open_new_branch\computation_results;

comm_loc=comm_loc(:,1:11);
bank86=find(comm_loc(:,1)==1);
banknew=find(comm_loc(:,1)==0 & comm_loc(:,11)==1);

seg=unique(segment(:,2));
nseg=length(seg);
seg_nmkt(1:nseg,1)=nan;
seg_n86bank(1:nseg,1)=nan;
seg_nnewbank_data(1:nseg,1)=nan;
for i_seg=1:nseg
    segmkt=find(segment(:,2)==seg(i_seg));
    seg_nmkt(i_seg)=length(segmkt);
    seg_n86bank(i_seg)=length(intersect(segmkt,bank86));
    seg_nnewbank_data(i_seg)=length(intersect(segmkt,banknew));
%     eval(['segmkt',num2str(i_seg),'=segmkt;']);
end
seg_nnewbank_bound(:,1)=max(seg_nnewbank_data-10,0);
seg_nnewbank_bound(:,2)=min(seg_nnewbank_data+10,seg_nmkt-seg_n86bank);
save('seg_nnewbank_bound.mat','seg_nnewbank_bound','comm_loc','segment');


segment=[segment,profit];

% for each segment, rank the profit, and see how many are group 1
seg_g1top(1:nseg,1)=0;
seg_open(1:nseg,1:nmkt)=-1; % > 0 for top group 1; = 0 for 86 branch; = -1 for those unsure; = -2 for markets not in this segment

for i_seg=1:nseg
    segmkt=find(segment(:,2)==seg(i_seg));
    seg_open(i_seg,setdiff(segment(:,1),segmkt))=-2;
    seg_open(i_seg,segment(:,2)==i_seg & segment(:,3)==0)=0;

    seg_profit=segment(segmkt,:);
    seg_profit=sortrows(seg_profit,4,'descend');

    seg_g1top(i_seg)=find(seg_profit(:,3) ~= 1, 1 ,'first')-1;
    seg_open(i_seg,seg_profit(1:seg_g1top(i_seg),1)) = 1:1:seg_g1top(i_seg);
end

seg_nnewbank_bound(:,1)=seg_nnewbank_bound(:,1)-seg_g1top;

seg_nnewbank_bound(:,2)=seg_nnewbank_bound(:,2)-seg_g1top;

seg_nnewbank_bound(seg_nnewbank_bound<0)=0;
save('seg_open.mat','seg_open');

fileID = fopen('combination_input\seg.dat','w');
print_int1(fileID,segment(:,2)-1,nmkt);
fileID = fopen('combination_input\seg_open.dat','w');
print_int2(fileID,seg_open(:,1:nmkt),nseg,nmkt);
fileID = fopen('combination_input\seg_nnewbank_bound.dat','w');
print_int2(fileID,seg_nnewbank_bound,nseg,2);
