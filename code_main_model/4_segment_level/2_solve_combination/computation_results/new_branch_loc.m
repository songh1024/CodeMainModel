clear all;

cd ..; cd ..;
load('1_solve_highest_possible_1year_netprofit\open_new_branch\computation_results\seg_nnewbank_bound.mat')
load('1_solve_highest_possible_1year_netprofit\open_new_branch\computation_results\seg_open.mat')
cd 2_solve_combination\computation_results;
nmkt=size(comm_loc,1);
nyear=size(comm_loc,2);

bank86=find(comm_loc(:,1)==1);
seg=unique(segment(:,2));
nseg=length(seg);
seg_year_newbranch(1:nseg,1:nyear)=0;
for i_seg=1:nseg
    for i_year=1:nyear-1
        seg_year_newbranch(i_seg, i_year+1)=length(find(segment(:,2)==seg(i_seg) & comm_loc(:,i_year)==0 & comm_loc(:,i_year+1)==1));
    end
end
nnewbank=diff(sum(comm_loc));
fileID = fopen('knapsack_input\nnewbank.dat','w');
print_int1(fileID,[0,nnewbank],nyear);
fileID = fopen('knapsack_input\seg_year_newbranch.dat','w');
print_int2(fileID,seg_year_newbranch',nyear,nseg);
fileID = fopen('knapsack_input\seg_nnewbank_bound.dat','w');
print_int2(fileID,seg_nnewbank_bound,nseg,2);

% read location and profit
for i_seg=1:nseg
    eval(['fileID = fopen(''location' num2str(i_seg) '.txt'',''r'');']);
    eval(['location' num2str(i_seg) '=fscanf(fileID,''%f,'',[1, inf]);']);
    eval(['ncase=length(location' num2str(i_seg) ')/nmkt;']);
    eval(['location' num2str(i_seg) '=read2(location' num2str(i_seg) ', ncase, nmkt);']);

    eval(['fileID = fopen(''profit' num2str(i_seg) '.txt'',''r'');']);
    eval(['profit' num2str(i_seg) '=fscanf(fileID,''%f,'',[1, inf]);']);
end


segment_openorder(1:nseg, 1:nmkt)=0;
for i_seg=1:nseg
    temp=seg_open(i_seg,:);
    ttemp=unique(temp);
    n=max(ttemp); % top group 1
    if n>0
        for i=1:n
            segment_openorder(i_seg,i)=find(seg_open(i_seg,:)==i);
        end
    end
    
    cum=n;
    if n<seg_nnewbank_bound(i_seg,2)
        eval(['location = location' num2str(i_seg) ';']);
        ncase=size(location,1);
        for i_case=2:ncase
            open_index=find(location(i_case,:)==1);
            new_open_index=open_index(~ismember(open_index,segment_openorder(i_seg,1:cum)));
            segment_openorder(i_seg,cum+1)=new_open_index;
            cum=cum+1;
        end

    end
    
end


fileID = fopen('knapsack_input\bank86.dat','w');
print_int1(fileID,comm_loc(:,1),nmkt);
fileID = fopen('knapsack_input\seg.dat','w');
print_int1(fileID,segment(:,2)-1,nmkt);
fileID = fopen('knapsack_input\segment_openorder.dat','w');
print_int2(fileID,segment_openorder(:,1:nmkt),nseg,nmkt);