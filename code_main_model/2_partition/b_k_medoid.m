clear all;
load('input.mat','distance');
load('initial_guess.mat');

n_guess=1e6;
max_it=20;
n_clus=35;
n_loc=1428;
D=30;

location_index_final(1:n_guess,1:n_loc)=0.0;
mediod(1:max_it,1:n_clus)=0.0;
mediod_final(1:n_guess,1:n_clus)=0.0;
loss_final(1:n_guess)=0.0;

tic;
for i=1:n_guess
    loc_index(1:max_it,1:n_loc)=0.0;
    mediod(1:max_it+1,1:n_clus)=0.0;
    it=1;
    mediod(it,:)=guess_mediod(i,:);
    for it=1:max_it
        loc_index(it,mediod(it,:))=1:1:n_clus;
        % assign
        dist_temp=distance(:,mediod(it,:));
        [~,index]=min(dist_temp');
        loc_index(it,:)=index;

        % in each partition, find new medoid
        for i_par=1:n_clus
            loc_temp=find(loc_index(it,:)==i_par);
            n_loc_temp=length(loc_temp);
            dist_temp=distance(loc_temp,loc_temp);
            loss=sum(dist_temp);
            [~,index]=min(loss);
            mediod(it+1,i_par)=loc_temp(index);
        end
        if (sum(mediod(it+1,:)==mediod(it,:))==n_clus)
            location_index_final(i,:)=loc_index(it,:);
            mediod_final(i,:)=mediod(it,:);
            loss_temp(1:n_clus,1:n_clus)=0.0;
            for ii=1:n_clus
                for jj=1:n_clus
                    loc_temp_ii=find(location_index_final(i,:)==ii);
                    loc_temp_jj=find(location_index_final(i,:)==jj);
                    dist_temp=distance(loc_temp_ii,loc_temp_jj);
                    loss_temp(ii,jj)=min(dist_temp,[],'all');
                end
            end
            a(1:n_clus)=0.0;
            for ii=1:n_clus
                temp=loss_temp(ii,:);
                a(ii)=min(temp(temp~=0));
            end
            loss_final(i)=sum(a<D);
            break;
        end
    end
end
toc

[~,index]=min(loss_final);
location_index_final=location_index_final(index,:);
mediod_final=mediod_final(index,:);

save('final_choice.mat','location_index_final','mediod_final')