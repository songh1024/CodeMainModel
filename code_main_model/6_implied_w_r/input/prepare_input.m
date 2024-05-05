clear all;

load('input.mat')

%% write Pi
fileID = fopen('Pi.dat','w');
print1(fileID,Pi,nmkt);

%% write Z
fileID = fopen('Z.dat','w');
print1(fileID,Z,nmkt);

%% write d
cd ..; cd ..; load('5_country_level\computation_results\data.mat','bank_opt');
cd ..; load('data\comm_loc.mat'); comm_loc=comm_loc(:,1:11);
cd code_main_model\6_implied_w_r\input;
nyears=11;
comm_loc_model(1:nyears,1:nmkt)=0;
d(1:nyears,1:nmkt)=0;
for i_t=1:nyears
    for i_mkt=1:nmkt
        d(i_t,i_mkt)=1e10;
        for j_mkt=1:nmkt
            if (bank_opt(j_mkt)<i_t && (Pi(j_mkt)<h || i_mkt==j_mkt) )
                d(i_t,i_mkt)=min(d(i_t,i_mkt),distance(i_mkt,j_mkt));
            end
        end
    end
end


fileID = fopen('d.dat','w');
print2(fileID,d,nyears,nmkt);

comm_loc_model=comm_loc_model';
save('comm_loc_model.mat','comm_loc_model','comm_loc');