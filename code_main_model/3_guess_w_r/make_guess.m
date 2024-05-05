clear all;

cd ..; load('1_solve_ss\computation_results\data.mat','w_iterate','r_iterate','pdf_B','nmkt','nz','nb_pdf','ntheta_pdf'); cd 3_guess_w_r;
period_trans=200;

r_end_index=min(find(r_iterate==0))-1;
w_end_index(1:nmkt,1)=min(find(w_iterate(r_end_index,:,1)==0))-1;

r(1:period_trans,1)=r_iterate(r_end_index);
w(1:period_trans,1:nmkt)=nan;
for i=1:nmkt
    w(:,i)=w_iterate(r_end_index,w_end_index,i);
end

save('ini.mat','pdf_B','w','r');

%% write w
fileID = fopen('w.dat','w');
print2(fileID,w,period_trans,nmkt);

%% write r
fileID = fopen('r.dat','w');
print1(fileID,r,period_trans);

%% write pdf_ini
fileID = fopen('pdf_ini.dat','w');
print4(fileID,pdf_ini,nmkt,nz,nb_pdf,ntheta_pdf);