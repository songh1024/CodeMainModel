clear all;
%% Generate initial guess 
n_guess=1e6;
n_partition=35;
n_location=1428;
guess_mediod(1:n_guess,1:n_partition)=0;
for i=1:n_guess
    guess_mediod(i,:)=randsample(n_location,n_partition);
end

save('initial_guess.mat','guess_mediod');