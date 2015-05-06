function mu = mean_shape_procrustes(fnt)
% calculate mean shape mu from input data fnt.
% fnt: mxn complex matrix. m is the number of repetition. n is the number of landmarks.  

num_data = size(fnt,1);
num_landmark = size(fnt,2);

fnt = fnt - repmat(mean(fnt,2),1,num_landmark);
fnt = fnt ./ repmat(sqrt(sum(fnt.*conj(fnt),2)),1,num_landmark);
S = zeros(num_landmark,num_landmark);
for ii = 1:num_data 
    S = S + fnt(ii,:).'*conj(fnt(ii,:));
end
[U S1 V] = svds(S,10);
mu = U(:,1);
