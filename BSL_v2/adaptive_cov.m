function covariance = adaptive_cov(thetas,C,n)
    sd = 4*2.4^2/4;
    epsilon = 0;
    mean_old = mean(thetas(:,1:end-1),2);
    mean_new = mean(thetas(:,1:end),2); 
    covariance = ((n-1)/n)*C+(sd/n)*(n*mean_old*mean_old'-(n+1)*mean_new*mean_new' + thetas(:,end)*thetas(:,end)'+epsilon*eye(4));
end