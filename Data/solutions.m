load s0a
load s0b
load s0c
load s0d
load s0e

load s1a
load s1b
load s1c
load s1d
load s1e

load s2a
load s2b
load s2c
load s2d
load s2e

load s3a
load s3b
load s3c
load s3d
load s3e

%figure(1);
%plot( s0d(:,1), s0d(:,2),'k',... 
%      s0a(:,1), s0a(:,3),'r', s0b(:,1), s0b(:,3),'m', s0c(:,1), s0c(:,3),'g', s0d(:,1), s0d(:,3),'c' , s0e(:,1), s0e(:,3),'b' );
%axis([0,20,-0.12,0.28]);
%
%figure(2);
%plot( s1d(:,1), s1d(:,2),'k',... 
%      s1a(:,1), s1a(:,3),'r', s1b(:,1), s1b(:,3),'m', s1c(:,1), s1c(:,3),'g', s1d(:,1), s1d(:,3),'c' , s1e(:,1), s1e(:,3),'b' );
%axis([0,20,-0.12,0.28]);
%
%figure(3);
%plot( s2d(:,1), s2d(:,2),'k',... 
%      s2a(:,1), s2a(:,3),'r', s2b(:,1), s2b(:,3),'m', s2c(:,1), s2c(:,3),'g', s2d(:,1), s2d(:,3),'c' , s2e(:,1), s2e(:,3),'b' );
%axis([0,20,-0.12,0.28]);
%
%figure(4);
%plot( s3d(:,1), s3d(:,2),'k',... 
%      s3a(:,1), s3a(:,3),'r', s3b(:,1), s3b(:,3),'m', s3c(:,1), s3c(:,3),'g', s3d(:,1), s3d(:,3),'c' , s3e(:,1), s3e(:,3),'b' );
%axis([0,20,-0.12,0.28]);

plot(s3d(:,1), s3d(:,2),'k',... 
     s0b(:,1), s0b(:,3),'r', s1c(:,1), s1c(:,3),'m', s2d(:,1), s2d(:,3),'c', s3e(:,1), s3e(:,3),'b' );
axis([0,20,-0.12,0.28]);
legend('exact','M0-T1','M1-T2','M2-T3', 'M3-T4');


ax = gca;
ax.FontSize = 12;


%cols=space  rows=time  

% maximum errors
err(1,1) = max( abs(s0a(:,2)-s0a(:,3)));
err(2,1) = max( abs(s0b(:,2)-s0b(:,3)));
err(3,1) = max( abs(s0c(:,2)-s0c(:,3)));
err(4,1) = max( abs(s0d(:,2)-s0d(:,3)));
err(5,1) = max( abs(s0e(:,2)-s0e(:,3)));
       
err(1,2) = max( abs(s1a(:,2)-s1a(:,3)));
err(2,2) = max( abs(s1b(:,2)-s1b(:,3)));
err(3,2) = max( abs(s1c(:,2)-s1c(:,3)));
err(4,2) = max( abs(s1d(:,2)-s1d(:,3)));
err(5,2) = max( abs(s1e(:,2)-s1e(:,3)));
       
err(1,3) = max( abs(s2a(:,2)-s2a(:,3)));
err(2,3) = max( abs(s2b(:,2)-s2b(:,3)));
err(3,3) = max( abs(s2c(:,2)-s2c(:,3)));
err(4,3) = max( abs(s2d(:,2)-s2d(:,3)));
err(5,3) = max( abs(s2e(:,2)-s2e(:,3)));
       
err(1,4) = max( abs(s3a(:,2)-s3a(:,3)));
err(2,4) = max( abs(s3b(:,2)-s3b(:,3)));
err(3,4) = max( abs(s3c(:,2)-s3c(:,3)));
err(4,4) = max( abs(s3d(:,2)-s3d(:,3)));
err(5,4) = max( abs(s3e(:,2)-s3e(:,3)));


% additional information about the computation of the matrix entries
% to get this information run welle with the flags j=2 and -v=1 
%

% number of nonzeros in all A[d] matrices
nzA=[
4.900800e+04  4.411680e+05  4.144224e+06  4.263605e+07  
9.398400e+04  8.124960e+05  7.124352e+06  4.263605e+07  
1.836960e+05  1.549968e+06  1.308216e+07  1.142695e+08  
3.624960e+05  3.025008e+06  2.498136e+07  2.097090e+08  
7.204800e+05  5.974608e+06  4.877976e+07  4.005943e+08
]; 

% number of times a $B_d$ is calculated
nB =[
5.236800e+04  5.099040e+05  5.303904e+06  6.137947e+07 
9.734400e+04  8.812320e+05  8.284032e+06  6.137947e+07 
1.870560e+05  1.618704e+06  1.424184e+07  1.330640e+08 
3.658560e+05  3.093744e+06  2.614104e+07  2.285035e+08 
7.238400e+05  6.043344e+06  4.993944e+07  4.193888e+08 
];
             
% total number of convex hull integrals that have been computed
nI = [
2.093760e+05  1.726752e+06  1.438776e+07  1.249805e+08 
4.206560e+05  3.377184e+06  2.760768e+07  1.249805e+08 
8.362400e+05  6.661296e+06  5.413397e+07  4.445256e+08 
1.651040e+06  1.321800e+07  1.069501e+08  8.702290e+08 
3.280976e+06  2.630429e+07  2.125647e+08  1.721697e+09 
];
                             
% cpu time to setup all A[d] matrices
cpuA = [
4.715125e-01  3.407511e+00  2.468513e+01 2.392575e+02 
9.022170e-01  8.251651e+00  6.108267e+01 2.367662e+02 
1.813027e+00  1.586463e+01  1.294217e+02 1.040855e+03 
3.641324e+00  3.097609e+01  2.521483e+02 1.959415e+03 
6.838564e+00  5.621650e+01  4.863241e+02 3.849451e+03 
];
                       
% cpu time of block forward solver
cpuS = [
3.850184e-03  3.980307e+00  1.611002e+01 3.929237e+02
1.575495e-02  3.988405e+00  3.679882e+01 3.970349e+02
6.260416e-02  5.867530e+00  5.583610e+01 2.138550e+02
3.787725e-01  1.520235e+01  1.261528e+02 6.609430e+02
9.620714e-01  4.790974e+01  3.908093e+02 2.740853e+03
];


%number of integral / nonzeros
ratio = nI./nzA;
