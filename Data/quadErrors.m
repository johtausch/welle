load q01;
load q03;
load q05;
load q07;
load q09;
load q11;
load q13;
load q15;
load q17;
load q19;
load q21;

n = length(q03); % the other data files have the same length

i1 = 0;
i2 = 0;
i3 = 0;
i4 = 0;
i5 = 0;
i6 = 0;

for  i=1:n
  if ( q03(i,1) == 1 )
    i1 = i1 + 1;
    err1(i1,:) = [q01(i,2) q03(i,2) q05(i,2) q07(i,2) q09(i,2) q11(i,2) q13(i,2) q15(i,2) q17(i,2) q19(i,2)];
    err1(i1,:) = abs(err1(i1,:) - q21(i,2));
  elseif ( q03(i,1) == 2 )
    i2 = i2 + 1;
    err2(i2,:) = [q01(i,2) q03(i,2) q05(i,2) q07(i,2) q09(i,2) q11(i,2) q13(i,2) q15(i,2) q17(i,2) q19(i,2)];
    err2(i2,:) = abs(err2(i2,:) - q21(i,2));
  elseif ( q03(i,1) == 3 )
    i3 = i3 + 1;
    err3(i3,:) = [q01(i,2) q03(i,2) q05(i,2) q07(i,2) q09(i,2) q11(i,2) q13(i,2) q15(i,2) q17(i,2) q19(i,2)];
    err3(i3,:) = abs(err3(i3,:) - q21(i,2));
  elseif ( q03(i,1) == 4 )
    i4 = i4 + 1;
    err4(i4,:) = [q01(i,2) q03(i,2) q05(i,2) q07(i,2) q09(i,2) q11(i,2) q13(i,2) q15(i,2) q17(i,2) q19(i,2)];
    err4(i4,:) = abs(err4(i4,:) - q21(i,2));
  elseif ( q03(i,1) == 5 )
    i5 = i1 + 1;
    err5(i5,:) = [q01(i,2) q03(i,2) q05(i,2) q07(i,2) q09(i,2) q11(i,2) q13(i,2) q15(i,2) q17(i,2) q19(i,2)];
    err5(i5,:) = abs(err5(i5,:) - q21(i,2));
  elseif ( q03(i,1) == 6 )
    i6 = i6 + 1;
    err6(i6,:) = [q01(i,2) q03(i,2) q05(i,2) q07(i,2) q09(i,2) q11(i,2) q13(i,2) q15(i,2) q17(i,2) q19(i,2)];
    err6(i6,:) = abs(err6(i6,:) - q21(i,2));
  else
    disp("unknown type");
  end
end
  
max1 = max(err1); 
max2 = max(err2); 
max3 = max(err3); 
max4 = max(err4); 
max5 = max(err5); 
max6 = max(err6); 

ord = [1,3,5,7,9,11,13,15,17,19];

figure(1)
semilogy(ord,max1,'k+-', ord,max2,'k-o', ord,max3,'ks-', 'LineWidth', 1 );
axis([1,19,1e-16,1e-2]);
yticks( [1e-15 1e-13 1e-11 1e-9 1e-7 1e-5 1e-3]);
xticks(ord)
legend('regular apices','singular apices','mixed apices');


figure(2)
semilogy(ord,max4,'k+-', ord,max5,'ko-', ord,max6,'ks-', 'LineWidth', 1 ) ;
axis([1,19,1e-16,1e-2]);
yticks( [1e-15 1e-13 1e-11 1e-9 1e-7 1e-5 1e-3]);
legend('regular apices','singular apices','mixed apices');
xticks(ord)

ax = gca;
ax.FontSize = 12;





