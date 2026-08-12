load u0b
load u0c
load u0d
load u0e
load u0f
load u0g
load u0h
load u0i
load u0k

load u1b
load u1c
load u1d
load u1e
load u1f
load u1g
load u1h
load u1i
load u1k

%figure(1);
%plot( u0d(:,1), u0d(:,2),'k', u0a(:,1), u0a(:,3),'r', u0b(:,1), u0b(:,3),'m', u0c(:,1), u0c(:,3),'c', u0d(:,1), u0d(:,3),'b' );
%axis([0,20,-0.12,0.28]);

%figure(2);
%plot( u1d(:,1), u1d(:,2),'k', u1a(:,1), u1a(:,3),'r', u1b(:,1), u1b(:,3),'m', u1c(:,1), u1c(:,3),'c', u1d(:,1), u1d(:,3),'b' );
%axis([0,20,-0.12,0.28]);

ord = [3,5,7,9,11,13,15,17];

err0(1) = max( abs( u0b(:,3) - u0k(:,3)) );
err0(2) = max( abs( u0c(:,3) - u0k(:,3)) );
err0(3) = max( abs( u0d(:,3) - u0k(:,3)) );
err0(4) = max( abs( u0e(:,3) - u0k(:,3)) );
err0(5) = max( abs( u0f(:,3) - u0k(:,3)) );
err0(6) = max( abs( u0g(:,3) - u0k(:,3)) );
err0(7) = max( abs( u0h(:,3) - u0k(:,3)) );
err0(8) = max( abs( u0i(:,3) - u0k(:,3)) );


err1(1) = max( abs( u1b(:,3) - u1k(:,3)) );
err1(2) = max( abs( u1c(:,3) - u1k(:,3)) );
err1(3) = max( abs( u1d(:,3) - u1k(:,3)) );
err1(4) = max( abs( u1e(:,3) - u1k(:,3)) );
err1(5) = max( abs( u1f(:,3) - u1k(:,3)) );
err1(6) = max( abs( u1g(:,3) - u1k(:,3)) );
err1(7) = max( abs( u1h(:,3) - u1k(:,3)) );
err1(8) = max( abs( u1i(:,3) - u1k(:,3)) );

semilogy( ord, err0, 'k+-', ord, err1, 'ko-', 'LineWidth', 1 );
axis([3,17,1e-8,1e-1]);
legend('M0-T0','M1-T1');
xticks(ord)

ax = gca;
ax.FontSize = 12;

