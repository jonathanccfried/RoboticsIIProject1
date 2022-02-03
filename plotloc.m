%plot trajectory on the ground

figure(1)
plot(qdlist(1,:),qdlist(2,:),qlist(1,:),qlist(2,:),qhatKFlist(1,:),qhatKFlist(2,:),qhatlist(1,:),qhatlist(2,:),qhatNLlist(1,:),qhatNLlist(2,:))
grid on
xlabel('x')
ylabel('y')
title('traj')
legend('qd','q','qKF','qhat','qNL','Location','southeast')

%plot orientation across time

figure(2)
plot(tlist,qdlist(3,:),tlist,qlist(3,:),tlist,qhatKFlist(3,:),tlist,qhatlist(3,:),tlist,qhatNLlist(3,:))
grid on
xlabel('t')
ylabel('theta')
title('trajtheta')
legend('thetad','theta','thetaKF','thetahat','thetaNL','Location','southwest')


%calculate least square error

errhat = immse(qlist,qhatlist);
errNL = immse(qlist,qhatNLlist);
errKF = immse(qlist,qhatKFlist);