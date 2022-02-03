%plot landmark trajectory on the ground

figure(1)
if nb+nr==4
plot(pL(1,1),pL(2,1),'o',pL(1,2),pL(2,2),'o',pL(1,3),pL(2,3),'o',pB(1,1),pB(2,1),'o',LStatelist(1,:),LStatelist(2,:),LStatelist(3,:),LStatelist(4,:),LStatelist(5,:),LStatelist(6,:),LStatelist(7,:),LStatelist(8,:))
grid on
xlabel('x')
ylabel('y')
title('landrepo')
legend('pL1','pL2','pL3','pB1','hpL1','hpL2','hpL3','hpB1','Location','southeastoutside')
elseif nb+nr==6
plot(pL(1,1),pL(2,1),'o',pL(1,2),pL(2,2),'o',pL(1,3),pL(2,3),'o',pL(1,4),pL(2,4),'o',pB(1,1),pB(2,1),'o',pB(1,2),pB(2,2),'o',LStatelist(1,:),LStatelist(2,:),LStatelist(3,:),LStatelist(4,:),LStatelist(5,:),LStatelist(6,:),LStatelist(7,:),LStatelist(8,:),LStatelist(9,:),LStatelist(10,:),LStatelist(11,:),LStatelist(12,:))
grid on
xlabel('x')
ylabel('y')
title('landrepo')
legend('pL1','pL2','pL3','pL4','pB1','pB2','hpL1','hpL2','hpL3','hpL4','hpB1','hpB2','Location','southeastoutside')
else
plot(pL(1,1),pL(2,1),'o',pL(1,2),pL(2,2),'o',pL(1,3),pL(2,3),'o',pL(1,4),pL(2,4),'o',pL(1,5),pL(2,5),'o',pB(1,1),pB(2,1),'o',pB(1,2),pB(2,2),'o',pB(1,3),pB(2,3),'o',LStatelist(1,:),LStatelist(2,:),LStatelist(3,:),LStatelist(4,:),LStatelist(5,:),LStatelist(6,:),LStatelist(7,:),LStatelist(8,:),LStatelist(9,:),LStatelist(10,:),LStatelist(11,:),LStatelist(12,:),LStatelist(13,:),LStatelist(14,:),LStatelist(15,:),LStatelist(16,:))
grid on
xlabel('x')
ylabel('y')
title('landrepo')
legend('pL1','pL2','pL3','pL4','pL5','pB1','pB2','pB3','hpL1','hpL2','hpL3','hpL4','hpL5','hpB1','hpB2','hpB3','Location','southeastoutside')
end

%plot landmark error across time

figure(2)
plot(tlist,LStateTrue*ones(size(tlist))-LStatelist)
grid on
xlabel('t')
ylabel('e')
title('lande')
if nb+nr ==4
legend('eL_1','eL_2','eL_3','eB_1','Location','southwest')
elseif nb+nr == 6
legend('eL_1','eL_2','eL_3','eL_4','eB_1','eB_2','Location','southwest')
else
legend('eL_1','eL_2','eL_3','eL_4','eL_5','eB_1','eB_2','eB_3','Location','southwest')
end

%calculate least square error

errland = immse(LStateTrue*ones(size(tlist)),LStatelist);
