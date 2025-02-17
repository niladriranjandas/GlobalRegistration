
function a = plotReadings(rd)

 x=[0 0.002 0.004 0.006 0.008 0.010];
hold on ; 
title('Comparison of Methods under Noisy consitions');
xlabel('Noise');
ylabel('RMS Error');
ICP=plot(x,rd(2,2:7)','r');
SDP=plot(x,rd(3,2:7)','b'); 
EIG=plot(x,rd(4,2:7)','m');
hold off
legend([ICP EIG SDP],'Sequential-ICP','EIG-SYNC','SDP-SYNC');

end