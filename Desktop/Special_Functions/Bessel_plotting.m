%% S.Aksimsek, 2011
% The Graph of Bessel Function

z=0:0.01:20;
for v=0:6
    for c=1:length(z);
        g(v+1,c)=BesselFunction(v,z(c));
    end
end
figure
plot(z,g');grid;
% axis([0 2000 -0.5 1 ])
title('Bessel Function')
ylabel('J_\nu(z)','fontsize',18,'fontweight','b')
xlabel('z','fontsize',18,'fontweight','b')
legend('J_0','J_1','J_2','J_3','J_4','J_5','J_6')
set(gca,'FontName','Times New Roman','FontSize',24)

z=0:0.01:20;
for v=0:6
    for c=1:length(z);
        g(v+1,c)=NeumannFunction(v,z(c));
    end
end
figure
plot(z,g');grid;
axis([0 12 -1.5 1 ])
%axis([0 2000 -1.5 0.5 ])
title('Neumann Function')
ylabel('Y_\nu(z)','fontsize',18,'fontweight','b')
xlabel('z','fontsize',18,'fontweight','b')
legend('Y_0','Y_1','Y_2','Y_3','Y_4','Y_5','Y_6')
set(gca,'FontName','Times New Roman','FontSize',24)



