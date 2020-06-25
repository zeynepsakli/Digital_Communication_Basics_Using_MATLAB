clear all
close all
clc
to = 2000* (0: (200)*10 - 1) / 10000;
hStr = RandStream('mt19937ar', 'Seed', 0);
x = 2*randi(hStr, [0 1], 200, 1)-1;

NRZ_out=[];
Vp=5;
load('RRC_filtre1.mat');
for index=1:size(x,1)
 if x(index)==1
 NRZ_out=[NRZ_out [1 1 1 1 1 1 1 1 1 1]*Vp];
 elseif x(index)==-1
 NRZ_out=[NRZ_out [1 1 1 1 1 1 1 1 1 1]*(-Vp)];
 end
end


plot(to, NRZ_out, 'r-'); 
title('Kutuplu NRZ transmit data')
% % Set axes and labels.
 axis([0 168 -10 10]);  xlabel('Time (ms)'); ylabel('Amplitude');
 a=upsample(NRZ_out,8);


darbe=conv(a,filtre);
tb = 1000* (0: (1608)*10 - 1) / 10000;
hold on
plot(tb/4,darbe);
title('RRC-filtre1.mat ile convole edildi');
legend('orjinal sinyal','darbe bicimlenmiþ')

% % ---------------------------------------
% % subplot(312)
% % load('RRC_filtre2.mat');
% % 
% % plot(to, NRZ_out, 'r-'); ;
% % title('Kutuplu NRZ transmit data')
% % % Set axes and labels.
% %  axis([0 168 -10 10]);  xlabel('Time (ms)'); ylabel('Amplitude');
% %  a=upsample(NRZ_out,8);
% % darbe=conv(a,filtre);
% % tb = 1000* (0: (1608)*10 - 1) / 10000;
% % hold on
% % plot(tb/4,darbe,'Color','k');
% % title('RRC-filtre2.mat ile convole edildi');
% % legend('orjinal sinyal','darbe bicimlenmiþ')
% 
% % %% --------------------------------------------
% % subplot(313)
% % load('RRC_filtre3.mat');
% % 
% % plot(to, NRZ_out, 'r-'); ;
% % title('Kutuplu NRZ transmit data')
% % % % Set axes and labels.
% %  axis([0 168 -10 10]);  xlabel('Time (ms)'); ylabel('Amplitude');
% %  a=upsample(NRZ_out,8);
% % darbe=conv(a,filtre);
% % tb = 1000* (0: (1608)*10 - 1) / 10000;
% % hold on
% % plot(tb/4,darbe,'Color','k');
% % title('RRC-filtre3.mat ile convole edildi');
% % legend('orjinal sinyal','darbe bicimlenmiþ')
% % % 
% 



