% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This script provides:
% - A way to look into the TPSF 'tpsfs' data generated.
% - Create simulation spectral emission profiles (gaussian for example).
% 
% Jason T. Smith, 10/17/2019, Rensselaer Polytechnic Institute
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% Plot 3D profiles of TPSFs across all channels for visualization
% Specify (x, y)
xX = ; 
yY = ;
figure;
for i = 1:16
    hold on;
    % 3D line-plot
    plot3((1:256),ones([1 256]).*i,squeeze(tpsfs(xX,yY,i,:)),'LineWidth',1.5);
    % 3D fill-plot with transparency
%     fill3((1:256),ones([1 256]).*i,squeeze(tpsfs(xX,yY,i,:)),'b','FaceAlpha',.75);
end

%% Create simulation spectra (gaussians)

fun = @(x,cC)(1./sqrt(2*pi*(cC(1)^2))).*exp(-((x-cC(2)).^2)./(2*cC(1)^2));
s1 = fun((1:16)-6,[1.5 .8]).*.75;
s2 = fun((1:16)-4,[1.25 .7]).*.8;
s3 = fun((1:16)-10,[1.6 .99]).*1.25;

figure; plot(s1); 
hold on;
plot(s2);
hold on;
plot(s3);
xlim([1 16]);
xlabel('channels')
ylabel('intensity (a.u.)')

c1 = [];
c2 = [];
c3 = [];
for i = 1:1000
    c1 = [c1 (s1 + rand(1,16).*.025)'];
    c2 = [c2 (s2 + rand(1,16).*.025)'];
    c3 = [c3 (s3 + rand(1,16).*.015)'];
end
d1 = max(c1.*.5e5-750,0);
d2 = max(c2.*.5e5-750,0);
d3 = max(c3.*.5e5-750,0);

% Visualize
figure; plot((1:16),d1);
hold on; 
plot((1:16),d2);
hold on; 
plot((1:16),d3);
xlim([1 16]);
xlabel('channels')
ylabel('intensity (CW)')

d1N = mean(d1,2)./max(mean(d1,2));
d2N = mean(d2,2)./max(mean(d2,2));
d3N = mean(d3,2)./max(mean(d3,2));

figure; plot((1:16),d1N);
hold on; 
plot((1:16),d2N);
hold on; 
plot((1:16),d3N);
xlim([1 16]);
xlabel('channels')
ylabel('intensity (normalized)')
