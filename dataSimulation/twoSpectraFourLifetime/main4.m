% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This script is for generating data for hyperspectral macroscopic
% fluorescence lifetime unmixing (without RET correction).
% 
% Jason T. Smith, 10/17/2019, Rensselaer Polytechnic Institute
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

load train_binary
load irf
load simulationSpectra

% How many data samples to simulate?
num = 500;
% Time-points
nTG = 256;
% Gate width (32.6e-3 ns)
width = 32.6e-3;
time = [0.5:nTG-0.5]*width;

% Mono-exponential setup
eMon = @(a,x) a(1).*exp(-x./a(2));
% Target (x,y)
nX = 16;
nY = 16;

for q = 1:num
    
    im = train_images(:,:,round(rand()*2000)); % Choose random MNIST
    im = imresize(im,[nX nY]); % Resize to (16,16) for memory efficiency
    im = imbinarize(im,.1); % Bring back to binary
    cw = zeros(nX,nY,16); % pre-allocation of CW spatially
    tpsfs = zeros(nX,nY,16,nTG); % pre-allocation of TPSF data (4D)
    
    r1All_1 = zeros(nX,nY,16); % pre-allocation of individual c1 profiles
    r1All_2 = zeros(nX,nY,16); % pre-allocation of individual c2 profiles
    r2All_1 = zeros(nX,nY,16); % pre-allocation of individual c3 profiles
    r2All_2 = zeros(nX,nY,16); % pre-allocation of individual c3 profiles

    t1_1 = im.*(rand()*.1 + 1); % Fluorophore #1_1
    t1_2 = im.*(rand()*.1 + .45); % Fluorophore #1_2
    t2_1 = im.*(rand()*.1 + .25); % Fluorophore #2_1
    t2_2 = im.*(rand()*.1 + 1.5); % Fluorophore #2_2
    
    c1Im_1 = zeros(size(im)); 
    c1Im_2 = zeros(size(im));
    c2Im_1 = zeros(size(im));
    c2Im_2 = zeros(size(im));

    for i = 1:nX
        for j = 1:nY
            if im(i,j) == 0
                continue
            end
            
%% Spectral profiles assigned are dictated by multiples with max-normalized emission spectra
                c1_1 = d1N(:,max(round(rand()*size(d1N,2)),1)).*(rand()*500 + 75);
                c1_2 = d1N(:,max(round(rand()*size(d1N,2)),1)).*(rand()*500 + 75);
                c2_1 = d2N(:,max(round(rand()*size(d1N,2)),1)).*(rand()*500 + 75);
                c2_2 = d2N(:,max(round(rand()*size(d1N,2)),1)).*(rand()*500 + 75);
                rV = rand();
                % Mono-exponential c1_1
                if rV >= 0 && rV < .1
%                     c1_1 = zeros(size(c1_1));
                    c1_2 = zeros(size(c1_2));
                    c2_1 = zeros(size(c2_1));
                    c2_2 = zeros(size(c2_2));
                    r1All_1(i,j,:) = c1_1;
                    r1All_2(i,j,:) = c1_2;
                    r2All_1(i,j,:) = c2_1;
                    r2All_2(i,j,:) = c2_2;
                    c1Im_1(i,j) = 1;
                    
                    for k = 1:16
                        e1C_1 = eMon([c1_1(k) t1_1(i,j)],time);
                        e1C_1 = conv(e1C_1,irf./sum(irf));
                        e1C_1 = e1C_1(1:nTG);
                        e1C_1 = poissrnd(round(e1C_1'));
                        
                        e1C_2 = eMon([c1_2(k) t1_2(i,j)],time);
                        e1C_2 = conv(e1C_2,irf./sum(irf));
                        e1C_2 = e1C_2(1:nTG);
                        e1C_2 = poissrnd(round(e1C_2'));                       
                        
                        e2C_1 = eMon([c2_1(k) t2_1(i,j)],time);
                        e2C_1 = conv(e2C_1,irf./sum(irf));
                        e2C_1 = e2C_1(1:nTG);
                        e2C_1 = poissrnd(round(e2C_1'));
                        
                        e2C_2 = eMon([c2_2(k) t2_2(i,j)],time);
                        e2C_2 = conv(e2C_2,irf./sum(irf));
                        e2C_2 = e2C_2(1:nTG);
                        e2C_2 = poissrnd(round(e2C_2'));
                       
                        
                        
                        
                        eTC = poissrnd(round(e1C_1'+e1C_2'+e2C_1'+e2C_2'));
                        eTC = eTC + max(rand(256,1).*6 - 3,0);
                        cw(i,j,k) = sum(eTC);
                        rVal = rand();
                        if rVal > .75
                            rG = round(rand()*2 +1);
                            eTC = [zeros(rG,1);eTC];
                        elseif rVal < .25
                            rG = round(rand()*2 +1);
                            eTC = [eTC; zeros(rG,1)];
                        end
                        eTC = eTC(1:256);
                        tpsfs(i,j,k,:) = eTC;
                    end
                
                % Mono-exponential c1_2
                elseif rV >= .1 && rV < .2
                    c1_1 = zeros(size(c1_1));
%                     c1_2 = zeros(size(c1_2));
                    c2_1 = zeros(size(c2_1));
                    c2_2 = zeros(size(c2_2));
                    r1All_1(i,j,:) = c1_1;
                    r1All_2(i,j,:) = c1_2;
                    r2All_1(i,j,:) = c2_1;
                    r2All_2(i,j,:) = c2_2;
                    c1Im_2(i,j) = 1;
                    
                    for k = 1:16
                        e1C_1 = eMon([c1_1(k) t1_1(i,j)],time);
                        e1C_1 = conv(e1C_1,irf./sum(irf));
                        e1C_1 = e1C_1(1:nTG);
                        e1C_1 = poissrnd(round(e1C_1'));
                        
                        e1C_2 = eMon([c1_2(k) t1_2(i,j)],time);
                        e1C_2 = conv(e1C_2,irf./sum(irf));
                        e1C_2 = e1C_2(1:nTG);
                        e1C_2 = poissrnd(round(e1C_2'));                       
                        
                        e2C_1 = eMon([c2_1(k) t2_1(i,j)],time);
                        e2C_1 = conv(e2C_1,irf./sum(irf));
                        e2C_1 = e2C_1(1:nTG);
                        e2C_1 = poissrnd(round(e2C_1'));
                        
                        e2C_2 = eMon([c2_2(k) t2_2(i,j)],time);
                        e2C_2 = conv(e2C_2,irf./sum(irf));
                        e2C_2 = e2C_2(1:nTG);
                        e2C_2 = poissrnd(round(e2C_2'));
                       
                        
                        
                        
                        eTC = poissrnd(round(e1C_1'+e1C_2'+e2C_1'+e2C_2'));
                        eTC = eTC + max(rand(256,1).*6 - 3,0);
                        cw(i,j,k) = sum(eTC);
                        rVal = rand();
                        if rVal > .75
                            rG = round(rand()*2 +1);
                            eTC = [zeros(rG,1);eTC];
                        elseif rVal < .25
                            rG = round(rand()*2 +1);
                            eTC = [eTC; zeros(rG,1)];
                        end
                        eTC = eTC(1:256);
                        tpsfs(i,j,k,:) = eTC;
                    end
                                                       
                % Mono-exponential c2_1
                elseif rV >= .2 && rV < .3
                    c1_1 = zeros(size(c1_1));
                    c1_2 = zeros(size(c1_2));
%                     c2_1 = zeros(size(c2_1));
                    c2_2 = zeros(size(c2_2));
                    r1All_1(i,j,:) = c1_1;
                    r1All_2(i,j,:) = c1_2;
                    r2All_1(i,j,:) = c2_1;
                    r2All_2(i,j,:) = c2_2;
                    c2Im_1(i,j) = 1;
                    
                    for k = 1:16
                        e1C_1 = eMon([c1_1(k) t1_1(i,j)],time);
                        e1C_1 = conv(e1C_1,irf./sum(irf));
                        e1C_1 = e1C_1(1:nTG);
                        e1C_1 = poissrnd(round(e1C_1'));
                        
                        e1C_2 = eMon([c1_2(k) t1_2(i,j)],time);
                        e1C_2 = conv(e1C_2,irf./sum(irf));
                        e1C_2 = e1C_2(1:nTG);
                        e1C_2 = poissrnd(round(e1C_2'));                       
                        
                        e2C_1 = eMon([c2_1(k) t2_1(i,j)],time);
                        e2C_1 = conv(e2C_1,irf./sum(irf));
                        e2C_1 = e2C_1(1:nTG);
                        e2C_1 = poissrnd(round(e2C_1'));
                        
                        e2C_2 = eMon([c2_2(k) t2_2(i,j)],time);
                        e2C_2 = conv(e2C_2,irf./sum(irf));
                        e2C_2 = e2C_2(1:nTG);
                        e2C_2 = poissrnd(round(e2C_2'));
                       
                        
                        
                        
                        eTC = poissrnd(round(e1C_1'+e1C_2'+e2C_1'+e2C_2'));
                        eTC = eTC + max(rand(256,1).*6 - 3,0);
                        cw(i,j,k) = sum(eTC);
                        rVal = rand();
                        if rVal > .75
                            rG = round(rand()*2 +1);
                            eTC = [zeros(rG,1);eTC];
                        elseif rVal < .25
                            rG = round(rand()*2 +1);
                            eTC = [eTC; zeros(rG,1)];
                        end
                        eTC = eTC(1:256);
                        tpsfs(i,j,k,:) = eTC;
                    end
%                     Mono c2_2
                elseif rV >= .3 && rV < .4
                    c1_1 = zeros(size(c1_1));
                    c1_2 = zeros(size(c1_2));
                    c2_1 = zeros(size(c2_1));
%                     c2_2 = zeros(size(c2_2));
                    r1All_1(i,j,:) = c1_1;
                    r1All_2(i,j,:) = c1_2;
                    r2All_1(i,j,:) = c2_1;
                    r2All_2(i,j,:) = c2_2;
                    c2Im_2(i,j) = 1;
                    
                    for k = 1:16
                        e1C_1 = eMon([c1_1(k) t1_1(i,j)],time);
                        e1C_1 = conv(e1C_1,irf./sum(irf));
                        e1C_1 = e1C_1(1:nTG);
                        e1C_1 = poissrnd(round(e1C_1'));
                        
                        e1C_2 = eMon([c1_2(k) t1_2(i,j)],time);
                        e1C_2 = conv(e1C_2,irf./sum(irf));
                        e1C_2 = e1C_2(1:nTG);
                        e1C_2 = poissrnd(round(e1C_2'));                       
                        
                        e2C_1 = eMon([c2_1(k) t2_1(i,j)],time);
                        e2C_1 = conv(e2C_1,irf./sum(irf));
                        e2C_1 = e2C_1(1:nTG);
                        e2C_1 = poissrnd(round(e2C_1'));
                        
                        e2C_2 = eMon([c2_2(k) t2_2(i,j)],time);
                        e2C_2 = conv(e2C_2,irf./sum(irf));
                        e2C_2 = e2C_2(1:nTG);
                        e2C_2 = poissrnd(round(e2C_2'));

                        eTC = poissrnd(round(e1C_1'+e1C_2'+e2C_1'+e2C_2'));
                        eTC = eTC + max(rand(256,1).*6 - 3,0);
                        cw(i,j,k) = sum(eTC);
                        rVal = rand();
                        if rVal > .75
                            rG = round(rand()*2 +1);
                            eTC = [zeros(rG,1);eTC];
                        elseif rVal < .25
                            rG = round(rand()*2 +1);
                            eTC = [eTC; zeros(rG,1)];
                        end
                        eTC = eTC(1:256);
                        tpsfs(i,j,k,:) = eTC;
                    end
                    
%                 bi-exp c1_1 & c1_2
                elseif rV >= .4 && rV < .5
%                     c1_1 = zeros(size(c1_1));
%                     c1_2 = zeros(size(c1_2));
                    c2_1 = zeros(size(c2_1));
                    c2_2 = zeros(size(c2_2));
                    r1All_1(i,j,:) = c1_1;
                    r1All_2(i,j,:) = c1_2;
                    r2All_1(i,j,:) = c2_1;
                    r2All_2(i,j,:) = c2_2;
                    c1Im_1(i,j) = 1;
                    c1Im_2(i,j) = 1;
                    
                    for k = 1:16
                        e1C_1 = eMon([c1_1(k) t1_1(i,j)],time);
                        e1C_1 = conv(e1C_1,irf./sum(irf));
                        e1C_1 = e1C_1(1:nTG);
                        e1C_1 = poissrnd(round(e1C_1'));
                        
                        e1C_2 = eMon([c1_2(k) t1_2(i,j)],time);
                        e1C_2 = conv(e1C_2,irf./sum(irf));
                        e1C_2 = e1C_2(1:nTG);
                        e1C_2 = poissrnd(round(e1C_2'));                       
                        
                        e2C_1 = eMon([c2_1(k) t2_1(i,j)],time);
                        e2C_1 = conv(e2C_1,irf./sum(irf));
                        e2C_1 = e2C_1(1:nTG);
                        e2C_1 = poissrnd(round(e2C_1'));
                        
                        e2C_2 = eMon([c2_2(k) t2_2(i,j)],time);
                        e2C_2 = conv(e2C_2,irf./sum(irf));
                        e2C_2 = e2C_2(1:nTG);
                        e2C_2 = poissrnd(round(e2C_2'));
                       
                        
                        
                        
                        eTC = poissrnd(round(e1C_1'+e1C_2'+e2C_1'+e2C_2'));
                        eTC = eTC + max(rand(256,1).*6 - 3,0);
                        cw(i,j,k) = sum(eTC);
                        rVal = rand();
                        if rVal > .75
                            rG = round(rand()*2 +1);
                            eTC = [zeros(rG,1);eTC];
                        elseif rVal < .25
                            rG = round(rand()*2 +1);
                            eTC = [eTC; zeros(rG,1)];
                        end
                        eTC = eTC(1:256);
                        tpsfs(i,j,k,:) = eTC;
                    end
                % bi-exp c2_1 & c2_2    
                elseif rV >= .5 && rV < .6
                    c1_1 = zeros(size(c1_1));
                    c1_2 = zeros(size(c1_2));
%                     c2_1 = zeros(size(c2_1));
%                     c2_2 = zeros(size(c2_2));
                    r1All_1(i,j,:) = c1_1;
                    r1All_2(i,j,:) = c1_2;
                    r2All_1(i,j,:) = c2_1;
                    r2All_2(i,j,:) = c2_2;
                    c2Im_1(i,j) = 1;
                    c2Im_2(i,j) = 1;
                    
                    for k = 1:16
                        e1C_1 = eMon([c1_1(k) t1_1(i,j)],time);
                        e1C_1 = conv(e1C_1,irf./sum(irf));
                        e1C_1 = e1C_1(1:nTG);
                        e1C_1 = poissrnd(round(e1C_1'));
                        
                        e1C_2 = eMon([c1_2(k) t1_2(i,j)],time);
                        e1C_2 = conv(e1C_2,irf./sum(irf));
                        e1C_2 = e1C_2(1:nTG);
                        e1C_2 = poissrnd(round(e1C_2'));                       
                        
                        e2C_1 = eMon([c2_1(k) t2_1(i,j)],time);
                        e2C_1 = conv(e2C_1,irf./sum(irf));
                        e2C_1 = e2C_1(1:nTG);
                        e2C_1 = poissrnd(round(e2C_1'));
                        
                        e2C_2 = eMon([c2_2(k) t2_2(i,j)],time);
                        e2C_2 = conv(e2C_2,irf./sum(irf));
                        e2C_2 = e2C_2(1:nTG);
                        e2C_2 = poissrnd(round(e2C_2'));
                       
                        
                        
                        
                        eTC = poissrnd(round(e1C_1'+e1C_2'+e2C_1'+e2C_2'));
                        eTC = eTC + max(rand(256,1).*6 - 3,0);
                        cw(i,j,k) = sum(eTC);
                        rVal = rand();
                        if rVal > .75
                            rG = round(rand()*2 +1);
                            eTC = [zeros(rG,1);eTC];
                        elseif rVal < .25
                            rG = round(rand()*2 +1);
                            eTC = [eTC; zeros(rG,1)];
                        end
                        eTC = eTC(1:256);
                        tpsfs(i,j,k,:) = eTC;
                    end
                    % bi-exp c1_2 & c2_1 
                elseif rV >= .6 && rV < .7
                    c1_1 = zeros(size(c1_1));
%                     c1_2 = zeros(size(c1_2));
%                     c2_1 = zeros(size(c2_1));
                    c2_2 = zeros(size(c2_2));
                    r1All_1(i,j,:) = c1_1;
                    r1All_2(i,j,:) = c1_2;
                    r2All_1(i,j,:) = c2_1;
                    r2All_2(i,j,:) = c2_2;
                    c1Im_2(i,j) = 1;
                    c2Im_1(i,j) = 1;
                    
                    for k = 1:16
                        e1C_1 = eMon([c1_1(k) t1_1(i,j)],time);
                        e1C_1 = conv(e1C_1,irf./sum(irf));
                        e1C_1 = e1C_1(1:nTG);
                        e1C_1 = poissrnd(round(e1C_1'));
                        
                        e1C_2 = eMon([c1_2(k) t1_2(i,j)],time);
                        e1C_2 = conv(e1C_2,irf./sum(irf));
                        e1C_2 = e1C_2(1:nTG);
                        e1C_2 = poissrnd(round(e1C_2'));                       
                        
                        e2C_1 = eMon([c2_1(k) t2_1(i,j)],time);
                        e2C_1 = conv(e2C_1,irf./sum(irf));
                        e2C_1 = e2C_1(1:nTG);
                        e2C_1 = poissrnd(round(e2C_1'));
                        
                        e2C_2 = eMon([c2_2(k) t2_2(i,j)],time);
                        e2C_2 = conv(e2C_2,irf./sum(irf));
                        e2C_2 = e2C_2(1:nTG);
                        e2C_2 = poissrnd(round(e2C_2'));
                       
                        
                        
                        
                        eTC = poissrnd(round(e1C_1'+e1C_2'+e2C_1'+e2C_2'));
                        eTC = eTC + max(rand(256,1).*6 - 3,0);
                        cw(i,j,k) = sum(eTC);
                        rVal = rand();
                        if rVal > .75
                            rG = round(rand()*2 +1);
                            eTC = [zeros(rG,1);eTC];
                        elseif rVal < .25
                            rG = round(rand()*2 +1);
                            eTC = [eTC; zeros(rG,1)];
                        end
                        eTC = eTC(1:256);
                        tpsfs(i,j,k,:) = eTC;
                    end
%                 tri-exp c1_1 c1_2 c2_1
                elseif rV >= .7 && rV < .8
%                     c1_1 = zeros(size(c1_1));
%                     c1_2 = zeros(size(c1_2));
%                     c2_1 = zeros(size(c2_1));
                    c2_2 = zeros(size(c2_2));
                    r1All_1(i,j,:) = c1_1;
                    r1All_2(i,j,:) = c1_2;
                    r2All_1(i,j,:) = c2_1;
                    r2All_2(i,j,:) = c2_2;
                    c1Im_1(i,j) = 1;
                    c1Im_2(i,j) = 1;
                    c2Im_1(i,j) = 1;
                    
                    for k = 1:16
                        e1C_1 = eMon([c1_1(k) t1_1(i,j)],time);
                        e1C_1 = conv(e1C_1,irf./sum(irf));
                        e1C_1 = e1C_1(1:nTG);
                        e1C_1 = poissrnd(round(e1C_1'));
                        
                        e1C_2 = eMon([c1_2(k) t1_2(i,j)],time);
                        e1C_2 = conv(e1C_2,irf./sum(irf));
                        e1C_2 = e1C_2(1:nTG);
                        e1C_2 = poissrnd(round(e1C_2'));                       
                        
                        e2C_1 = eMon([c2_1(k) t2_1(i,j)],time);
                        e2C_1 = conv(e2C_1,irf./sum(irf));
                        e2C_1 = e2C_1(1:nTG);
                        e2C_1 = poissrnd(round(e2C_1'));
                        
                        e2C_2 = eMon([c2_2(k) t2_2(i,j)],time);
                        e2C_2 = conv(e2C_2,irf./sum(irf));
                        e2C_2 = e2C_2(1:nTG);
                        e2C_2 = poissrnd(round(e2C_2'));
                       
                        
                        
                        
                        eTC = poissrnd(round(e1C_1'+e1C_2'+e2C_1'+e2C_2'));
                        eTC = eTC + max(rand(256,1).*6 - 3,0);
                        cw(i,j,k) = sum(eTC);
                        rVal = rand();
                        if rVal > .75
                            rG = round(rand()*2 +1);
                            eTC = [zeros(rG,1);eTC];
                        elseif rVal < .25
                            rG = round(rand()*2 +1);
                            eTC = [eTC; zeros(rG,1)];
                        end
                        eTC = eTC(1:256);
                        tpsfs(i,j,k,:) = eTC;
                    end
%                 tri-exp c1_1 c2_1 c2_2
                elseif rV >= .8 && rV < .9
%                     c1_1 = zeros(size(c1_1));
                    c1_2 = zeros(size(c1_2));
%                     c2_1 = zeros(size(c2_1));
%                     c2_2 = zeros(size(c2_2));
                    r1All_1(i,j,:) = c1_1;
                    r1All_2(i,j,:) = c1_2;
                    r2All_1(i,j,:) = c2_1;
                    r2All_2(i,j,:) = c2_2;
                    c1Im_1(i,j) = 1;
                    c2Im_1(i,j) = 1;
                    c2Im_2(i,j) = 1;
                    
                    for k = 1:16
                        e1C_1 = eMon([c1_1(k) t1_1(i,j)],time);
                        e1C_1 = conv(e1C_1,irf./sum(irf));
                        e1C_1 = e1C_1(1:nTG);
                        e1C_1 = poissrnd(round(e1C_1'));
                        
                        e1C_2 = eMon([c1_2(k) t1_2(i,j)],time);
                        e1C_2 = conv(e1C_2,irf./sum(irf));
                        e1C_2 = e1C_2(1:nTG);
                        e1C_2 = poissrnd(round(e1C_2'));                       
                        
                        e2C_1 = eMon([c2_1(k) t2_1(i,j)],time);
                        e2C_1 = conv(e2C_1,irf./sum(irf));
                        e2C_1 = e2C_1(1:nTG);
                        e2C_1 = poissrnd(round(e2C_1'));
                        
                        e2C_2 = eMon([c2_2(k) t2_2(i,j)],time);
                        e2C_2 = conv(e2C_2,irf./sum(irf));
                        e2C_2 = e2C_2(1:nTG);
                        e2C_2 = poissrnd(round(e2C_2'));
                       
                        
                        
                        
                        eTC = poissrnd(round(e1C_1'+e1C_2'+e2C_1'+e2C_2'));
                        eTC = eTC + max(rand(256,1).*6 - 3,0);
                        cw(i,j,k) = sum(eTC);
                        rVal = rand();
                        if rVal > .75
                            rG = round(rand()*2 +1);
                            eTC = [zeros(rG,1);eTC];
                        elseif rVal < .25
                            rG = round(rand()*2 +1);
                            eTC = [eTC; zeros(rG,1)];
                        end
                        eTC = eTC(1:256);
                        tpsfs(i,j,k,:) = eTC;
                    end
%                     ALL
                else
%                     c1_1 = zeros(size(c1_1));
%                     c1_2 = zeros(size(c1_2));
%                     c2_1 = zeros(size(c2_1));
%                     c2_2 = zeros(size(c2_2));
                    r1All_1(i,j,:) = c1_1;
                    r1All_2(i,j,:) = c1_2;
                    r2All_1(i,j,:) = c2_1;
                    r2All_2(i,j,:) = c2_2;
                    c1Im_1(i,j) = 1;
                    c1Im_2(i,j) = 1;
                    c2Im_1(i,j) = 1;
                    c2Im_2(i,j) = 1;
                    
                    for k = 1:16
                        e1C_1 = eMon([c1_1(k) t1_1(i,j)],time);
                        e1C_1 = conv(e1C_1,irf./sum(irf));
                        e1C_1 = e1C_1(1:nTG);
                        e1C_1 = poissrnd(round(e1C_1'));
                        
                        e1C_2 = eMon([c1_2(k) t1_2(i,j)],time);
                        e1C_2 = conv(e1C_2,irf./sum(irf));
                        e1C_2 = e1C_2(1:nTG);
                        e1C_2 = poissrnd(round(e1C_2'));                       
                        
                        e2C_1 = eMon([c2_1(k) t2_1(i,j)],time);
                        e2C_1 = conv(e2C_1,irf./sum(irf));
                        e2C_1 = e2C_1(1:nTG);
                        e2C_1 = poissrnd(round(e2C_1'));
                        
                        e2C_2 = eMon([c2_2(k) t2_2(i,j)],time);
                        e2C_2 = conv(e2C_2,irf./sum(irf));
                        e2C_2 = e2C_2(1:nTG);
                        e2C_2 = poissrnd(round(e2C_2'));
                       
                        
                        
                        
                        eTC = poissrnd(round(e1C_1'+e1C_2'+e2C_1'+e2C_2'));
                        eTC = eTC + max(rand(256,1).*6 - 3,0);
                        cw(i,j,k) = sum(eTC);
                        rVal = rand();
                        if rVal > .75
                            rG = round(rand()*2 +1);
                            eTC = [zeros(rG,1);eTC];
                        elseif rVal < .25
                            rG = round(rand()*2 +1);
                            eTC = [eTC; zeros(rG,1)];
                        end
                        eTC = eTC(1:256);
                        tpsfs(i,j,k,:) = eTC;
                    end                    
                end
        end
    end
    %% Coefficient determination
    c = zeros(nX,nY,4); % pre-allocate
    for i = 1:nX
        for j = 1:nY
            curCW=squeeze(cw(i,j,:));
            if sum(curCW) == 0
                continue
            end
                valsC = lsqnonneg([mean(d1,2) mean(d2,2)],curCW);
                c(i,j,1) = valsC(1)*(max(squeeze(r1All_1(i,j,:)))/max(squeeze(r1All_1(i,j,:))+squeeze(r1All_2(i,j,:))));
                c(i,j,2) = valsC(1)*(max(squeeze(r1All_2(i,j,:)))/max(squeeze(r1All_1(i,j,:))+squeeze(r1All_2(i,j,:))));
                c(i,j,3) = valsC(2)*(max(squeeze(r2All_1(i,j,:)))/max(squeeze(r2All_1(i,j,:))+squeeze(r2All_2(i,j,:))));
                c(i,j,4) = valsC(2)*(max(squeeze(r2All_2(i,j,:)))/max(squeeze(r2All_1(i,j,:))+squeeze(r2All_2(i,j,:))));
                
                if isnan(max(squeeze(r1All_1(i,j,:)))/max(squeeze(r1All_1(i,j,:))+squeeze(r1All_2(i,j,:)))) == 1
                    c(i,j,1:2) = zeros([2 1]);
                elseif isnan(max(squeeze(r1All_2(i,j,:)))/max(squeeze(r1All_1(i,j,:))+squeeze(r1All_2(i,j,:)))) == 1
                    c(i,j,1:2) = zeros([2 1]);
                elseif isnan(max(squeeze(r2All_1(i,j,:)))/max(squeeze(r2All_1(i,j,:))+squeeze(r2All_2(i,j,:)))) == 1
                    c(i,j,3:4) = zeros([2 1]);
                elseif isnan(max(squeeze(r2All_2(i,j,:)))/max(squeeze(r2All_1(i,j,:))+squeeze(r2All_2(i,j,:)))) == 1
                    c(i,j,3:4) = zeros([2 1]);
                end
        end
    end
    if q >=0 && q < 10
        n = ['0000' num2str(q)];
    elseif q >=10 && q<100
        n = ['000' num2str(q)];
    elseif q >=100 && q<1000
        n = ['00' num2str(q)];
    elseif q >=1000 && q<10000
        n = ['0' num2str(q)];
    else
        n = num2str(q);
    end
    c(:,:,1) = c(:,:,1).*c1Im_1;
    c(:,:,2) = c(:,:,2).*c1Im_2;
    c(:,:,3) = c(:,:,3).*c2Im_1;
    c(:,:,4) = c(:,:,4).*c2Im_2;
    
    t1_1 = t1_1.*c1Im_1;
    t1_2 = t1_2.*c1Im_2;
    t2_1 = t2_1.*c2Im_1;
    t2_2 = t2_2.*c2Im_2;
    
    fileDir = '';
    save([fileDir '\a_' n], 'c','cw','tpsfs','t1_1','t1_2','t2_1', 't2_2', '-v7.3')
end

