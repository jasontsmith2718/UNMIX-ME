% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% * * * * * * RESONANCE ENERGY TRANSFER CORRECTION* * * * * * *
% This script is for generating data for hyperspectral macroscopic
% fluorescence lifetime unmixing. Should be adaptable for other FRET pairs.
% 
% Jason T. Smith, 10/16/2019, Rensselaer Polytechnic Institute
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Load MNIST files (complete set too large for GitHub)
load train_binary
% Load system IRF
load irf
% Load spectra
load FRETpair_spectralD


% How many data samples to simulate?
num = 5000;
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
    r1All = zeros(nX,nY,16); % pre-allocation of individual c1 profiles
    r2All = zeros(nX,nY,16); % pre-allocation of individual c2 profiles
    r3All = zeros(nX,nY,16); % pre-allocation of individual c3 profiles

    %% FRET is THREE-LIFETIME (adapt accordingly)
    % AF700 non-quenched
    t1 = (rand(nX,nY).*.15 + .9).*im;
    % AF700 quenched
    t2 = (rand(nX,nY).*.15 + .25).*im;
    % AF750
    t3 = (rand(nX,nY).*.15 + .55).*im;

    rN = im;
    c_im = zeros([nX nY 3]);
    fretIm = zeros(size(im));
    for i = 1:nX
        for j = 1:nY
            if im(i,j) == 0
                continue
            end
            
            
%% Spectral profiles assigned are dictated by multiples with max-normalized emission spectra
%                 c1 = d1N(:,max(round(rand()*size(d1N,2)),1)).*(rand()*500 + 50);
%                 c2 = c1;
%                 c3 = d3N(:,max(round(rand()*size(d3N,2)),1)).*(rand()*500 + 50);

% Higher noise model
                  c1 = d1N(:,max(round(rand()*size(d1N,2)),1)).*(rand()*225 + 25);
                  c2 = c1;
                  c3 = d3N(:,max(round(rand()*size(d3N,2)),1)).*(rand()*225 + 25);

                  
                rV = rand(); % Which case will it be?
                
                % Mono-exponential (AF700_NQ)
                if rV >= 0 && rV < .2
                    c2 = zeros(size(c2));
                    c3 = zeros(size(c3));
                    r1All(i,j,:) = c1;
                    r2All(i,j,:) = c2;
                    r3All(i,j,:) = c3;
                    c_im(i,j,1) = 1;
                    for k = 1:16
                        e1C = eMon([c1(k) t1(i,j)],time);
                        e1C = conv(e1C,irf./sum(irf));
                        e1C = e1C(1:nTG);
                        e1C = poissrnd(round(e1C'));
                        e2C = eMon([c2(k) t2(i,j)],time);
                        e2C = conv(e2C,irf./sum(irf));
                        e2C = e2C(1:nTG);
                        e2C = poissrnd(round(e2C'));
                        e3C = eMon([c3(k) t3(i,j)],time);
                        e3C = conv(e3C,irf./sum(irf));
                        e3C = e3C(1:nTG);
                        e3C = poissrnd(round(e3C'));
                        eTC = poissrnd(round(e1C'+e2C'+e3C'));
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
                    
                % Mono-exponential (AF750)
                elseif rV >= .2 && rV < .4
                    c1 = zeros(size(c1));
                    c2 = zeros(size(c2));
                    r1All(i,j,:) = c1;
                    r2All(i,j,:) = c2;
                    r3All(i,j,:) = c3;
                    c_im(i,j,3) = 1;
                    for k = 1:16
                        e1C = eMon([c1(k) t1(i,j)],time);
                        e1C = conv(e1C,irf./sum(irf));
                        e1C = e1C(1:nTG);
                        e1C = poissrnd(round(e1C'));
                        e2C = eMon([c2(k) t2(i,j)],time);
                        e2C = conv(e2C,irf./sum(irf));
                        e2C = e2C(1:nTG);
                        e2C = poissrnd(round(e2C'));
                        e3C = eMon([c3(k) t3(i,j)],time);
                        e3C = conv(e3C,irf./sum(irf));
                        e3C = e3C(1:nTG);
                        e3C = poissrnd(round(e3C'));
                        eTC = poissrnd(round(e1C'+e2C'+e3C'));
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
                    
                % Bi-exponential (AF700_NQ, AF750)
                elseif rV >= .4 && rV < .65
                    c2 = zeros(16,1);
                    r1All(i,j,:) = c1;
                    r2All(i,j,:) = c2;
                    r3All(i,j,:) = c3;
                    c_im(i,j,:) = [1;0;1];
                    for k = 1:16
                        e1C = eMon([c1(k) t1(i,j)],time);
                        e1C = conv(e1C,irf./sum(irf));
                        e1C = e1C(1:nTG);
                        e1C = poissrnd(round(e1C'));
                        e2C = eMon([c2(k) t2(i,j)],time);
                        e2C = conv(e2C,irf./sum(irf));
                        e2C = e2C(1:nTG);
                        e2C = poissrnd(round(e2C'));
                        e3C = eMon([c3(k) t3(i,j)],time);
                        e3C = conv(e3C,irf./sum(irf));
                        e3C = e3C(1:nTG);
                        e3C = poissrnd(round(e3C'));
                        eTC = poissrnd(round(e1C'+e2C'+e3C'));
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
                    
                % Tri-exponential (AF700_NQ, AF700_Q, AF750)
                else
                    rFRET = rand()*.7 + .1;
                    r1All(i,j,:) = c1.*(1-rFRET);
                    r2All(i,j,:) = c1.*rFRET;
                    r3All(i,j,:) = c3;
                    fretIm(i,j) = rFRET;
                    c2 = c1.*(rFRET);
                    c1 = c1.*(1-rFRET);
                    c_im(i,j,:) = ones(3,1);
                    for k = 1:16
                        e1C = eMon([c1(k) t1(i,j)],time);
                        e1C = conv(e1C,irf./sum(irf));
                        e1C = e1C(1:nTG);
                        e1C = poissrnd(round(e1C'));
                        e2C = eMon([c2(k) t2(i,j)],time);
                        e2C = conv(e2C,irf./sum(irf));
                        e2C = e2C(1:nTG);
                        e2C = poissrnd(round(e2C'));
                        e3C = eMon([c3(k) t3(i,j)],time);
                        e3C = conv(e3C,irf./sum(irf));
                        e3C = e3C(1:nTG);
                        e3C = poissrnd(round(e3C'));
                        eTC = poissrnd(round(e1C'+e2C'+e3C'));
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
    c = zeros(nX,nY,3); % pre-allocate
    for i = 1:nX
        for j = 1:nY
            curCW=squeeze(cw(i,j,:));
            if sum(curCW) == 0
                continue
            end
            
            if fretIm(i,j) ~=0 % FRET CORRECTION
                fretE = .37; % Determined from literature
                qY_AF700 = .25; % Determined from literature
                qY_AF750 = .12; % Determined from literature
                kL = max(mean(r1,2))./max(mean(r3,2));
                valsC = lsqnonneg([mean(r1N,2).*(1-fretE/1.5) mean(r3N,2)],curCW);
                c(i,j,1) = valsC(1)*(1-fretIm(i,j));
                c(i,j,2) = valsC(1)*fretIm(i,j);
                c(i,j,3) = valsC(2)./(fretE*(qY_AF750/qY_AF700)*kL);
            
            else % ALL OTHER CASES
                valsC = lsqnonneg([mean(r1,2) mean(r3,2)],curCW);
                  c(i,j,1) = valsC(1);
                  c(i,j,3) = valsC(2);
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
    c = c.*c_im;
    t1 = t1.*c_im(:,:,1);
    t2 = t2.*c_im(:,:,2);
    t3 = t3.*c_im(:,:,3);
    fileDir = '';
    save([fileDir 'a_' n], 'c','cw','tpsfs', 't1','t2','t3', '-v7.3')
end