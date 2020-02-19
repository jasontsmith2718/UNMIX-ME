% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This script is for generating data for hyperspectral macroscopic
% fluorescence lifetime unmixing (without RET correction).
% 
% Jason T. Smith, 02/19/2019, Rensselaer Polytechnic Institute
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

load train_binary
load irf
load simulationSpectra

% # time-points for TPSFs
nTG = 256;
% Gate width (32.6e-3 ns)
width = 32.6e-3;
time = [0.5:nTG-0.5]*width;

% # of wavelength channels
channels = 16;
% # of expected fluorescent species
specieN = 4;

% Mono-exponential setup
eMon = @(a,x) a(1).*exp(-x./a(2));

% Target (x,y)
nX = 16;
nY = 16;

% How many data samples to simulate?
num = 500;

for q = 1:num
    
    im = train_images(:,:,max(round(rand()*2000),1)); % Choose random MNIST
    
    if sum(sum(im)) < 250   %1100
        continue
    end
    
    im = imresize(im,[nX nY]); % Resize to (16,16) for memory efficiency
    im = imbinarize(im,.1); % Bring back to binary
    cw = zeros(nX,nY,channels); % pre-allocation of CW spatially
    tpsfs = zeros(nX,nY,channels,nTG); % pre-allocation of TPSF data (4D)
    rAll = zeros(nX,nY,channels,specieN); % pre-allocation of emission profiles (4D)

    t1_1 = im.*(rand()*.15 + 1); % Fluorophore #1_1 (between 1-1.15 ns)
    t1_2 = im.*(rand()*.15 + .45); % Fluorophore #1_2 (between 0.45-0.6 ns)
    t2_1 = im.*(rand()*.15 + .25); % Fluorophore #2_1 (between 0.35-0.5 ns)
    t2_2 = im.*(rand()*.15 + 1.5); % Fluorophore #2_2 (between 1.5-1.65 ns)
     
    % Pre-allocate binary abundance coefficient maps (notating which
    % fluorohores are present at each spatial location)
    cIm = zeros([size(im), specieN]);

    for i = 1:nX
        for j = 1:nY
            if im(i,j) == 0
                continue
            end
            
%% Spectral profiles assigned are dictated by multiples with max-normalized emission spectra
% - A single, experimentally representative emission profile will also do.
% - Each emission profile is multiplied by a number between 75-500 to
% randomize overall contribution.

            % Random emission profile of fluorophore # 1
            c1_1 = d1N(:,max(round(rand()*size(d1N,2)),1)).*(rand()*500 + 75);
            % Random emission profile of fluorophore # 2
            c1_2 = d1N(:,max(round(rand()*size(d1N,2)),1)).*(rand()*500 + 75);
            % Random emission profile of fluorophore # 3
            c2_1 = d2N(:,max(round(rand()*size(d1N,2)),1)).*(rand()*500 + 75);
            % Random emission profile of fluorophore # 4
            c2_2 = d2N(:,max(round(rand()*size(d1N,2)),1)).*(rand()*500 + 75);
            
            % ^^ Concatenated ^^
            cS = [c1_1 c1_2 c2_1 c2_2];
            
            % Random binary assignment of size [specieN, 1]
            cF = round(rand(specieN,1));
            % Assign cF to its spatial location
            cIm(i,j,:) = cF;
            
            % Pre-allocate for TPSF data per fluorophore per channel.
            tpsfC1_1 = zeros(nTG,channels);
            tpsfC1_2 = zeros(nTG,channels);
            tpsfC2_1 = zeros(nTG,channels);
            tpsfC2_2 = zeros(nTG,channels);

            for v = 1:specieN
                % Skip when fluorophores aren't assigned
                if cF(v) == 0
                    continue
                end
                % Fluorophore #1
                if v == 1
                    for ch = 1:channels
                        % Create mono-exponential for fluorophore #1
                        eC = eMon([c1_1(ch) t1_1(i,j)],time);
                        % IRF convolution
                        eC = conv(eC,irf);
                        % Sample back to total # time-points
                        eC = eC(1:nTG);
                        eC = round(eC');
                        % Assign to TPSF array
                        tpsfC1_1(:,ch) = eC;
                    end
                    % Sum over time to get CW for fluorophore #1
                    rAll(i,j,:,v) = sum(tpsfC1_1)';
                % Fluorophore #2
                elseif v == 2
                    for ch = 1:channels
                        % Create mono-exponential for fluorophore #2
                        eC = eMon([c1_2(ch) t1_2(i,j)],time);
                        eC = conv(eC,irf);
                        % Sample back to total # time-points
                        eC = eC(1:nTG);
                        eC = round(eC');
                        % Assign to TPSF array
                        tpsfC1_2(:,ch) = eC;
                    end
                    % Sum over time to get CW for fluorophore #2
                    rAll(i,j,:,v) = sum(tpsfC1_2)';
                % etc., etc., etc...
                elseif v == 3
                    for ch = 1:channels
                        eC = eMon([c2_1(ch) t2_1(i,j)],time);
                        eC = conv(eC,irf);
                        eC = eC(1:nTG);
                        eC = round(eC');
                        tpsfC2_1(:,ch) = eC;
                    end
                    rAll(i,j,:,v) = sum(tpsfC2_1)';
                elseif v == 4
                    for ch = 1:channels
                        eC = eMon([c2_2(ch) t2_2(i,j)],time);
                        eC = conv(eC,irf);
                        eC = eC(1:nTG);
                        eC = round(eC');
                        tpsfC2_2(:,ch) = eC;
                    end
                    rAll(i,j,:,v) = sum(tpsfC2_2)';
                % **Add/subtract 'elseif' statements as needed for target 
                % number of fluorophores **
                end
            end
        % Sum all HFLI data from each fluorophore
        tpsfC = round(tpsfC1_1 + tpsfC1_2 + tpsfC2_1 + tpsfC2_2);
        % Poisson noise
        tpsfs(i,j,:,:) = poissrnd(permute(tpsfC,[2 1]));
        % Assign spatially generated CW over all TPSFs
        cw(i,j,:) = sum(rAll(i,j,:,:),4);
        end
    end
    
    % Pre-allocate for abundance coefficient values
    c = zeros(nX,nY,specieN);
    
    for i = 1:nX
        for j = 1:nY
            curCW=squeeze(cw(i,j,:));
            if sum(curCW) == 0
                continue
            end
                % Nonlinear spectral decomposition
                valsC = lsqnonneg([mean(d1,2) mean(d2,2)],curCW);
                % Obtain fluorophore/spatially-specific CW values
                r1_1 = squeeze(rAll(i,j,:,1));
                r1_2 = squeeze(rAll(i,j,:,2));
                r2_1 = squeeze(rAll(i,j,:,3));
                r2_2 = squeeze(rAll(i,j,:,4));
                
                % Ensure abundance coefficients of fluorophores with same 
                % emission spectra are correctly retrieved.
                c(i,j,1) = valsC(1)*(max(r1_1)/max(r1_1+r1_2));
                c(i,j,2) = valsC(1)*(max(r1_2)/max(r1_1+r1_2));
                c(i,j,3) = valsC(2)*(max(r2_1)/max(r2_1+r2_2));
                c(i,j,4) = valsC(2)*(max(r2_2)/max(r2_1+r2_2));
                
                % Ensure NaN values aren't left around
                if isnan(max(r1_1)/max(r1_1+r1_2)) == 1
                    c(i,j,1:2) = zeros([2 1]);
                elseif isnan(max(r1_2)/max(r1_1+r1_2)) == 1
                    c(i,j,1:2) = zeros([2 1]);
                elseif isnan(max(r2_1)/max(r2_1+r2_2)) == 1
                    c(i,j,3:4) = zeros([2 1]);
                elseif isnan(max(r2_2)/max(r2_1+r2_2)) == 1
                    c(i,j,3:4) = zeros([2 1]);
                end
        end
    end
    
    % Label simulation data
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
    
    % Apply binary mask to abundance coefficients (takes care of spectral
    % bleedthrough)
    c = c.*cIm;
    
    % Create tau matrix and do same
    tau = zeros(nX, nY, specieN);
    tau(:,:,1) = t1_1;
    tau(:,:,2) = t1_2;
    tau(:,:,3) = t2_1;
    tau(:,:,4) = t2_2;
    tau = tau.*cIm;
    
    % Assign directory for saving
    fileDir = '';
    % Save
    save([fileDir '\a_' n], 'c','cw','tpsfs','tau', '-v7.3');
    q = q+1;
end

