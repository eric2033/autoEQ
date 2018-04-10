%Zhiguang Eric Zhang - NYU Music Technology - Surgical EQ
%
%primary = priority signal file
%secondary = secondary signal file
%threshold = dB threshold for primary signal for effect
%atten = dB level for primary signal to realize over secondary signal through
%attenuation of secondary signal
%winlength = analysis / synthesis window length (hop size is half)
%equal_loudness = iso226:2003 equal loudness contour scaling (0 or 1)
%
%autoEQnew3stereo('Audio 12_01.wav','Audio 16_01.wav',5,4,2048,1);
%
%function output = autoEQnew3stereo(primary,secondary,threshold,atten,winlength,equal_loudness)
function output = autoEQnew3stereo(primary,secondary,threshold,atten,winlength,equal_loudness)

%hop size
hop = round(winlength/2);

%read signals
[wave1,fs]=audioread(primary);
[wave2,fs]=audioread(secondary);

%check that they are equal length
if size(wave1) ~= size(wave2)
    
    error('Input signals must be the same length and dimension.');
    
end
    
%check mono/stereo
size1 = size(wave1);
size2 = size(wave2);

%allow stereo signals only for this implementation
if size1(2) < 2
    
    error('Stereo signals only.');
    
end

%starting sample
startingSamp = 1;

%create the sine window function
win = sin(pi/winlength*[1:winlength]');

%create a half cosine window
coswin = win(winlength/2+1:end);

%create a half sine window
sinwin = win(1:winlength/2);

%interpolated iso226 equal loudness curve with 40 phon reference (A-weighting)
[~,Vq] = equal_loudness(40,winlength);
[~,Vq_half] = equal_loudness(40,hop);

%number of loops + 1
numblocks = floor(length(wave1)/hop);

%determine how many samples to pad
pad = length(wave1)-(numblocks*hop);

%pad the end of the signals with zeros
wave1 = [wave1;zeros(pad,2)];
wave2 = [wave2;zeros(pad,2)];

%output vectors
output1 = zeros(size(wave1));
output2 = zeros(size(wave2));

%first block
%create a vector containing a block of signal
chunk1aL = wave1(startingSamp : startingSamp + hop - 1, 1);
chunk1aR = wave1(startingSamp : startingSamp + hop - 1, 2);
chunk2aL = wave2(startingSamp : startingSamp + hop - 1, 1);
chunk2aR = wave2(startingSamp : startingSamp + hop - 1, 2);
    
%apply window
winchunk1L = coswin .* chunk1aL;
winchunk1R = coswin .* chunk1aR;
winchunk2L = coswin .* chunk2aL;
winchunk2R = coswin .* chunk2aR;
    
%take fft
WINCHUNK1L = fft(winchunk1L);
WINCHUNK1R = fft(winchunk1R);
WINCHUNK2L = fft(winchunk2L);
WINCHUNK2R = fft(winchunk2R);
    
%auto eq algorithm
for k = 1:length(WINCHUNK1L)
        
        %if primary freq bin magnitude is less than XdB
        %if abs(WINCHUNK1(i)) <= db2mag(atten)
            
            %set secondary freq bin to zero
        %    WINCHUNK2(i) = 0;
        
        %if primary freq L bin dB is greater than threshold
        if mag2db(abs(WINCHUNK1L(k))) > threshold
            
            %get phases
            PHASE1L(k) = angle(WINCHUNK1L(k));
            PHASE2L(k) = angle(WINCHUNK2L(k));
            
            %get magnitudes
            MAG1L(k) = abs(WINCHUNK1L(k));
            MAG2L(k) = abs(WINCHUNK2L(k));
            
            %adjust phase to be from 0 to 2 * pi
            if PHASE1L(k) < 0
                
                PHASE1xL(k) = PHASE1L(k) + (2*pi);
                
            else
                
                PHASE1xL(k) = PHASE1L(k);
                
            end
            
            %adjust phase to be from 0 to 2 * pi
            if PHASE2L(k) < 0
                
                PHASE2xL(k) = PHASE2L(k) + (2*pi);
                
            else
                
                PHASE2xL(k) = PHASE2L(k);
                
            end
            
            %if vectors point in exact opposite directions
            if abs(PHASE1xL(k) - PHASE2xL(k)) == pi || PHASE1xL(k) == PHASE2xL(k)
                
                %if primary freq bin is less than XdB above secondary
                if mag2db(MAG1L(k)) - mag2db(MAG2L(k)) < atten 
                    
                    %attenuate secondary freq bin such that primary is XdB
                    %above secondary
                    
                    if (equal_loudness)
                    
                        MAG2L(k) = MAG2L(k) - (MAG2L(k) - db2mag(mag2db(MAG1L(k)) - (atten - (atten*Vq_half(k)))));
                        
                    else
                        
                        MAG2L(k) = MAG2L(k) - (MAG2L(k) - db2mag(mag2db(MAG1L(k)) - atten));
                    
                    end
                        
                    %polar to cartesian
                    [re, im] = pol2cart(PHASE2L(k), MAG2L(k));
            
                    %rectangular to complex
                    WINCHUNK2L(k) = re + 1i*im;
                
                end    
            
            %otherwise
            else
                
                %find resultant magnitude
                MAGNITUDEL(k) = sqrt((real(WINCHUNK1L(k))+real(WINCHUNK2L(k)))^2 + (imag(WINCHUNK1L(k))+imag(WINCHUNK2L(k)))^2);
                
                %find resultant angle
                RESULTANTL(k) = atan2((imag(WINCHUNK1L(k))+imag(WINCHUNK2L(k))),(real(WINCHUNK1L(k))+real(WINCHUNK2L(k))));
                
                %adjust the phase of resultant angle to be from 0 to 2 * pi
                if RESULTANTL(k) < 0
                
                    RESULTANTxL(k) = RESULTANTL(k) + (2*pi);
                
                else
                
                    RESULTANTxL(k) = RESULTANTL(k);
                
                end
                
                %calculate primary theta
                THETA1L(k) = abs(PHASE1xL(k) - RESULTANTxL(k));
                
                %calculate secondary theta
                THETA2L(k) = abs(PHASE2xL(k) - RESULTANTxL(k));
                
                %calculate ratio of hypotenuse to missing side for primary
                RATIO1L(k) = abs(sin(THETA1L(k)));
                
                %calculate ratio of hypotenuse to missing side for
                %secondary
                RATIO2L(k) = abs(sin(THETA2L(k)));
                
                %calculate primary missing side
                DROP1L(k) = RATIO1L(k) * MAG1L(k);
                
                %calculate secondary missing side
                DROP2L(k) = RATIO2L(k) * MAG2L(k);
                               
                %if primary contribution is less than XdB above secondary
                %contribution
                if mag2db(sqrt(MAG1L(k)^2 - DROP1L(k)^2)) - mag2db(sqrt(MAG2L(k)^2 - DROP2L(k)^2)) < atten
                    
                    %calculate the difference in contribution necessary to
                    %achieve XdB gain of primary over secondary
                    
                    if (equal_loudness)
                        
                        DIFFL(k) = sqrt(MAG2L(k)^2 - DROP2L(k)^2) - db2mag(mag2db(sqrt(MAG1L(k)^2 - DROP1L(k)^2)) - (atten - (atten*Vq_half(k))));
                    
                    else
                        
                        DIFFL(k) = sqrt(MAG2L(k)^2 - DROP2L(k)^2) - db2mag(mag2db(sqrt(MAG1L(k)^2 - DROP1L(k)^2)) - atten);
                        
                    end
                        
                    %calculate ratio of drop to difference in
                    %contribution for primary
                    DELTARATIO1L(k) = abs(tan(THETA1L(k)));
                    
                    %calculate ratio of drop to difference in contribution
                    %for secondary
                    DELTARATIO2L(k) = abs(tan(THETA2L(k)));
                    
                    %if primary contribution is less than XdB
                    %if mag2db(sqrt(MAG1(k)^2 - DROP1(k)^2)) < atten
                    
                        %boost primary
                    %    MAG1x(k) = MAG1(k) + sqrt((DELTARATIO1(k) * DIFF(k))^2 + DIFF(k)^2);
                        
                        %polar to cartesian
                    %    [re, im] = pol2cart(PHASE1(k), MAG1x(k));
            
                        %rectangular to complex
                    %    WINCHUNK1(k) = re + i*im;
                        
                    %otherwise
                    %else
                        
                        %attenuate secondary
                        MAG2xL(k) = MAG2L(k) - sqrt((DELTARATIO2L(k) * DIFFL(k))^2 + DIFFL(k)^2);
                        
                        %polar to cartesian
                        [re, im] = pol2cart(PHASE2L(k), MAG2xL(k));
            
                        %rectangular to complex
                        WINCHUNK2L(k) = re + 1i*im;
                        
                    %end
                    
                end
                
            end
            
        end
        
        %if primary freq R bin dB is greater than threshold
        if mag2db(abs(WINCHUNK1R(k))) > threshold
            
            %get phases
            PHASE1R(k) = angle(WINCHUNK1R(k));
            PHASE2R(k) = angle(WINCHUNK2R(k));
        
            %get magnitudes
            MAG1R(k) = abs(WINCHUNK1R(k));
            MAG2R(k) = abs(WINCHUNK2R(k));
            
            %adjust phase to be from 0 to 2 * pi
            if PHASE1R(k) < 0
                
                PHASE1xR(k) = PHASE1R(k) + (2*pi);
                
            else
                
                PHASE1xR(k) = PHASE1R(k);
                
            end
            
            %adjust phase to be from 0 to 2 * pi
            if PHASE2R(k) < 0
                
                PHASE2xR(k) = PHASE2R(k) + (2*pi);
                
            else
                
                PHASE2xR(k) = PHASE2R(k);
                
            end
            
            %if vectors point in exact opposite directions or same
            %direction
            if abs(PHASE1xR(k) - PHASE2xR(k)) == pi || PHASE1xR(k) == PHASE2xR(k)
                
                %if primary freq bin is less than XdB above secondary
                if mag2db(MAG1R(k)) - mag2db(MAG2R(k)) < atten 
                    
                    %attenuate secondary freq bin such that primary is XdB
                    %above secondary
                    
                    if (equal_loudness)
                        
                        MAG2R(k) = MAG2R(k) - (MAG2R(k) - db2mag(mag2db(MAG1R(k)) - (atten - (atten*Vq_half(k)))));
                        
                    else
                        
                        MAG2R(k) = MAG2R(k) - (MAG2R(k) - db2mag(mag2db(MAG1R(k)) - atten));
                        
                    end
                    
                    %polar to cartesian
                    [re, im] = pol2cart(PHASE2R(k), MAG2R(k));
            
                    %rectangular to complex
                    WINCHUNK2R(k) = re + 1i*im;
                
                end    
            
            %otherwise
            else
                
                %find resultant magnitude
                MAGNITUDER(k) = sqrt((real(WINCHUNK1R(k))+real(WINCHUNK2R(k)))^2 + (imag(WINCHUNK1R(k))+imag(WINCHUNK2R(k)))^2);
                
                %find resultant angle
                RESULTANTR(k) = atan2((imag(WINCHUNK1R(k))+imag(WINCHUNK2R(k))),(real(WINCHUNK1R(k))+real(WINCHUNK2R(k))));
                
                %adjust the phase of resultant angle to be from 0 to 2 * pi
                if RESULTANTR(k) < 0
                
                    RESULTANTxR(k) = RESULTANTR(k) + (2*pi);
                
                else
                
                    RESULTANTxR(k) = RESULTANTR(k);
                
                end
                
                %calculate primary theta
                THETA1R(k) = abs(PHASE1xR(k) - RESULTANTxR(k));
                
                %calculate secondary theta
                THETA2R(k) = abs(PHASE2xR(k) - RESULTANTxR(k));
                
                %calculate ratio of hypotenuse to missing side for primary
                RATIO1R(k) = abs(sin(THETA1R(k)));
                
                %calculate ratio of hypotenuse to missing side for
                %secondary
                RATIO2R(k) = abs(sin(THETA2R(k)));
                
                %calculate primary missing side
                DROP1R(k) = RATIO1R(k) * MAG1R(k);
                
                %calculate secondary missing side
                DROP2R(k) = RATIO2R(k) * MAG2R(k);
                               
                %if primary contribution is less than XdB above secondary
                %contribution
                if mag2db(sqrt(MAG1R(k)^2 - DROP1R(k)^2)) - mag2db(sqrt(MAG2R(k)^2 - DROP2R(k)^2)) < atten
                    
                    %calculate the difference in contribution necessary to
                    %achieve XdB gain of primary over secondary
                    
                    if (equal_loudness)
                        
                        DIFFR(k) = sqrt(MAG2R(k)^2 - DROP2R(k)^2) - db2mag(mag2db(sqrt(MAG1R(k)^2 - DROP1R(k)^2)) - (atten - (atten*Vq_half(k))));
                        
                    else
                        
                        DIFFR(k) = sqrt(MAG2R(k)^2 - DROP2R(k)^2) - db2mag(mag2db(sqrt(MAG1R(k)^2 - DROP1R(k)^2)) - atten);
                        
                    end
                    
                    %calculate ratio of drop to difference in
                    %contribution for primary
                    DELTARATIO1R(k) = abs(tan(THETA1R(k)));
                    
                    %calculate ratio of drop to difference in contribution
                    %for secondary
                    DELTARATIO2R(k) = abs(tan(THETA2R(k)));
                    
                    %if primary contribution is less than XdB
                    %if mag2db(sqrt(MAG1(k)^2 - DROP1(k)^2)) < atten
                    
                        %boost primary
                    %    MAG1x(k) = MAG1(k) + sqrt((DELTARATIO1(k) * DIFF(k))^2 + DIFF(k)^2);
                        
                        %polar to cartesian
                    %    [re, im] = pol2cart(PHASE1(k), MAG1x(k));
            
                        %rectangular to complex
                    %    WINCHUNK1(k) = re + i*im;
                        
                    %otherwise
                    %else
                        
                        %attenuate secondary
                        MAG2xR(k) = MAG2R(k) - sqrt((DELTARATIO2R(k) * DIFFR(k))^2 + DIFFR(k)^2);
                        
                        %polar to cartesian
                        [re, im] = pol2cart(PHASE2R(k), MAG2xR(k));
            
                        %rectangular to complex
                        WINCHUNK2R(k) = re + 1i*im;
                        
                    %end
                    
                end
                
            end
            
        end
          
    end
    
%take inverse fft
winchunk1L = real(ifft(WINCHUNK1L));
winchunk1R = real(ifft(WINCHUNK1R));
winchunk2L = real(ifft(WINCHUNK2L));
winchunk2R = real(ifft(WINCHUNK2R));
    
%apply window
winchunk1L = coswin .* winchunk1L;
winchunk1R = coswin .* winchunk1R;
winchunk2L = coswin .* winchunk2L;
winchunk2R = coswin .* winchunk2R;
    
%overlap and add
output1(startingSamp:startingSamp + hop - 1, 1) = output1(startingSamp:startingSamp + hop - 1, 1) + winchunk1L;
output1(startingSamp:startingSamp + hop - 1, 2) = output1(startingSamp:startingSamp + hop - 1, 2) + winchunk1R;
output2(startingSamp:startingSamp + hop - 1, 1) = output2(startingSamp:startingSamp + hop - 1, 1) + winchunk2L;
output2(startingSamp:startingSamp + hop - 1, 2) = output2(startingSamp:startingSamp + hop - 1, 2) + winchunk2R;

%iterate over blocks
for j = 1:numblocks - 1

    %create a vector containing a block of signal
    chunk1aL = wave1(startingSamp : startingSamp + hop - 1, 1);
    chunk1aR = wave1(startingSamp : startingSamp + hop - 1, 2);
    chunk2aL = wave2(startingSamp : startingSamp + hop - 1, 1);
    chunk2aR = wave2(startingSamp : startingSamp + hop - 1, 2);
    
    %retain the starting sample position
    startingSamp1 = startingSamp;
        
    %set the next starting sample position
    startingSamp = (startingSamp + winlength) - hop;
    
    %segment second block
    chunk1bL = wave1(startingSamp : startingSamp + hop - 1, 1);
    chunk1bR = wave1(startingSamp : startingSamp + hop - 1, 2);
    chunk2bL = wave2(startingSamp : startingSamp + hop - 1, 1);
    chunk2bR = wave2(startingSamp : startingSamp + hop - 1, 2);

    %concatenate chunks
    chunk1L = [chunk1aL;chunk1bL];
    chunk1R = [chunk1aR;chunk1bR];
    chunk2L = [chunk2aL;chunk2bL];
    chunk2R = [chunk2aR;chunk2bR];
    
    %apply window
    winchunk1L = win .* chunk1L;
    winchunk1R = win .* chunk1R;
    winchunk2L = win .* chunk2L;
    winchunk2R = win .* chunk2R;
    
    %take fft
    WINCHUNK1L = fft(winchunk1L);
    WINCHUNK1R = fft(winchunk1R);
    WINCHUNK2L = fft(winchunk2L);
    WINCHUNK2R = fft(winchunk2R);
    
    %auto eq algorithm
    for k = 1:length(WINCHUNK1L)
        
        %if primary freq bin magnitude is less than XdB
        %if abs(WINCHUNK1(i)) <= db2mag(atten)
            
            %set secondary freq bin to zero
        %    WINCHUNK2(i) = 0;
        
        %if primary freq bin L dB is greater than threshold
        if mag2db(abs(WINCHUNK1L(k))) > threshold
            
            %get phases
            PHASE1L(k) = angle(WINCHUNK1L(k));
            PHASE2L(k) = angle(WINCHUNK2L(k));
                
            %get magnitudes
            MAG1L(k) = abs(WINCHUNK1L(k));
            MAG2L(k) = abs(WINCHUNK2L(k));
            
            %adjust phase to be from 0 to 2 * pi
            if PHASE1L(k) < 0
                
                PHASE1xL(k) = PHASE1L(k) + (2*pi);
                
            else
                
                PHASE1xL(k) = PHASE1L(k);
                
            end
            
            %adjust phase to be from 0 to 2 * pi
            if PHASE2L(k) < 0
                
                PHASE2xL(k) = PHASE2L(k) + (2*pi);
                
            else
                
                PHASE2xL(k) = PHASE2L(k);
                
            end
            
            %if vectors point in exact opposite directions or same
            %direction
            if abs(PHASE1xL(k) - PHASE2xL(k)) == pi || PHASE1xL(k) == PHASE2xL(k)
                
                %if primary freq bin is less than XdB above secondary
                if mag2db(MAG1L(k)) - mag2db(MAG2L(k)) < atten 
                    
                    %attenuate secondary freq bin such that primary is XdB
                    %above secondary
                    
                    if (equal_loudness)
                        
                        MAG2L(k) = MAG2L(k) - (MAG2L(k) - db2mag(mag2db(MAG1L(k)) - (atten - (atten*Vq(k)))));
                        
                    else
                        
                        MAG2L(k) = MAG2L(k) - (MAG2L(k) - db2mag(mag2db(MAG1L(k)) - atten));
                        
                    end
                    
                    %polar to cartesian
                    [re, im] = pol2cart(PHASE2L(k), MAG2L(k));
            
                    %rectangular to complex
                    WINCHUNK2L(k) = re + 1i*im;
                
                end    
            
            %otherwise
            else
                
                %find resultant magnitude
                MAGNITUDEL(k) = sqrt((real(WINCHUNK1L(k))+real(WINCHUNK2L(k)))^2 + (imag(WINCHUNK1L(k))+imag(WINCHUNK2L(k)))^2);
                
                %find resultant angle
                RESULTANTL(k) = atan2((imag(WINCHUNK1L(k))+imag(WINCHUNK2L(k))),(real(WINCHUNK1L(k))+real(WINCHUNK2L(k))));
                
                %adjust the phase of resultant angle to be from 0 to 2 * pi
                if RESULTANTL(k) < 0
                
                    RESULTANTxL(k) = RESULTANTL(k) + (2*pi);
                
                else
                
                    RESULTANTxL(k) = RESULTANTL(k);
                
                end
                
                %calculate primary theta
                THETA1L(k) = abs(PHASE1xL(k) - RESULTANTxL(k));
                
                %calculate secondary theta
                THETA2L(k) = abs(PHASE2xL(k) - RESULTANTxL(k));
                
                %calculate ratio of hypotenuse to missing side for primary
                RATIO1L(k) = abs(sin(THETA1L(k)));
                
                %calculate ratio of hypotenuse to missing side for
                %secondary
                RATIO2L(k) = abs(sin(THETA2L(k)));
                
                %calculate primary missing side
                DROP1L(k) = RATIO1L(k) * MAG1L(k);
                
                %calculate secondary missing side
                DROP2L(k) = RATIO2L(k) * MAG2L(k);
                               
                %if primary contribution is less than XdB above secondary
                %contribution
                if mag2db(sqrt(MAG1L(k)^2 - DROP1L(k)^2)) - mag2db(sqrt(MAG2L(k)^2 - DROP2L(k)^2)) < atten
                    
                    %calculate the difference in contribution necessary to
                    %achieve XdB gain of primary over secondary
                    
                    if (equal_loudness)
                        
                        DIFFL(k) = sqrt(MAG2L(k)^2 - DROP2L(k)^2) - db2mag(mag2db(sqrt(MAG1L(k)^2 - DROP1L(k)^2)) - (atten - (atten*Vq(k))));
                        
                    else
                        
                        DIFFL(k) = sqrt(MAG2L(k)^2 - DROP2L(k)^2) - db2mag(mag2db(sqrt(MAG1L(k)^2 - DROP1L(k)^2)) - atten);
                        
                    end
                    
                    %calculate ratio of drop to difference in
                    %contribution for primary
                    DELTARATIO1L(k) = abs(tan(THETA1L(k)));
                   
                    %calculate ratio of drop to difference in contribution
                    %for secondary
                    DELTARATIO2L(k) = abs(tan(THETA2L(k)));
                    
                    %if primary contribution is less than XdB (this is an
                    %erroneous conditional as dB is relative)
                    %if mag2db(sqrt(MAG1(k)^2 - DROP1(k)^2)) < atten
                    
                        %boost primary (hypotenuse)
                        %MAG1x(k) = MAG1(k) + sqrt((DELTARATIO1(k) * DIFF(k))^2 + DIFF(k)^2);
                        
                        %polar to cartesian
                        %[re, im] = pol2cart(PHASE1(k), MAG1x(k));
            
                        %rectangular to complex
                        %WINCHUNK1(k) = re + i*im;
                        
                    %otherwise
                    %else
                        
                        %attenuate secondary
                        MAG2xL(k) = MAG2L(k) - sqrt((DELTARATIO2L(k) * DIFFL(k))^2 + DIFFL(k)^2);
                        
                        %polar to cartesian
                        [re, im] = pol2cart(PHASE2L(k), MAG2xL(k));
            
                        %rectangular to complex
                        WINCHUNK2L(k) = re + 1i*im;
                        
                    %end
                    
                end
                
            end
            
        end
        
        %if primary freq bin L dB is greater than threshold
        if mag2db(abs(WINCHUNK1R(k))) > threshold
            
            %get phases
            PHASE1R(k) = angle(WINCHUNK1R(k));
            PHASE2R(k) = angle(WINCHUNK2R(k));
        
            %get magnitudes
            MAG1R(k) = abs(WINCHUNK1R(k));
            MAG2R(k) = abs(WINCHUNK2R(k));
            
            %adjust phase to be from 0 to 2 * pi
            if PHASE1R(k) < 0
                
                PHASE1xR(k) = PHASE1R(k) + (2*pi);
                
            else
                
                PHASE1xR(k) = PHASE1R(k);
                
            end
            
            %adjust phase to be from 0 to 2 * pi
            if PHASE2R(k) < 0
                
                PHASE2xR(k) = PHASE2R(k) + (2*pi);
                
            else
                
                PHASE2xR(k) = PHASE2R(k);
                
            end
            
            %if vectors point in exact opposite directions or same
            %direction
            if abs(PHASE1xR(k) - PHASE2xR(k)) == pi || PHASE1xR(k) == PHASE2xR(k)
                
                %if primary freq bin is less than XdB above secondary
                if mag2db(MAG1R(k)) - mag2db(MAG2R(k)) < atten 
                    
                    %attenuate secondary freq bin such that primary is XdB
                    %above secondary
                    
                    if (equal_loudness)
                        
                        MAG2R(k) = MAG2R(k) - (MAG2R(k) - db2mag(mag2db(MAG1R(k)) - (atten - (atten*Vq(k)))));
                    
                    else
                        
                        MAG2R(k) = MAG2R(k) - (MAG2R(k) - db2mag(mag2db(MAG1R(k)) - atten));
                        
                    end
                    
                    %polar to cartesian
                    [re, im] = pol2cart(PHASE2R(k), MAG2R(k));
            
                    %rectangular to complex
                    WINCHUNK2R(k) = re + 1i*im;
                
                end    
            
            %otherwise
            else
                
                %find resultant magnitude
                MAGNITUDER(k) = sqrt((real(WINCHUNK1R(k))+real(WINCHUNK2R(k)))^2 + (imag(WINCHUNK1R(k))+imag(WINCHUNK2R(k)))^2);
                
                %find resultant angle
                RESULTANTR(k) = atan2((imag(WINCHUNK1R(k))+imag(WINCHUNK2R(k))),(real(WINCHUNK1R(k))+real(WINCHUNK2R(k))));
                
                %adjust the phase of resultant angle to be from 0 to 2 * pi
                if RESULTANTR(k) < 0
                
                    RESULTANTxR(k) = RESULTANTR(k) + (2*pi);
                
                else
                
                    RESULTANTxR(k) = RESULTANTR(k);
                
                end
                
                %calculate primary theta
                THETA1R(k) = abs(PHASE1xR(k) - RESULTANTxR(k));
                
                %calculate secondary theta
                THETA2R(k) = abs(PHASE2xR(k) - RESULTANTxR(k));
                
                %calculate ratio of hypotenuse to missing side for primary
                RATIO1R(k) = abs(sin(THETA1R(k)));
                
                %calculate ratio of hypotenuse to missing side for
                %secondary
                RATIO2R(k) = abs(sin(THETA2R(k)));
                
                %calculate primary missing side
                DROP1R(k) = RATIO1R(k) * MAG1R(k);
                
                %calculate secondary missing side
                DROP2R(k) = RATIO2R(k) * MAG2R(k);
                               
                %if primary contribution is less than XdB above secondary
                %contribution
                if mag2db(sqrt(MAG1R(k)^2 - DROP1R(k)^2)) - mag2db(sqrt(MAG2R(k)^2 - DROP2R(k)^2)) < atten
                    
                    %calculate the difference in contribution necessary to
                    %achieve XdB gain of primary over secondary
                    
                    if (equal_loudness)
                                           
                        DIFFR(k) = sqrt(MAG2R(k)^2 - DROP2R(k)^2) - db2mag(mag2db(sqrt(MAG1R(k)^2 - DROP1R(k)^2)) - (atten - (atten*Vq(k))));
                        
                    else
                                                
                        DIFFR(k) = sqrt(MAG2R(k)^2 - DROP2R(k)^2) - db2mag(mag2db(sqrt(MAG1R(k)^2 - DROP1R(k)^2)) - atten);
                        
                    end
                    
                    %calculate ratio of drop to difference in
                    %contribution for primary
                    DELTARATIO1R(k) = abs(tan(THETA1R(k)));
                   
                    %calculate ratio of drop to difference in contribution
                    %for secondary
                    DELTARATIO2R(k) = abs(tan(THETA2R(k)));
                    
                    %if primary contribution is less than XdB (this is an
                    %erroneous conditional as dB is relative)
                    %if mag2db(sqrt(MAG1(k)^2 - DROP1(k)^2)) < atten
                    
                        %boost primary (hypotenuse)
                        %MAG1x(k) = MAG1(k) + sqrt((DELTARATIO1(k) * DIFF(k))^2 + DIFF(k)^2);
                        
                        %polar to cartesian
                        %[re, im] = pol2cart(PHASE1(k), MAG1x(k));
            
                        %rectangular to complex
                        %WINCHUNK1(k) = re + i*im;
                        
                    %otherwise
                    %else
                        
                        %attenuate secondary
                        MAG2xR(k) = MAG2R(k) - sqrt((DELTARATIO2R(k) * DIFFR(k))^2 + DIFFR(k)^2);
                        
                        %polar to cartesian
                        [re, im] = pol2cart(PHASE2R(k), MAG2xR(k));
            
                        %rectangular to complex
                        WINCHUNK2R(k) = re + 1i*im;
                        
                    %end
                    
                end
                
            end
            
        end
          
    end
    
    %take inverse fft
    winchunk1L = real(ifft(WINCHUNK1L));
    winchunk1R = real(ifft(WINCHUNK1R));
    winchunk2L = real(ifft(WINCHUNK2L));
    winchunk2R = real(ifft(WINCHUNK2R));
    
    %apply window
    winchunk1L = win .* winchunk1L;
    winchunk1R = win .* winchunk1R;
    winchunk2L = win .* winchunk2L;
    winchunk2R = win .* winchunk2R;
    
    %overlap and add
    output1(startingSamp1:startingSamp + hop - 1, 1) = output1(startingSamp1:startingSamp + hop - 1, 1) + winchunk1L;
    output1(startingSamp1:startingSamp + hop - 1, 2) = output1(startingSamp1:startingSamp + hop - 1, 2) + winchunk1R;
    output2(startingSamp1:startingSamp + hop - 1, 1) = output2(startingSamp1:startingSamp + hop - 1, 1) + winchunk2L;
    output2(startingSamp1:startingSamp + hop - 1, 2) = output2(startingSamp1:startingSamp + hop - 1, 2) + winchunk2R;

end

%final block
%create a vector containing a block of audio
chunk1aL = wave1(startingSamp : startingSamp + hop - 1, 1);
chunk1aR = wave1(startingSamp : startingSamp + hop - 1, 2);
chunk2aL = wave2(startingSamp : startingSamp + hop - 1, 1);
chunk2aR = wave2(startingSamp : startingSamp + hop - 1, 2);
    
%retain the starting sample position
startingSamp1 = startingSamp;
        
%set the next starting sample position
startingSamp = (startingSamp + winlength) - hop;
    
%apply window
winchunk1L = sinwin .* chunk1aL;
winchunk1R = sinwin .* chunk1aR;
winchunk2L = sinwin .* chunk2aL;
winchunk2R = sinwin .* chunk2aR;
    
%take fft
WINCHUNK1L = fft(winchunk1L);
WINCHUNK1R = fft(winchunk1R);
WINCHUNK2L = fft(winchunk2L);
WINCHUNK2R = fft(winchunk2R);
   
%auto eq algorithm
for k = 1:length(WINCHUNK1L)
        
        %if primary freq bin magnitude is less than XdB
        %if abs(WINCHUNK1(i)) <= db2mag(atten)
            
            %set secondary freq bin to zero
        %    WINCHUNK2(i) = 0;
        
        %if primary freq bin dB is greater than threshold
        if mag2db(abs(WINCHUNK1L(k))) > threshold
            
            %get phases
            PHASE1L(k) = angle(WINCHUNK1L(k));
            PHASE2L(k) = angle(WINCHUNK2L(k));
        
            %get magnitudes
            MAG1L(k) = abs(WINCHUNK1L(k));
            MAG2L(k) = abs(WINCHUNK2L(k));
            
            %adjust phase to be from 0 to 2 * pi
            if PHASE1L(k) < 0
                
                PHASE1xL(k) = PHASE1L(k) + (2*pi);
                
            else
                
                PHASE1xL(k) = PHASE1L(k);
                
            end
            
            %adjust phase to be from 0 to 2 * pi
            if PHASE2L(k) < 0
                
                PHASE2xL(k) = PHASE2L(k) + (2*pi);
                
            else
                
                PHASE2xL(k) = PHASE2L(k);
                
            end
            
            %if vectors point in exact opposite directions or same
            %direction
            if abs(PHASE1xL(k) - PHASE2xL(k)) == pi || PHASE1xL(k) == PHASE2xL(k)
                
                %if primary freq bin is less than XdB above secondary
                if mag2db(MAG1L(k)) - mag2db(MAG2L(k)) < atten 
                    
                    %attenuate secondary freq bin such that primary is XdB
                    %above secondary
                    
                    if (equal_loudness)
                    
                        MAG2L(k) = MAG2L(k) - (MAG2L(k) - db2mag(mag2db(MAG1L(k)) - (atten - (atten*Vq_half(k)))));
                        
                    else
                        
                        MAG2L(k) = MAG2L(k) - (MAG2L(k) - db2mag(mag2db(MAG1L(k)) - atten));
                        
                    end
                    
                    %polar to cartesian
                    [re, im] = pol2cart(PHASE2L(k), MAG2L(k));
            
                    %rectangular to complex
                    WINCHUNK2L(k) = re + 1i*im;
                
                end    
            
            %otherwise
            else
                
                %find resultant magnitude
                MAGNITUDEL(k) = sqrt((real(WINCHUNK1L(k))+real(WINCHUNK2L(k)))^2 + (imag(WINCHUNK1L(k))+imag(WINCHUNK2L(k)))^2);
                
                %find resultant angle
                RESULTANTL(k) = atan2((imag(WINCHUNK1L(k))+imag(WINCHUNK2L(k))),(real(WINCHUNK1L(k))+real(WINCHUNK2L(k))));
                
                %adjust the phase of resultant angle to be from 0 to 2 * pi
                if RESULTANTL(k) < 0
                
                    RESULTANTxL(k) = RESULTANTL(k) + (2*pi);
                
                else
                
                    RESULTANTxL(k) = RESULTANTL(k);
                
                end
                
                %calculate primary theta
                THETA1L(k) = abs(PHASE1xL(k) - RESULTANTxL(k));
                
                %calculate secondary theta
                THETA2L(k) = abs(PHASE2xL(k) - RESULTANTxL(k));
                
                %calculate ratio of hypotenuse to missing side for primary
                RATIO1L(k) = abs(sin(THETA1L(k)));
                
                %calculate ratio of hypotenuse to missing side for
                %secondary
                RATIO2L(k) = abs(sin(THETA2L(k)));
                
                %calculate primary missing side
                DROP1L(k) = RATIO1L(k) * MAG1L(k);
                
                %calculate secondary missing side
                DROP2L(k) = RATIO2L(k) * MAG2L(k);
                               
                %if primary contribution is less than XdB above secondary
                %contribution
                if mag2db(sqrt(MAG1L(k)^2 - DROP1L(k)^2)) - mag2db(sqrt(MAG2L(k)^2 - DROP2L(k)^2)) < atten
                    
                    %calculate the difference in contribution necessary to
                    %achieve XdB gain of primary over secondary
                    
                    if (equal_loudness)
                        
                        DIFFL(k) = sqrt(MAG2L(k)^2 - DROP2L(k)^2) - db2mag(mag2db(sqrt(MAG1L(k)^2 - DROP1L(k)^2)) - (atten - (atten*Vq_half(k))));
                        
                    else
                        
                        DIFFL(k) = sqrt(MAG2L(k)^2 - DROP2L(k)^2) - db2mag(mag2db(sqrt(MAG1L(k)^2 - DROP1L(k)^2)) - atten);
                        
                    end
                    
                    %calculate ratio of drop to difference in
                    %contribution for primary
                    DELTARATIO1L(k) = abs(tan(THETA1L(k)));
                    
                    %calculate ratio of drop to difference in contribution
                    %for secondary
                    DELTARATIO2L(k) = abs(tan(THETA2L(k)));
                    
                    %if primary contribution is less than XdB
                    %if mag2db(sqrt(MAG1(k)^2 - DROP1(k)^2)) < atten
                    
                        %boost primary (hypotenuse)
                    %    MAG1x(k) = MAG1(k) + sqrt((DELTARATIO1(k) * DIFF(k))^2 + DIFF(k)^2);
                        
                        %polar to cartesian
                    %    [re, im] = pol2cart(PHASE1(k), MAG1x(k));
            
                        %rectangular to complex
                    %    WINCHUNK1(k) = re + i*im;
                        
                    %otherwise
                    %else
                        
                        %attenuate secondary
                        MAG2xL(k) = MAG2L(k) - sqrt((DELTARATIO2L(k) * DIFFL(k))^2 + DIFFL(k)^2);
                        
                        %polar to cartesian
                        [re, im] = pol2cart(PHASE2L(k), MAG2xL(k));
            
                        %rectangular to complex
                        WINCHUNK2L(k) = re + 1i*im;
                        
                    %end
                    
                end
                
            end
            
        end
        
        %if primary freq bin dB is greater than threshold
        if mag2db(abs(WINCHUNK1R(k))) > threshold
            
            %get phases
            PHASE1R(k) = angle(WINCHUNK1R(k));
            PHASE2R(k) = angle(WINCHUNK2R(k));
        
            %get magnitudes
            MAG1R(k) = abs(WINCHUNK1R(k));
            MAG2R(k) = abs(WINCHUNK2R(k));
            
            %adjust phase to be from 0 to 2 * pi
            if PHASE1R(k) < 0
                
                PHASE1xR(k) = PHASE1R(k) + (2*pi);
                
            else
                
                PHASE1xR(k) = PHASE1R(k);
                
            end
            
            %adjust phase to be from 0 to 2 * pi
            if PHASE2R(k) < 0
                
                PHASE2xR(k) = PHASE2R(k) + (2*pi);
                
            else
                
                PHASE2xR(k) = PHASE2R(k);
                
            end
            
            %if vectors point in exact opposite directions or same
            %direction
            if abs(PHASE1xR(k) - PHASE2xR(k)) == pi || PHASE1xR(k) == PHASE2xR(k)
                
                %if primary freq bin is less than XdB above secondary
                if mag2db(MAG1R(k)) - mag2db(MAG2R(k)) < atten 
                    
                    %attenuate secondary freq bin such that primary is XdB
                    %above secondary
                    
                    if (equal_loudness)
                        
                        MAG2R(k) = MAG2R(k) - (MAG2R(k) - db2mag(mag2db(MAG1R(k)) - (atten - (atten*Vq_half(k)))));
                        
                    else
                        
                        MAG2R(k) = MAG2R(k) - (MAG2R(k) - db2mag(mag2db(MAG1R(k)) - atten));
                                              
                    end
                    
                    %polar to cartesian
                    [re, im] = pol2cart(PHASE2R(k), MAG2R(k));
            
                    %rectangular to complex
                    WINCHUNK2R(k) = re + 1i*im;
                
                end    
            
            %otherwise
            else
                
                %find resultant magnitude
                MAGNITUDER(k) = sqrt((real(WINCHUNK1R(k))+real(WINCHUNK2R(k)))^2 + (imag(WINCHUNK1R(k))+imag(WINCHUNK2R(k)))^2);
                
                %find resultant angle
                RESULTANTR(k) = atan2((imag(WINCHUNK1R(k))+imag(WINCHUNK2R(k))),(real(WINCHUNK1R(k))+real(WINCHUNK2R(k))));
                
                %adjust the phase of resultant angle to be from 0 to 2 * pi
                if RESULTANTR(k) < 0
                
                    RESULTANTxR(k) = RESULTANTR(k) + (2*pi);
                
                else
                
                    RESULTANTxR(k) = RESULTANTR(k);
                
                end
                
                %calculate primary theta
                THETA1R(k) = abs(PHASE1xR(k) - RESULTANTxR(k));
                
                %calculate secondary theta
                THETA2R(k) = abs(PHASE2xR(k) - RESULTANTxR(k));
                
                %calculate ratio of hypotenuse to missing side for primary
                RATIO1R(k) = abs(sin(THETA1R(k)));
                
                %calculate ratio of hypotenuse to missing side for
                %secondary
                RATIO2R(k) = abs(sin(THETA2R(k)));
                
                %calculate primary missing side
                DROP1R(k) = RATIO1R(k) * MAG1R(k);
                
                %calculate secondary missing side
                DROP2R(k) = RATIO2R(k) * MAG2R(k);
                               
                %if primary contribution is less than XdB above secondary
                %contribution
                if mag2db(sqrt(MAG1R(k)^2 - DROP1R(k)^2)) - mag2db(sqrt(MAG2R(k)^2 - DROP2R(k)^2)) < atten
                    
                    %calculate the difference in contribution necessary to
                    %achieve XdB gain of primary over secondary
                    
                    if (equal_loudness)
                        
                        DIFFR(k) = sqrt(MAG2R(k)^2 - DROP2R(k)^2) - db2mag(mag2db(sqrt(MAG1R(k)^2 - DROP1R(k)^2)) - (atten - (atten*Vq_half(k))));
                        
                    else
                        
                        DIFFR(k) = sqrt(MAG2R(k)^2 - DROP2R(k)^2) - db2mag(mag2db(sqrt(MAG1R(k)^2 - DROP1R(k)^2)) - atten);
                        
                    end
                                          
                    %calculate ratio of drop to difference in
                    %contribution for primary
                    DELTARATIO1R(k) = abs(tan(THETA1R(k)));
                    
                    %calculate ratio of drop to difference in contribution
                    %for secondary
                    DELTARATIO2R(k) = abs(tan(THETA2R(k)));
                   
                    %if primary contribution is less than XdB
                    %if mag2db(sqrt(MAG1(k)^2 - DROP1(k)^2)) < atten
                    
                        %boost primary (hypotenuse)
                    %    MAG1x(k) = MAG1(k) + sqrt((DELTARATIO1(k) * DIFF(k))^2 + DIFF(k)^2);
                        
                        %polar to cartesian
                    %    [re, im] = pol2cart(PHASE1(k), MAG1x(k));
            
                        %rectangular to complex
                    %    WINCHUNK1(k) = re + i*im;
                        
                    %otherwise
                    %else
                        
                        %attenuate secondary
                        MAG2xR(k) = MAG2R(k) - sqrt((DELTARATIO2R(k) * DIFFR(k))^2 + DIFFR(k)^2);
                       
                        %polar to cartesian
                        [re, im] = pol2cart(PHASE2R(k), MAG2xR(k));
            
                        %rectangular to complex
                        WINCHUNK2R(k) = re + 1i*im;
                        
                    %end
                    
                end
                
            end
            
        end
          
    end
    
%take inverse fft
winchunk1L = real(ifft(WINCHUNK1L));
winchunk1R = real(ifft(WINCHUNK1R));
winchunk2L = real(ifft(WINCHUNK2L));
winchunk2R = real(ifft(WINCHUNK2R));
    
%apply window
winchunk1L = sinwin .* winchunk1L;
winchunk1R = sinwin .* winchunk1R;
winchunk2L = sinwin .* winchunk2L;
winchunk2R = sinwin .* winchunk2R;
    
%overlap and add
output1(startingSamp1:startingSamp - 1, 1) = output1(startingSamp1:startingSamp - 1, 1) + winchunk1L;
output1(startingSamp1:startingSamp - 1, 2) = output1(startingSamp1:startingSamp - 1, 2) + winchunk1R;
output2(startingSamp1:startingSamp - 1, 1) = output2(startingSamp1:startingSamp - 1, 1) + winchunk2L;
output2(startingSamp1:startingSamp - 1, 2) = output2(startingSamp1:startingSamp - 1, 2) + winchunk2R;

%trim outputs
output1=output1(1:end-pad, :);
output2=output2(1:end-pad, :);

%add the two signals in the time domain
sum = output1 + output2;

%output the result
%plot(output1);
%plot(output2);
audiowrite('output3.wav',output2,fs);
audiowrite('sum2new3.wav',sum,fs);

end