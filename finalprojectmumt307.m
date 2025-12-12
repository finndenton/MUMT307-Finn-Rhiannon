% 
%
fs = 44100;
nBits = 16;
nChannels = 1;
recObj = audiorecorder(fs, nBits, nChannels);
choice = input("Would you like to sing a note? Y/N: ", "s");
if choice == "Y"
    disp('Sing a note!')
    recordblocking(recObj, 5);
    disp('End of recording.');
    fprintf("\n")
else
    disp("See you next time!")
    return;
end

y = getaudiodata(recObj);

% yin estimator for pitch tracking
% yin_estimator.m by Orchisama Das

function [time, f0] = yin_estimator(x, fs, varargin)

%function that implements YIN algorithm for
%fundamental pitch tracking
%x - input audio signal
%fs - sampling rate
%time,f0 - time vector and associated fundamental frequencies estimated

%window size -  we assume the minimum f0 to be 1/0.025 = 40Hz
win = round(0.025*fs);
N = length(x);
nframes = ceil(N/win);
%zero pad signal to have enough frames
x = [x, zeros(1,win*nframes - N)];
x_frame = zeros(nframes, win);
start = 1;
%break into windows
for i = 1:nframes
    x_frame(i,:) = x(start:start + win - 1);
    start = start + win;
end

%step 1 - calculate difference function 
d = zeros(nframes,win);
x_temp = [x_frame, zeros(nframes,win)];
for tau = 0:win-1
    for j = 1:win  
         d(:,tau+1) = d(:,tau+1) + (x_temp(:,j) - x_temp(:,j+tau)).^2;         
    end
end


%step 2 - cumulative mean normalised difference function
d_norm = zeros(nframes,win);
d_norm(:,1) = 1;

for i = 1:nframes
    for tau = 1:win-1
        d_norm(i,tau+1) = d(i,tau+1)/((1/tau) * sum(d(i,1:tau+1)));
    end
end

%step 3 - absolute thresholding
lag = zeros(1,nframes);
th = 0.1;
for i = 1:nframes
    l = find(d_norm(i,:) < th,1);
    if(isempty(l) == 1)
        [v,l] = min(d_norm(i,:));
    end
    lag(i) = l;
    
end

%step 4 - parabolic interpolation
period = zeros(1,nframes);
time = zeros(nframes,win);
f0 = zeros(nframes,win);
start = 1;

for i = 1:nframes
    if(lag(i) > 1 && lag(i) < win)
        alpha = d_norm(i,lag(i)-1);
        beta = d_norm(i,lag(i));
        gamma = d_norm(i,lag(i)+1);
        peak = 0.5*(alpha - gamma)/(alpha - 2*beta + gamma);
        %ordinate needs to be calculated from d and not d_norm - see paper
        %ordinate = d(i,lag(i)) - 0.25*(d(i,lag(i)-1) - d(i,lag(i)+1))*peak;
    else
        peak = 0;
    end
    %1 needs to be subtracted from 1 due to matlab's indexing nature
    period(i) = (lag(i)-1) + peak;
    f0(i,:) = fs/period(i)*ones(1,win);
    time(i,:) = ((i-1)*win:i*win-1)/fs;
    
end

%for silent frames estimated frequency should be 0Hz
if ~isempty(varargin)
    [f0] = silent_frame_classification(x_frame, f0);
end

f0 = reshape(f0',1,nframes*win);
time = reshape(time',1,nframes*win);

end


% calling the yin_estimator function
[time, f0] = yin_estimator(y', fs);

% number of samples
%fprintf('Samples: %d\n', numel(f0));

% avoid discontinuities in sound recording
samplesToSkip = round(0.2 * fs);
f0(1:samplesToSkip) = NaN;sm

% reject impossible pitch values
f0(f0 < 40 | f0 > 2000) = NaN;

f0_clean = f0;
for k = 2:length(f0_clean)
    if f0_clean(k) <= 0 || isnan(f0_clean(k))
        f0_clean(k) = f0_clean(k-1);   % hold previous pitch
    end
end

f0_smooth = movmedian(f0_clean, 50);

validF0 = f0_smooth(~isnan(f0_smooth));
f0_note = mean(validF0);          % ex. ~260 Hz

fprintf('Estimated note frequency in Hz: %.2f Hz\n', f0_note);

midi = round(69 + 12*log2(f0_note/440));
root = 440 * 2^((midi - 69)/12);

fprintf("Your note: %s\n", freqToNote(root));
fprintf("\n")


function nameOfNote = freqToNote(freq)
    midi = round(12*log2(freq/440) + 69);

    noteNames = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"];
    
    idx = mod(midi, 12) + 1;

    nameOfNote = noteNames(idx);
end


% major chord
major3 = root*5/4;
perfect5 = root*3/2;

majorChord = [root, major3, perfect5];

noteNames = strings(1,3);
for i = 1:3
    noteNames(i) = freqToNote(majorChord(i));
end

fprintf("Major triad with your note as root: %s-%s-%s\n", noteNames(1), noteNames(2), noteNames(3));
disp("The intervals are: Major 3rd, Perfect 5th");
fprintf("\n")

% minor chord
minor3 = root*6/5;

minorChord = [root, minor3, perfect5];

noteNames = strings(1,3);
for i = 1:3
    noteNames(i) = freqToNote(minorChord(i));
end

fprintf("Minor triad with your note as root: %s-%s-%s\n", noteNames(1), noteNames(2), noteNames(3));
disp("The intervals are: Minor 3rd, Perfect 5th");
fprintf("\n")


% augmented
augment5 = root*11/7;

augmentedChord = [root, major3, augment5];

noteNames = strings(1,3);
for i = 1:3
    noteNames(i) = freqToNote(augmentedChord(i));
end

fprintf("Augmented triad with your note as root: %s-%s-%s\n", noteNames(1), noteNames(2), noteNames(3));
disp("The intervals are: Major 3rd, Augmented 5th");
fprintf("\n")

% dimished
dim5 = root*1024/729;

diminishedChord = [root, minor3, dim5];

noteNames = strings(1,3);
for i = 1:3
    noteNames(i) = freqToNote(diminishedChord(i));
end

fprintf("Diminshed triad with your note as root: %s-%s-%s\n", noteNames(1), noteNames(2), noteNames(3));
disp("The intervals are: Minor 3rd, Diminished 5th");
fprintf("\n")

% time vector length of recording
N = length(y);
t = (0:N-1)/fs;

% Generate a pure tone at that pitch
y_tone = 0.8 * sin(2*pi*f0_note * t);

% Pure tone samples
% fprintf('Pure tone samples: %d\n', numel(y_tone));

sound(y_tone, fs);

combined = y_tone;   % Start with the original tone

while true

    combined = y_tone;
    
    disp("0 - Play Sound / Exit");
    disp("1 - Harmonize");
    disp("2 - Sound Effect");
    sound_choice = input("Choose a method from above (0-2): ");

    if sound_choice == 1
        while true
            disp("0 - Play Sound and Exit")
            disp("1 - Octave");
            disp("2 - Major 2nd");
            disp("3 - Minor 2nd");
            disp("4 - Major 3rd");
            disp("5 - Minor 3rd");
            disp("6 - Perfect 4th");
            disp("7 - Perfect 5th");
            disp("8 - Major 6th");
            disp("9 - Minor 6th");
            disp("10 - Major 7th");
            disp("11 - Minor 7th");
            disp("12 - Augmented 5th");
            disp("13 - Diminished 5th");

            harmonyStep = input("Choose a harmony number (1-13) or Exit(0): ");
        
            if harmonyStep == 0
                sound(combined, fs);
                disp("Playing sound...");
                break;
            end
        
            % Generate the harmonic frequency
            if harmonyStep == 1
                octave = 2*f0_note;
                y_harmony = 0.8 * sin(2*pi*octave * t);
            elseif harmonyStep == 2
                maj2 = 9/8*f0_note;
                y_harmony = 0.8 * sin(2*pi*maj2 * t);
            elseif harmonyStep == 3
                min2 = 16/15*f0_note;
                y_harmony = 0.8 * sin(2*pi*min2 * t);
            elseif harmonyStep == 4
                third = 5/4*f0_note;
                y_harmony = 0.8 * sin(2*pi*third * t);
            elseif harmonyStep == 5
                minorthird = 6/5*f0_note;
                y_harmony = 0.8 * sin(2*pi*minorthird * t);
            elseif harmonyStep == 6
                pfourth = 4/3*f0_note;
                y_harmony = 0.8 * sin(2*pi*pfourth * t);
            elseif harmonyStep == 7
                pfifth = 3/2*f0_note;
                y_harmony = 0.8 * sin(2*pi*pfifth * t);
            elseif harmonyStep == 8
                majsix = 5/3*f0_note;
                y_harmony = 0.8 * sin(2*pi*majsix * t);
            elseif harmonyStep == 9
                minsix = 8/5*f0_note;
                y_harmony = 0.8 * sin(2*pi*minsix * t);
            elseif harmonyStep == 10
                maj7 = 15/8*f0_note;
                y_harmony = 0.8 * sin(2*pi*maj7 * t);
            elseif harmonyStep == 11
                min7 = 7/4*f0_note;
                y_harmony = 0.8 * sin(2*pi*min7 * t);
            elseif harmonyStep == 12
                aug = f0_note * 2^(1/12);
                y_harmony = 0.8 * sin(2*pi*aug * t);
            elseif harmonyStep == 13
                dim = f0_note / 2^(1/12);
                y_harmony = 0.8 * sin(2*pi*dim * t);
            else
                disp("Please write the input correctly");
                continue;
            end
        
            % Add to layered signal
           combined = (combined + y_harmony)/2;
        
        
            disp("Harmony added!");
    
        end

    elseif sound_choice == 2
        
        while true

            %sound_effect = y_tone;
            sound_effect = y;

            disp("0 - Exit");
            disp("1 - Phaser");
            disp("2 - Ring Modulator");
            disp("3 - Distortion");
            disp("4 - Chorus");
            disp("5 - Vibrato");
            effect = input("Choose an effect(1-5) or Exit(0): ");

            if effect == 0
                    break;
            end

            if effect == 1
                lfo = 0.2; % sweep
                depth = 0.7;
                base_fc = 800;
                Q = 4; % "resonance" Q
               
                %N = length(sound_effect);
                x = y(:).'; 
                N = length(x);
                t = (0:N-1)/fs;
                
                Ts = 1/fs;
                r = exp(-pi*base_fc*Ts/Q); % pole radius
        
                modfc = base_fc.*(1+depth*sin(2*pi*lfo*t));
                
                %inputSound = y;
                sound_effect = zeros(1,N); % output sound vector
                zi = zeros(1,2); % filter function state vector
                a2 = r^2;
            
                for n = 1:N
                    a1 = -2*r*cos(2*pi*modfc(n)*Ts);
        
                    [tmp, zi] = filter([a2 a1 1],[1, a1, a2], x(n), zi);
        
                    sound_effect(n) = 0.7*x(n) + 0.7*tmp;
                end

                sound_effect = sound_effect / max(abs(sound_effect)+1e-12);
            
                disp("Playing original sound...")
                sound(y, fs);
                pause(6);
                disp("Playing tone with phaser...")
                sound(sound_effect, fs);

            elseif effect == 2

                mod_freq = 40;
                wet = 0.7;
            
                N = length(sound_effect);
                t = (0:N-1)/fs;
            
                % multiply by sin
                carrier = sin(2*pi*mod_freq * t);
            
                sound_effect = sound_effect(:).';  
                carrier = carrier(:).';             

                % apply modulation with carrier freq and original sound
                modulated = sound_effect .* carrier;
           
                sound_effect = (1-wet)*sound_effect + wet*modulated;
            
                %normalize
                sound_effect = sound_effect / max(abs(sound_effect));
            
                disp("Playing original sound...");
                sound(y, fs);
                pause(6);
                disp("Playing tone with ring modulation");
                sound(sound_effect, fs);

            elseif effect == 3

                x = y(:).';      
                gain   = 4;   
                clip = 0.6;

              
                d = gain * x;
                d = max(-clip, min(clip, d));

                d = d / clip;

                d = d / max(abs(d) + 1e-12);

                sound_effect = d;
        
                disp("Playing original sound...")
                sound(y, fs);
                pause(6);
                disp("Playing tone with distortion...")
                sound(sound_effect, fs);

            elseif effect == 4

                wet = 0.5;
                delay_ms = 20;
                depth_ms = 10;
                rate_hz = 1.5;
            
                delay_samples = round(delay_ms * fs / 1000);
                depth_samples = round(depth_ms * fs / 1000);
            
                lfo = depth_samples * sin(2*pi*rate_hz * (0:length(sound_effect)-1) / fs);
            
                chorus = zeros(size(sound_effect));
            
                for n = 1:length(sound_effect)
                    d = delay_samples + round(lfo(n));
                    if n - d > 0
                        chorus(n) = sound_effect(n - d);
                    end
                end
            
                sound_effect = sound_effect + wet * chorus;
                sound_effect = sound_effect / max(abs(sound_effect));

                disp("Playing original sound...")
                sound(y, fs);
                pause(6);
                disp("Playing tone with chorus...")
                sound(sound_effect, fs);

            elseif effect == 5
                vib_rate = 5;
                vib_depth = 3;
                wet = 0.8;

                x = y(:).';
                N = length(x);
                t = (0:N-1) / fs;

                depth_samples = vib_depth*fs*1e-3; % miliseconds to seconds

                Delay = depth_samples * (sin(2*pi*vib_rate * t));

                vibrato = zeros(1, N);

                for n = 1:N
                    idx = n-Delay(n);

                    if idx < 1
                        vibrato(n) = 0;
                    else
                        index = floor(idx);
                        frac = idx - index;

                        if index+1<=N
                            vibrato(n) = (1-frac)*x(index)+ frac*x(index+1);

                        else
                            vibrato(n) = x(index);
                        end

                    end
                end

                sound_effect = (1-wet)*x + wet*vibrato;
                sound_effect = sound_effect / max(abs(sound_effect)+1e-12);
                
                disp("Playing original sound...")
                sound(y, fs);
                pause(6);
                disp("Playing tone with vibrato...")
                sound(sound_effect, fs);

            end

        end

    elseif sound_choice == 0
        sound(y, fs);
        disp("Playing original sound...");
        break;
    else
        disp("Please write the input correctly: either 0, 1, or 2");
        continue;
    end

end

