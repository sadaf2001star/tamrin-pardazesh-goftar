% پارامترهای ضبط
fs = 44100; % نرخ نمونه‌برداری (هرتز)
nBits = 16; % تعداد بیت‌ها در هر نمونه
nChannels = 1; % تعداد کانال‌ها (1 برای مونو)
words = {'balee', 'nah', 'salam', 'khodahafez', 'lotfan', 'tashakor', 'bebakhshid', 'komak', 'tavaghof', 'boro', 'chap', 'rast', 'bala', 'paeen', 'shoroe', 'payan', 'baz', 'baste', 'roshan', 'khamosh'};
numRepetitions = 5;
numAugmentations = 10;
outputDirectory = 'recorded_audio';

% ایجاد دایرکتوری خروجی در صورت عدم وجود
if ~exist(outputDirectory, 'dir')
 mkdir(outputDirectory);
end

% ضبط ویس‌ها
for i = 1:length(words)
 for j = 1:numRepetitions
 recObj = audiorecorder(fs, nBits, nChannels);
 disp(['Please say ', words{i}, ' five times.']);
 recordblocking(recObj, 2);
 disp('End of Recording.');
 audioData = getaudiodata(recObj);
 filename = fullfile(outputDirectory, [words{i}, '_', sprintf('%02d', j), '.wav']);
 audiowrite(filename, audioData, fs);
 disp(['Audio recorded and saved to ', filename]);
 pause(2); % فاصله 2 ثانیه‌ای بین ضبط‌ها
 end
end
numAugmentations = 10;
outputDirectory = 'recorded_audio';

% افزایش داده‌ها
for i = 1:length(words)
 for j = 1:numRepetitions
 inputFilePath = fullfile(outputDirectory, [words{i}, '_', sprintf('%02d', j), '.wav']);
 augmentedFiles = augmentAudio(inputFilePath, numAugmentations, outputDirectory);
 disp(['Augmented files created for ', inputFilePath]);
 end
end

function augmentedAudioFiles = augmentAudio(filePath, numAugmentations, outputDir)
 [audioData, fs] = audioread(filePath);
 augmentationFunctions = {@changeSpeed, @addNoise, @shiftTime};
 augmentedAudioFiles = cell(numAugmentations, 1);
 for i = 1:numAugmentations
 funcIndex = mod(i,3)+1;
 augmentedAudio = augmentationFunctions{funcIndex}(audioData, fs);
 [~, name, ~] = fileparts(filePath);
 newFileName = fullfile(outputDir, [name, '_aug_', sprintf('%02d', i), '.wav']);
 audiowrite(newFileName, augmentedAudio, fs);
 augmentedAudioFiles{i} = newFileName;
 end
end
% توابع افزایش داده‌ها
function augmentedAudio = changeSpeed(audioData, fs)
    % تغییر سرعت پخش (با ضریب بین 0.8 تا 1.2)
    speedFactor = 0.8 + (1.2-0.8).*rand(1,1);
    augmentedAudio = resample(audioData, round(fs*speedFactor), fs);
end
function augmentedAudio = addNoise(audioData, fs)
    % افزودن نویز سفید با سطح بین 0.001 تا 0.005
    noiseLevel = 0.001 + (0.005-0.001).*rand(1,1);
    noise = noiseLevel * randn(size(audioData));
    augmentedAudio = audioData + noise;
end
function augmentedAudio = shiftTime(audioData, fs)
    % جابجایی زمانی (با حداکثر جابجایی 0.1 ثانیه)
    shiftTime = -0.1 + (0.1-(-0.1)).*rand(1,1);
    shiftSamples = round(shiftTime * fs);
    augmentedAudio = circshift(audioData, shiftSamples);
end

function mfccs = computeMFCC(signal, fs, numCoeffs)
 % تنظیم پارامترهای اولیه
 frameLength = round(0.025 * fs); % طول هر فریم 25 میلی‌ثانیه
 frameStep = round(0.01 * fs); % قدم هر فریم 10 میلی‌ثانیه
 signalLength = length(signal);
 % محاسبه تعداد فریم‌ها
 numFrames = 1 + floor((signalLength - frameLength) / frameStep);
 % ایجاد پنجره هموارکننده (Hamming)
 hammingWindow = hamming(frameLength);
 % آماده‌سازی ماتریس MFCC
 mfccs = zeros(numFrames, numCoeffs);
 % حلقه برای محاسبه MFCC هر فریم
 for i = 1:numFrames
 startIdx = (i - 1) * frameStep + 1;
 endIdx = startIdx + frameLength - 1;
 % اگر اندیس‌ها خارج از محدوده سیگنال باشند، حلقه را بشکنید
 if endIdx > signalLength
 break;
 end
 frame = signal(startIdx:endIdx) .* hammingWindow;
 % تبدیل فوریه سریع (FFT)
 frameFFT = fft(frame);
 magFrame = abs(frameFFT(1:frameLength/2+1));
 melSpectrum = melFilter(magFrame, fs);
 logMelSpectrum = log(melSpectrum + 1e-10);
 % تبدیل کسینوسی معکوس (DCT) برای به دست آوردن MFCC
 mfcc = dct(logMelSpectrum);
 % استخراج اولین numCoeffs ضرایب MFCC
 mfccs(i, :) = mfcc(1:numCoeffs);
 end
end

function melSpectrum = melFilter(spectrum, fs)
 numFilters = 26;
 NFFT = length(spectrum) * 2 - 2;
 % ایجاد فیلترهای Mel
 melPoints = linspace(0, hz2mel(fs / 2), numFilters + 2);
 hzPoints = mel2hz(melPoints);
 bin = floor((NFFT + 1) * hzPoints / fs);
 bin(bin < 1) = 1;
 bin(bin > NFFT / 2 + 1) = NFFT / 2 + 1;
 melFilterBank = zeros(numFilters, floor(NFFT / 2 + 1));
 for i = 1:numFilters
 melFilterBank(i, bin(i):bin(i+1)) = linspace(0, 1, bin(i+1) - bin(i) + 1);
 melFilterBank(i, bin(i+1):bin(i+2)) = linspace(1, 0, bin(i+2) - bin(i+1) + 1);
 end
 melSpectrum = melFilterBank * spectrum;
end

function mel = hz2mel(hz)
 mel = 2595 * log10(1 + hz / 700);
end

function hz = mel2hz(mel)
 hz = 700 * (10.^(mel / 2595) - 1);
end

function label = extractLabel(fileName)
 % این تابع برچسب را از نام فایل استخراج می‌کند
 % فرض بر این است که نام فایل‌ها به صورت word_number.wav هستند
 tokens = split(fileName, '_');
 word = tokens{1};
 % اینجا باید یک مپ از کلمات به برچسب‌های عددی داشته باشید
 % مثال:
 labelsMap = containers.Map({'balee', 'nah', 'salam', 'khodahafez', 'lotfan', 'tashakor', 'bebakhshid', 'komak', 'tavaghof', 'boro', 'chap', 'rast', 'bala', 'paeen', 'shoroe', 'payan', 'baz', 'baste', 'roshan', 'khamosh'}, 1:20);
 label = labelsMap(word);
end
% تعریف مسیر و خواندن فایل‌های صوتی
filePattern = fullfile('augmented_audio', '*.wav');
audioFiles = dir(filePattern);
numFiles = length(audioFiles);

% ایجاد دیتاست و شافل کردن داده‌ها
fileList = cell(numFiles, 1);
for k = 1:numFiles
 fileList{k} = fullfile('augmented_audio', audioFiles(k).name);
end
fileList = fileList(randperm(numFiles)); % شافل کردن فایل‌ها

% تقسیم دیتاست به داده‌های آموزشی و تست
numTrainFiles = 950;
trainFiles = fileList(1:numTrainFiles);
testFiles = fileList(numTrainFiles+1:end);

% استخراج ویژگی‌ها و برچسب‌ها
features = [];
labels = [];
for i = 1:length(trainFiles)
 [audioData, fs] = audioread(trainFiles{i});
 % پیش‌پردازش و استخراج ویژگی‌های MFCC
 mfccFeature = computeMFCC(audioData, fs, 13);
 features = [features; mfccFeature];
 % استخراج برچسب از نام فایل
 [~, name, ~] = fileparts(trainFiles{i});
 label = extractLabel(name); % تابعی که برچسب را از نام فایل استخراج می‌کند
 labels = [labels; label];
end

% تعریف معماری شبکه عصبی
inputLayerSize = size(features, 2);
numClasses = 20; % تعداد کلاس‌های مختلف
layers = [
 sequenceInputLayer(inputLayerSize)
 lstmLayer(100, 'OutputMode', 'sequence')
 fullyConnectedLayer(50)
 dropoutLayer(0.2)
 fullyConnectedLayer(numClasses)
 softmaxLayer
 classificationLayer];

% تنظیمات آموزش
options = trainingOptions('adam', ...
 'MaxEpochs', 30, ...
 'MiniBatchSize', 16, ...
 'InitialLearnRate', 1e-3, ...
 'Plots', 'training-progress');

% آموزش شبکه
net = trainNetwork(features, labels, layers, options);

% تست شبکه
testFeatures = [];
testLabels = [];
for i = 1:length(testFiles)
 [audioData, fs] = audioread(testFiles{i});
 mfccFeature = computeMFCC(audioData, fs, 13);
 testFeatures = [testFeatures; mfccFeature];
 [~, name, ~] = fileparts(testFiles{i});
 label = extractLabel(name);
 testLabels = [testLabels; label];
end

% ارزیابی شبکه
predictedLabels = classify(net, testFeatures);
accuracy = sum(predictedLabels == testLabels) / numel(testLabels);
disp(['Test Accuracy: ', num2str(accuracy)]);

