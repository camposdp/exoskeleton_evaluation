%% Processar os dados de CONTROLE (sem exoesqueleto)
clear; clc;

% Caminho base onde estão as pastas
base_path = "C:\_PESQUISA\JP\Coleta";

% Listar os 3 dias de controle
controles = ["Coleta-01-sem-exo","Coleta-02-sem-exo","Coleta-03-sem-exo"];

% Inicializar variáveis
control_dyn_ombro = [];
control_dyn_biceps = [];
control_p_dyn_ombro = [];
control_p_dyn_biceps = [];
control_stat_ombro = [];
control_stat_biceps = [];

% Parâmetros de processamento
fs = 1000;
low = 20;
high = 350;
[b,a] = butter(2, [low, high] / (fs/2), 'bandpass');
notch_freq = [59, 61];
[d,c] = butter(2, notch_freq / (fs/2), 'stop');

for i = 1:numel(controles)
    % Pasta do controle atual
    pasta = fullfile(base_path, controles(i));
    
    % Carregar MCV
    ombro_data  = readmatrix(fullfile(pasta, 'mcv-ombro.txt'));
    biceps_data = readmatrix(fullfile(pasta, 'mcv-biceps.txt'));
    ombro = filtfilt(b, a, filtfilt(d, c, ombro_data(:,6)));
    biceps = filtfilt(b, a, filtfilt(d, c, biceps_data(:,6)));
    MCV_ombro = max(suavizarEMG(ombro, fs));
    MCV_biceps = max(suavizarEMG(biceps, fs));
    
    % Carregar Dinâmico
    dinamico_data = readmatrix(fullfile(pasta, 'dinamica.txt'));
    sinal_biceps = filtfilt(b, a, filtfilt(d, c, dinamico_data(:,6)));
    sinal_ombro  = filtfilt(b, a, filtfilt(d, c, dinamico_data(:,7)));
    
    % Normalizar sinais (dividindo pelo MCV)
    norm_ombro = suavizarEMG(sinal_ombro, fs) / MCV_ombro;
    norm_biceps = suavizarEMG(sinal_biceps, fs) / MCV_biceps;
    
    % Calcular RMS por pulso
    [~, rms_ombro] = extrairMCVoRMS(sinal_ombro, fs, MCV_ombro, true);
    [~, rms_biceps] = extrairMCVoRMS(sinal_biceps, fs, MCV_biceps, true);

    % Calcular Percentis (5,25,50,75,95) dos sinais normalizados
    p_ombro = prctile(norm_ombro, [5 25 50 75 95]);
    p_biceps = prctile(norm_biceps, [5 25 50 75 95]);
    
    % Carregar Estático
    estatico_data = readmatrix(fullfile(pasta, 'isometrica.txt'));
    estatico_biceps = filtfilt(b, a, filtfilt(d, c, estatico_data(:,6)));
    estatico_ombro  = filtfilt(b, a, filtfilt(d, c, estatico_data(:,7)));

    [rms_estatico_ombro, ~] = extrairMCVoRMS(estatico_ombro, fs, MCV_ombro, false);
    [rms_estatico_biceps, ~] = extrairMCVoRMS(estatico_biceps, fs, MCV_biceps, false);


    % Normalizar estaticos
    %norm_estatico_ombro  = suavizarEMG(estatico_ombro, fs) / MCV_ombro;
    %norm_estatico_biceps = suavizarEMG(estatico_biceps, fs) / MCV_biceps;

    % Cálculo de RMS estático (em toda a duração do sinal)
    %rms_estatico_ombro = rms(norm_estatico_ombro);
    %rms_estatico_biceps = rms(norm_estatico_biceps);
    
    % Concatenar (adicionando aos vetores finais)
    control_dyn_ombro = [control_dyn_ombro, rms_ombro(1:10)];
    control_dyn_biceps = [control_dyn_biceps, rms_biceps(1:10)];
    
    control_p_dyn_ombro = [control_p_dyn_ombro; p_ombro];
    control_p_dyn_biceps = [control_p_dyn_biceps; p_biceps];

    control_stat_ombro = [control_stat_ombro, rms_estatico_ombro];
    control_stat_biceps = [control_stat_biceps, rms_estatico_biceps];
end

% Cálculo das médias dos percentis
mean_control_p_dyn_ombro = mean(control_p_dyn_ombro, 1);
mean_control_p_dyn_biceps = mean(control_p_dyn_biceps, 1);

% Salvar
save(fullfile(base_path, 'controle.mat'), ...
    'control_dyn_ombro', 'control_dyn_biceps', ...
    'control_p_dyn_ombro', 'control_p_dyn_biceps', ...
    'control_stat_ombro', 'control_stat_biceps', ...
    'mean_control_p_dyn_ombro', 'mean_control_p_dyn_biceps');

fprintf('Controle processado e salvo.\n');



%% --- Função auxiliar para suavizar
function sinal_suavizado = suavizarEMG(sinal_bruto, fs)
    sinal_retificado = abs(sinal_bruto);
    [b, a] = butter(4, 6 / (fs/2), 'low');
    sinal_suavizado = filtfilt(b, a, sinal_retificado);
end

%% --- Função auxiliar para extrair RMS de pulsos
function [valor_rms, rms_dinamico] = extrairMCVoRMS(sinal, fs, MCV, dinamico)
    if nargin < 4
        dinamico = false;
    end

    % Normalize imediatamente
    sinal_normalizado = suavizarEMG(sinal, fs) / MCV;
    
    baseline = sinal_normalizado(1:500);
    u = mean(baseline);
    s = std(baseline);
    th = u + 5*s;

    vetor_binario_bruto = aplicarLimiarDuplo(sinal_normalizado, fs, th, th, 500);
    pulso_dilatado = filtroDilatacao(vetor_binario_bruto, 300);
    vetor_binario = pulso_dilatado(:)';
    
    if dinamico
        try
            pulsos = identificarPulsos(vetor_binario);
            rms_dinamico = zeros(1, size(pulsos,1));
            for i = 1:size(pulsos,1)
                rms_dinamico(i) = rms(sinal_normalizado(pulsos(i,1):pulsos(i,2)));
            end
            valor_rms = NaN;
        catch
            warning('Erro na detecção dos pulsos. Retornando vetor NaN.');
            rms_dinamico = NaN(1,10);
            valor_rms = NaN;
        end
    else
        onset = vetor_binario;
        valor_rms = rms(sinal_normalizado(onset>0));
        rms_dinamico = [];
    end
end


%% --- Função auxiliar para aplicar limiar duplo
function vetor_binario_duplo = aplicarLimiarDuplo(sinal_suavizado, fs, threshold_low, threshold_high, minDurationMs)
    N = length(sinal_suavizado);
    vetor_binario_duplo = zeros(1, N);
    minDurationSamples = round(minDurationMs * fs / 1000);
    
    i = 1;
    while i <= N
        if sinal_suavizado(i) > threshold_high
            startIdx = i;
            j = i;
            while j <= N && (sinal_suavizado(j) > threshold_low)
                j = j + 1;
            end
            endIdx = j - 1;
            if (endIdx - startIdx + 1) >= minDurationSamples
                vetor_binario_duplo(startIdx:endIdx) = 1;
            end
            i = endIdx + 1;
        else
            i = i + 1;
        end
    end
end

%% --- Função auxiliar para dilatação morfológica
function pulso_dilatado = filtroDilatacao(pulso_binario, windowSize)
    pulso_binario = pulso_binario(:);
    elementoEstruturante = ones(windowSize, 1);
    sinal_dilatado = conv(double(pulso_binario), elementoEstruturante, 'same');
    pulso_dilatado = sinal_dilatado > 0;
end

%% --- Função para identificar pulsos
function pulsos = identificarPulsos(sinal_binario)
    if ~all(ismember(sinal_binario, [0 1]))
        error('O sinal de entrada deve ser binário.');
    end
    diferenca = diff([0, sinal_binario, 0]);
    inicio = find(diferenca == 1);
    fim = find(diferenca == -1) - 1;
    pulsos = [inicio(:), fim(:)];
end
