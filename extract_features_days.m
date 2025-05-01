%% Script para processar todas as coletas (1 a 30) com extração de percentis
clear
clc
close all

base_path = "C:\_PESQUISA\JP\Coleta";
fs = 1000;
low = 20;
high = 350;
[b,a] = butter(2, [low, high] / (fs/2), 'bandpass');
notch_freq = [59, 61];
[d,c] = butter(2, notch_freq / (fs/2), 'stop');

for i = 1:30
    coleta_id = sprintf('%02d', i);
    fprintf('Processando Coleta %s...\n', coleta_id);

    %% MCV Ombro
    path_mcv_ombro = fullfile(base_path, ['Coleta-' coleta_id], ['mcv-ombro-' coleta_id '.txt']);
    ombro_data = readmatrix(path_mcv_ombro);
    ombro = filtfilt(b, a, filtfilt(d, c, ombro_data(:,6)));
    MCV_ombro = max(suavizarEMG(ombro, fs));

    %% MCV Biceps
    path_mcv_biceps = fullfile(base_path, ['Coleta-' coleta_id], ['mcv-biceps-' coleta_id '.txt']);
    biceps_data = readmatrix(path_mcv_biceps);
    biceps = filtfilt(b, a, filtfilt(d, c, biceps_data(:,6)));
    MCV_biceps = max(suavizarEMG(biceps, fs));

    %% RMS Estático
    path_estatico = fullfile(base_path, ['Coleta-' coleta_id], "isometrica-" + coleta_id + ".txt");
    dados_estatico = readmatrix(path_estatico);
    ombro_estatico = filtfilt(b, a, filtfilt(d, c, dados_estatico(:,7)));
    biceps_estatico = filtfilt(b, a, filtfilt(d, c, dados_estatico(:,6)));

    [RMS_estatico_ombro, ~] = extrairMCVoRMS(ombro_estatico,fs,MCV_ombro,false);
    [RMS_estatico_biceps, ~] = extrairMCVoRMS(biceps_estatico,fs,MCV_biceps,false);

    %% RMS Dinâmico e normalização
    path_dinamico = fullfile(base_path, ['Coleta-' coleta_id], "dinamica-" + coleta_id + ".txt");
    dados_dinamico = readmatrix(path_dinamico);
    sinal1 = dados_dinamico(:,6);
    sinal2 = dados_dinamico(:,7);
    % Substituir inválidos
    sinal1(~isfinite(sinal1)) = 0;
    sinal2(~isfinite(sinal2)) = 0;

    sinal_notch1 = filtfilt(d, c, sinal1);
    biceps_dinamico = filtfilt(b, a, sinal_notch1);
    sinal_notch2 = filtfilt(d, c, sinal2);
    ombro_dinamico = filtfilt(b, a, sinal_notch2);


    %% Extração de Percentis (5,25,50,75,95)
    p_dyn_ombro   = prctile(suavizarEMG(ombro_dinamico, fs)/MCV_ombro,  [5 25 50 75 95]);
    p_dyn_biceps  = prctile(suavizarEMG(biceps_dinamico, fs)/MCV_biceps, [5 25 50 75 95]);

    %% Cálculo de RMS Dinâmico por pulso
    [~, rms_dinamico_ombro]  = extrairMCVoRMS(ombro_dinamico,  fs, MCV_ombro,true);
    [~, rms_dinamico_biceps] = extrairMCVoRMS(biceps_dinamico, fs, MCV_biceps,true);
    rms_dinamico_ombro  = rms_dinamico_ombro(1:10);
    rms_dinamico_biceps = rms_dinamico_biceps(1:10);



    %% Salvar resultados com percentis
    nome_saida = fullfile(base_path, ['Coleta-' coleta_id], ['resultados_emg_' coleta_id '.mat']);
    save(nome_saida, ...
        'MCV_ombro', 'MCV_biceps', ...
        'RMS_estatico_ombro', 'RMS_estatico_biceps', ...
        'rms_dinamico_ombro', 'rms_dinamico_biceps', ...
        'p_dyn_ombro', 'p_dyn_biceps');
end

fprintf('Processamento finalizado para todas as coletas\n');


%% --- Funções Auxiliares ---
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





function sinal_suavizado = suavizarEMG(sinal_bruto, fs)
    sinal_retificado = abs(sinal_bruto);
    [b, a] = butter(4, 6 / (fs / 2), 'low');
    sinal_suavizado = filtfilt(b, a, sinal_retificado);
end

function pulsos = identificarPulsos(sinal_binario)
    if ~all(ismember(sinal_binario, [0 1]))
        error('O sinal de entrada deve ser binário.');
    end
    diferenca = diff([0, sinal_binario, 0]);
    inicio = find(diferenca == 1);
    fim = find(diferenca == -1) - 1;
    pulsos = [inicio(:), fim(:)];
end


function vetor_binario_duplo = aplicarLimiarDuplo(sinal_suavizado, fs, threshold_low, threshold_high, minDurationMs)
% APLICARLIMIARDUPLO Aplica detecção de pulsos usando dois limiares e duração mínima.
%
% vetor_binario_duplo = aplicarLimiarDuplo(sinal_suavizado, fs, ...
%                                          threshold_low, threshold_high, ...
%                                          minDurationMs)
%
% Inputs:
%   - sinal_suavizado  : vetor de dados do sinal já retificado/suavizado (1D).
%   - fs               : frequência de amostragem em Hz.
%   - threshold_low    : limiar de amplitude inferior.
%   - threshold_high   : limiar de amplitude superior.
%   - minDurationMs    : duração mínima (em milissegundos) para que um pulso
%                        seja considerado válido.
%
% Output:
%   - vetor_binario_duplo : vetor binário (0/1) indicando onde o sinal 
%                           excedeu o limiar duplo por tempo suficiente.

    N = length(sinal_suavizado);
    vetor_binario_duplo = zeros(1, N);
    
    % Converte milissegundos em número de amostras
    minDurationSamples = round(minDurationMs * fs / 1000);
    
    i = 1;
    while i <= N
        % Se passou do limiar alto, começamos uma possível detecção
        if sinal_suavizado(i) > threshold_high
            startIdx = i;
            j = i;
            
            % Continua enquanto o sinal estiver acima do limiar baixo
            while j <= N && (sinal_suavizado(j) > threshold_low)
                j = j + 1;
            end
            
            endIdx = j - 1;  % saiu porque sinal_suavizado(j) <= threshold_low ou j > N
            
            % Verifica a duração: se for maior que minDurationSamples, validamos
            if (endIdx - startIdx + 1) >= minDurationSamples
                vetor_binario_duplo(startIdx:endIdx) = 1;
            end
            
            % Pula direto pro fim do pulso
            i = endIdx + 1;
        else
            i = i + 1;
        end
    end
end


function pulso_dilatado = filtroDilatacao(pulso_binario, windowSize)
% FILTRODILATACAO Aplica dilatação morfológica a um sinal binário.
%
%   pulso_dilatado = filtroDilatacao(pulso_binario, windowSize)
%
% Inputs:
%   pulso_binario : Vetor binário (0 ou 1) representando os pulsos detectados.
%   windowSize    : Tamanho da janela (em número de amostras) para a dilatação.
%
% Output:
%   pulso_dilatado : Sinal binário após dilatação, onde pulsos próximos foram unidos.
%
% Exemplo de uso:
%   % Suponha que 'sinal_binario' seja o vetor obtido pela detecção dos pulsos.
%   windowSize = 50; % Tamanho da janela de 50 amostras (ajuste conforme necessário)
%   sinal_dilatado = filtroDilatacao(sinal_binario, windowSize);
%

    % Verifica se a entrada é um vetor
    if ~isvector(pulso_binario)
        error('O sinal de entrada deve ser um vetor.');
    end

    % Converte o sinal para vetor coluna para facilitar o processamento
    pulso_binario = pulso_binario(:);

    % Cria o elemento estruturante: vetor de 'ones' com tamanho especificado.
    elementoEstruturante = ones(windowSize, 1);

    % Aplica a dilatação pela operação de convolução:
    % Quando a soma na vizinhança (definida pelo elemento estruturante) for >0,
    % significa que há pelo menos um '1' na região e o output será 1.
    sinal_dilatado = conv(double(pulso_binario), elementoEstruturante, 'same');

    % Binariza o sinal: atribui 1 se a soma for maior que zero, ou zero caso contrário.
    pulso_dilatado = sinal_dilatado > 0;

    % Caso deseje o output como double (0 e 1) em vez de lógico, descomente a linha abaixo:
    % pulso_dilatado = double(pulso_dilatado);
end
