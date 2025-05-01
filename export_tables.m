%% exporta_tabelas_emg.m – Versão Corrigida e Organizada
clear; clc;

%% CONFIGURAÇÕES BÁSICAS
BASE_PATH   = "C:\_PESQUISA\JP\Coleta";   % ajuste se necessário
RESULTS_PATH = fullfile("C:\_PESQUISA\JP\github", 'results');
if ~exist(RESULTS_PATH, 'dir')
    mkdir(RESULTS_PATH);
end

N_COLETAS   = 28;                        
GRP_SIZE    = 4;                         
N_SEMANAS   = N_COLETAS / GRP_SIZE;       

%% PRÉ‑ALOCAÇÃO
dyn_ombro  = cell(N_COLETAS,1);
dyn_biceps = cell(N_COLETAS,1);
stat_ombro = nan(N_COLETAS,1);
stat_bic   = nan(N_COLETAS,1);

p_dyn_ombro = nan(N_COLETAS,5);
p_dyn_bic   = nan(N_COLETAS,5);

%% CARREGAR DADOS DA INTERVENÇÃO
for k = 1:N_COLETAS
    id  = sprintf('%02d',k);
    f   = fullfile(BASE_PATH, ['Coleta-' id], ['resultados_emg_' id '.mat']);
    if exist(f,"file")
        S = load(f, ...
            "rms_dinamico_ombro","rms_dinamico_biceps", ...
            "RMS_estatico_ombro","RMS_estatico_biceps", ...
            "p_dyn_ombro","p_dyn_biceps");

        dyn_ombro{k}   = S.rms_dinamico_ombro(:);
        dyn_biceps{k}  = S.rms_dinamico_biceps(:);
        stat_ombro(k)  = S.RMS_estatico_ombro;
        stat_bic(k)    = S.RMS_estatico_biceps;
        p_dyn_ombro(k,:) = S.p_dyn_ombro;
        p_dyn_bic(k,:)   = S.p_dyn_biceps;
    else
        warning("Arquivo não encontrado: %s",f);
    end
end

%% CARREGAR DADOS DE CONTROLE
load(fullfile(BASE_PATH, 'controle.mat'), ...
     'control_dyn_ombro','control_dyn_biceps', ...
     'control_stat_ombro','control_stat_biceps', ...
     'mean_control_p_dyn_ombro','mean_control_p_dyn_biceps');

%% MÉDIAS DIÁRIAS DOS RMS DINÂMICOS
daily_ombro = nan(N_COLETAS,1);
daily_biceps = nan(N_COLETAS,1);

for k = 1:N_COLETAS
    if ~isempty(dyn_ombro{k})
        daily_ombro(k) = mean(dyn_ombro{k});
    end
    if ~isempty(dyn_biceps{k})
        daily_biceps(k) = mean(dyn_biceps{k});
    end
end

%% CÁLCULO DE MÉDIAS E SEM
weeks = (1:N_SEMANAS)';
mean_o_week = nan(N_SEMANAS,1);   sem_o_week = nan(N_SEMANAS,1);
mean_b_week = nan(N_SEMANAS,1);   sem_b_week = nan(N_SEMANAS,1);
mean_stat_o = nan(N_SEMANAS,1);   sem_stat_o = nan(N_SEMANAS,1);
mean_stat_b = nan(N_SEMANAS,1);   sem_stat_b = nan(N_SEMANAS,1);

for s = 1:N_SEMANAS
    idx = (s-1)*GRP_SIZE + (1:GRP_SIZE);
    mean_o_week(s) = mean(daily_ombro(idx),'omitnan');
    mean_b_week(s) = mean(daily_biceps(idx),'omitnan');
    sem_o_week(s)  = std(daily_ombro(idx),'omitnan')/sqrt(sum(~isnan(daily_ombro(idx))));
    sem_b_week(s)  = std(daily_biceps(idx),'omitnan')/sqrt(sum(~isnan(daily_biceps(idx))));
    mean_stat_o(s) = mean(stat_ombro(idx),'omitnan');
    sem_stat_o(s)  = std(stat_ombro(idx),'omitnan')/sqrt(sum(~isnan(stat_ombro(idx))));
    mean_stat_b(s) = mean(stat_bic(idx),'omitnan');
    sem_stat_b(s)  = std(stat_bic(idx),'omitnan')/sqrt(sum(~isnan(stat_bic(idx))));
end

%% CÁLCULO DE DIFERENÇA % RELATIVA AO CONTROLE
mean_control_dyn_ombro  = mean(control_dyn_ombro, 'omitnan');
mean_control_dyn_biceps = mean(control_dyn_biceps, 'omitnan');
mean_control_stat_ombro = mean(control_stat_ombro, 'omitnan');
mean_control_stat_biceps= mean(control_stat_biceps, 'omitnan');

diff_dyn_ombro  = 100 * (mean_o_week - mean_control_dyn_ombro) / mean_control_dyn_ombro;
diff_dyn_biceps = 100 * (mean_b_week - mean_control_dyn_biceps) / mean_control_dyn_biceps;
diff_stat_ombro = 100 * (mean_stat_o - mean_control_stat_ombro) / mean_control_stat_ombro;
diff_stat_biceps= 100 * (mean_stat_b - mean_control_stat_biceps) / mean_control_stat_biceps;

%% CÁLCULO DOS PERCENTIS (medianas semanais)
data_o = reshape(p_dyn_ombro, GRP_SIZE, [], 5);   % dias × semanas × quantis
data_b = reshape(p_dyn_bic  , GRP_SIZE, [], 5);
med_o  = squeeze(median(data_o,1));
med_b  = squeeze(median(data_b,1));

%% ANOVA ENTRE SEMANAS
labels = repelem(weeks, GRP_SIZE)';

[p_dyn_o, ~, stats_dyn_o] = anova1(daily_ombro, labels, 'off');
[p_dyn_b, ~, stats_dyn_b] = anova1(daily_biceps, labels, 'off');
[p_stat_o, ~, stats_stat_o] = anova1(stat_ombro, labels, 'off');
[p_stat_b, ~, stats_stat_b] = anova1(stat_bic, labels, 'off');

sig_dyn_o = []; sig_dyn_b = []; sig_stat_o = []; sig_stat_b = [];
if p_dyn_o < 0.05, sig_dyn_o = multcompare(stats_dyn_o, 'Display', 'off'); sig_dyn_o = sig_dyn_o(sig_dyn_o(:,6)<0.05,:); end
if p_dyn_b < 0.05, sig_dyn_b = multcompare(stats_dyn_b, 'Display', 'off'); sig_dyn_b = sig_dyn_b(sig_dyn_b(:,6)<0.05,:); end
if p_stat_o < 0.05, sig_stat_o = multcompare(stats_stat_o, 'Display', 'off'); sig_stat_o = sig_stat_o(sig_stat_o(:,6)<0.05,:); end
if p_stat_b < 0.05, sig_stat_b = multcompare(stats_stat_b, 'Display', 'off'); sig_stat_b = sig_stat_b(sig_stat_b(:,6)<0.05,:); end

%% CRIAÇÃO DAS TABELAS

% --- Tabela 1: RMS ± SEM + Controle + Diferença %
T1 = table(weeks, ...
           mean_o_week, sem_o_week, repmat(mean_control_dyn_ombro,N_SEMANAS,1), diff_dyn_ombro, ...
           mean_b_week, sem_b_week, repmat(mean_control_dyn_biceps,N_SEMANAS,1), diff_dyn_biceps, ...
           mean_stat_o, sem_stat_o, repmat(mean_control_stat_ombro,N_SEMANAS,1), diff_stat_ombro, ...
           mean_stat_b, sem_stat_b, repmat(mean_control_stat_biceps,N_SEMANAS,1), diff_stat_biceps, ...
    'VariableNames', {'Semana', ...
                      'RMS_Ombro_Din', 'SEM_Ombro_Din', 'Controle_Ombro_Din', 'DeltaPct_Ombro_Din', ...
                      'RMS_Biceps_Din','SEM_Biceps_Din','Controle_Biceps_Din','DeltaPct_Biceps_Din', ...
                      'RMS_Ombro_Est', 'SEM_Ombro_Est', 'Controle_Ombro_Est', 'DeltaPct_Ombro_Est', ...
                      'RMS_Biceps_Est','SEM_Biceps_Est','Controle_Biceps_Est','DeltaPct_Biceps_Est'});

% --- Tabela 2: ANOVA e Comparações significativas
T2 = table( ...
    {'Dinâmico Ombro'; 'Dinâmico Bíceps'; 'Estático Ombro'; 'Estático Bíceps'}, ...
    [p_dyn_o; p_dyn_b; p_stat_o; p_stat_b], ...
    {sig_dyn_o; sig_dyn_b; sig_stat_o; sig_stat_b}, ...
    'VariableNames', {'Variável', 'p_ANOVA', 'Pares_Significativos'});

% --- Tabela 3: Percentis Semanais Ombro
quantilLabels = [5; 25; 50; 75; 95];
T3 = array2table(med_o', 'VariableNames', {'Sem1','Sem2','Sem3','Sem4','Sem5','Sem6','Sem7'});
T3 = addvars(T3, quantilLabels, 'Before', 1, 'NewVariableNames', 'Quantil');

% --- Tabela 4: Percentis Semanais Bíceps
T4 = array2table(med_b', 'VariableNames', {'Sem1','Sem2','Sem3','Sem4','Sem5','Sem6','Sem7'});
T4 = addvars(T4, quantilLabels, 'Before', 1, 'NewVariableNames', 'Quantil');

%% EXPORTAÇÃO
writetable(T1, fullfile(RESULTS_PATH, 'Resultados_RMS_Semanais.csv'));
writetable(T2, fullfile(RESULTS_PATH, 'Resultados_ANOVA.csv'));
writetable(T3, fullfile(RESULTS_PATH, 'Percentis_Semanais_Ombro.csv'));
writetable(T4, fullfile(RESULTS_PATH, 'Percentis_Semanais_Biceps.csv'));

writetable(T1, fullfile(RESULTS_PATH, 'Relatorio_EMG.xlsx'), 'Sheet', 'RMS_Semanais');
writetable(T2, fullfile(RESULTS_PATH, 'Relatorio_EMG.xlsx'), 'Sheet', 'ANOVA');
writetable(T3, fullfile(RESULTS_PATH, 'Relatorio_EMG.xlsx'), 'Sheet', 'Percentis_Ombro');
writetable(T4, fullfile(RESULTS_PATH, 'Relatorio_EMG.xlsx'), 'Sheet', 'Percentis_Biceps');

disp('Todas as tabelas exportadas em "results" com sucesso.');
