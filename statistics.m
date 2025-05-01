%% Estatística – Boxplots Diários + Controle
clear; clc; close all;

%% Parâmetros
base_path   = "C:\_PESQUISA\JP\Coleta";
num_coletas = 28;       
group_size  = 4;        
num_groups  = num_coletas / group_size;

%% Pré-alocação
dyn_ombro    = cell(num_coletas,1);
dyn_biceps   = cell(num_coletas,1);
stat_ombro   = nan(num_coletas,1);
stat_biceps  = nan(num_coletas,1);

%% Carregar dados de intervenção
for i = 1:num_coletas
    id   = sprintf('%02d', i);
    file = fullfile(base_path, ['Coleta-' id], ['resultados_emg_' id '.mat']);
    if exist(file,'file')
        load(file, ...
            'rms_dinamico_ombro','rms_dinamico_biceps', ...
            'RMS_estatico_ombro','RMS_estatico_biceps');
        dyn_ombro{i}   = rms_dinamico_ombro;
        dyn_biceps{i}  = rms_dinamico_biceps;
        stat_ombro(i)  = RMS_estatico_ombro;
        stat_biceps(i) = RMS_estatico_biceps;
    end
end

%% Carregar dados de controle
load(fullfile(base_path, 'controle.mat'), ...
     'control_dyn_ombro', 'control_dyn_biceps', ...
     'control_stat_ombro', 'control_stat_biceps');

%% 1) Boxplots Diários RMS Dinâmico + Controle
% Ombro
all_o = control_dyn_ombro(:)'; 
grp_o = repmat({'C'}, 1, numel(control_dyn_ombro));
for i=1:num_coletas
    d = dyn_ombro{i};
    if ~isempty(d)
        all_o = [all_o, d];
        grp_o = [grp_o, repmat({sprintf('%d',i)}, 1, numel(d))];
    end
end
figure;
boxplot(all_o, grp_o, 'Colors', 'b', 'Symbol', 'r+');
xlabel('Dia'); ylabel('RMS Dinâmico – Ombro (%MCV)');
title('Boxplot Diário – Dinâmico (Ombro) + Controle');
grid on;

% Bíceps
all_b = control_dyn_biceps(:)'; 
grp_b = repmat({'C'}, 1, numel(control_dyn_biceps));
for i=1:num_coletas
    d = dyn_biceps{i};
    if ~isempty(d)
        all_b = [all_b, d];
        grp_b = [grp_b, repmat({sprintf('%d',i)}, 1, numel(d))];
    end
end
figure;
boxplot(all_b, grp_b, 'Colors', 'r', 'Symbol', 'b+');
xlabel('Dia'); ylabel('RMS Dinâmico – Bíceps (%MCV)');
title('Boxplot Diário – Dinâmico (Bíceps) + Controle');
grid on;

%% 2) Médias Diárias RMS Dinâmico
daily_o = nan(num_coletas,1);
daily_b = nan(num_coletas,1);
for i=1:num_coletas
    if ~isempty(dyn_ombro{i})
        daily_o(i) = mean(dyn_ombro{i});
    end
    if ~isempty(dyn_biceps{i})
        daily_b(i) = mean(dyn_biceps{i});
    end
end

% Controle (média e SEM)
mean_control_dyn_ombro  = mean(control_dyn_ombro, 'omitnan');
sem_control_dyn_ombro   = std(control_dyn_ombro, 'omitnan') / sqrt(numel(control_dyn_ombro));
mean_control_dyn_biceps = mean(control_dyn_biceps, 'omitnan');
sem_control_dyn_biceps  = std(control_dyn_biceps, 'omitnan') / sqrt(numel(control_dyn_biceps));
mean_control_stat_ombro = mean(control_stat_ombro, 'omitnan');
mean_control_stat_biceps= mean(control_stat_biceps, 'omitnan');

%% 3) Médias Semanais e SEM
weeks = (1:num_groups)';
mean_dyn_o_week = nan(num_groups,1);
mean_dyn_b_week = nan(num_groups,1);
sem_dyn_o_week  = nan(num_groups,1);
sem_dyn_b_week  = nan(num_groups,1);
mean_stat_o_week= nan(num_groups,1);
mean_stat_b_week= nan(num_groups,1);

for g = 1:num_groups
    idx = (g-1)*group_size + (1:group_size);
    mean_dyn_o_week(g) = mean(daily_o(idx), 'omitnan');
    mean_dyn_b_week(g) = mean(daily_b(idx), 'omitnan');
    mean_stat_o_week(g)= mean(stat_ombro(idx), 'omitnan');
    mean_stat_b_week(g)= mean(stat_biceps(idx), 'omitnan');
    sem_dyn_o_week(g)  = std(daily_o(idx), 'omitnan')/sqrt(sum(~isnan(daily_o(idx))));
    sem_dyn_b_week(g)  = std(daily_b(idx), 'omitnan')/sqrt(sum(~isnan(daily_b(idx))));
end

%% 4) Plot Linha – Médias Semanais DINÂMICO Ombro e Bíceps (com Controle)
figure;
hold on;
yline(mean_control_dyn_ombro, 'k-', 'LineWidth', 1.5);
yline(mean_control_dyn_ombro + sem_control_dyn_ombro, 'k--', 'LineWidth', 1);
yline(mean_control_dyn_ombro - sem_control_dyn_ombro, 'k--', 'LineWidth', 1);
errorbar(weeks, mean_dyn_o_week, sem_dyn_o_week, 'bo-', 'MarkerFaceColor', 'b', 'LineWidth',2);
xlabel('Semana');
ylabel('Média RMS Dinâmico – Ombro (%MCV)');
title('Evolução Semanal – Dinâmico Ombro');
legend({'Média Controle','+ SEM','- SEM','Semanas'}, 'Location', 'Best');
grid on;
xlim([0.5 num_groups+0.5]);
xticks(weeks);
ylim([0 0.3])
hold off;

figure;
hold on;
yline(mean_control_dyn_biceps, 'k-', 'LineWidth', 1.5);
yline(mean_control_dyn_biceps + sem_control_dyn_biceps, 'k--', 'LineWidth', 1);
yline(mean_control_dyn_biceps - sem_control_dyn_biceps, 'k--', 'LineWidth', 1);
errorbar(weeks, mean_dyn_b_week, sem_dyn_b_week, 'rs-', 'MarkerFaceColor', 'r', 'LineWidth',2);
xlabel('Semana');
ylabel('Média RMS Dinâmico – Bíceps (%MCV)');
title('Evolução Semanal – Dinâmico Bíceps');
legend({'Média Controle','+ SEM','- SEM','Semanas'}, 'Location', 'Best');
grid on;
xlim([0.5 num_groups+0.5]);
xticks(weeks);
ylim([0 0.3])
hold off;

%% 5) Plot Linha – Médias Semanais ESTÁTICO Ombro e Bíceps (sem SEM, linha direta)
figure;
plot(weeks, mean_stat_o_week, 'b-^', 'MarkerFaceColor', 'b', 'LineWidth',2);
hold on;
yline(mean_control_stat_ombro, 'k-', 'LineWidth', 1.5);
xlabel('Semana');
ylabel('Média RMS Estático – Ombro (%MCV)');
title('Evolução Semanal – Estático Ombro');
legend({'Semanas','Média Controle'}, 'Location', 'Best');
grid on;
xlim([0.5 num_groups+0.5]);
xticks(weeks);
ylim([0 0.25])
hold off;

figure;
plot(weeks, mean_stat_b_week, 'r-^', 'MarkerFaceColor', 'r', 'LineWidth',2);
hold on;
yline(mean_control_stat_biceps, 'k-', 'LineWidth', 1.5);
xlabel('Semana');
ylabel('Média RMS Estático – Bíceps (%MCV)');
title('Evolução Semanal – Estático Bíceps');
legend({'Semanas','Média Controle'}, 'Location', 'Best');
grid on;
xlim([0.5 num_groups+0.5]);
xticks(weeks);
ylim([0 0.25])
hold off;

%% 6) Teste Mann–Whitney – Dinâmico e Estático

fprintf('\n---- Comparação Controle vs Semanas (Dinâmico e Estático) ----\n');

pval_dyn_ombro = nan(num_groups,1);
pval_dyn_biceps = nan(num_groups,1);
pval_stat_ombro = nan(num_groups,1);
pval_stat_biceps = nan(num_groups,1);

for w = 1:num_groups
    idx = (w-1)*group_size + (1:group_size);
    semana_ombro_dyn = daily_o(idx);
    semana_biceps_dyn = daily_b(idx);
    semana_ombro_stat = stat_ombro(idx);
    semana_biceps_stat = stat_biceps(idx);

    % Teste Mann–Whitney para cada condição
    pval_dyn_ombro(w) = ranksum(semana_ombro_dyn, control_dyn_ombro);
    pval_dyn_biceps(w) = ranksum(semana_biceps_dyn, control_dyn_biceps);
    pval_stat_ombro(w) = ranksum(semana_ombro_stat, control_stat_ombro);
    pval_stat_biceps(w) = ranksum(semana_biceps_stat, control_stat_biceps);

    % Imprimir
    fprintf('Semana %d:\n', w);
    fprintf(' - Dinâmico Ombro: p = %.4f\n', pval_dyn_ombro(w));
    fprintf(' - Dinâmico Bíceps: p = %.4f\n', pval_dyn_biceps(w));
    fprintf(' - Estático Ombro: p = %.4f\n', pval_stat_ombro(w));
    fprintf(' - Estático Bíceps: p = %.4f\n\n', pval_stat_biceps(w));
end
fprintf('---------------------------------------------------------------\n');

%% 7) Exportar para Tabela Excel
T_controle_vs_semanas = table(...
    weeks, ...
    mean_dyn_o_week, mean_dyn_b_week, ...
    mean_stat_o_week, mean_stat_b_week, ...
    pval_dyn_ombro, pval_dyn_biceps, ...
    pval_stat_ombro, pval_stat_biceps, ...
    'VariableNames', {'Semana', ...
                      'RMS_Dinamico_Ombro', 'RMS_Dinamico_Biceps', ...
                      'RMS_Estatico_Ombro', 'RMS_Estatico_Biceps', ...
                      'p_Dinamico_Ombro', 'p_Dinamico_Biceps', ...
                      'p_Estatico_Ombro', 'p_Estatico_Biceps'});

writetable(T_controle_vs_semanas, fullfile(base_path, 'controle_vs_semanas_completo.xlsx'));

fprintf('\nTabela "controle_vs_semanas_completo.xlsx" salva.\n');
