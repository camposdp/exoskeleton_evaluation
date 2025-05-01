%% ergo_target_analysis.m – Versão Corrigida para 5 Quantis e Shading Ajustado com Controle
%clear; clc; close all;

%% Parâmetros
base_path   = "C:\_PESQUISA\JP\Coleta";
num_coletas = 28;
group_size  = 4;
num_weeks   = num_coletas / group_size;   % = 6 semanas
quantiles   = [5,25,50,75,95] / 100;       % 5 quantis (fração)

% Curvas de referência (do artigo)
x1 = [2,10,50]/100;
x2 = [5,14,70]/100;
y_ref = [5,50,95]/100;

%% Função auxiliar: ajusta curva logística
tmp = sprintf('');
function [xf, yf] = fitLogistic(xq, yq)
    ytf   = log((1 - yq) ./ yq);
    A     = [xq; ones(1,numel(xq))]';
    theta = (A' * A) \ (A' * ytf');
    xf    = linspace(0,1,200)';
    yf    = 1 ./ (1 + exp([xf, ones(size(xf))] * theta));
end

%% Pré-alocação
drms_o = nan(num_coletas,1);
drms_b = nan(num_coletas,1);
p_o    = nan(num_coletas,5);
p_b    = nan(num_coletas,5);

%% Carregar dados da intervenção
for i = 1:num_coletas
    id  = sprintf('%02d', i);
    fn  = fullfile(base_path, ['Coleta-' id], ['resultados_emg_' id '.mat']);
    S   = load(fn,'rms_dinamico_ombro','rms_dinamico_biceps','p_dyn_ombro','p_dyn_biceps');
    drms_o(i) = mean(S.rms_dinamico_ombro);
    drms_b(i) = mean(S.rms_dinamico_biceps);
    p_o(i,:)  = S.p_dyn_ombro*100;
    p_b(i,:)  = S.p_dyn_biceps*100;
end

% Carregar dados de controle
load(fullfile(base_path, 'controle.mat'), 'control_p_dyn_ombro', 'control_p_dyn_biceps');

% Média dos percentis de controle
mean_control_p_ombro  = mean(control_p_dyn_ombro, 1) * 100;
mean_control_p_biceps = mean(control_p_dyn_biceps, 1) * 100;

% Ajuste das curvas de referência
[xf1, yf1] = fitLogistic(x1, y_ref);
[xf2, yf2] = fitLogistic(x2, y_ref);
xf_plot = xf1 * 100;
y1_plot = yf1;
y2_plot = yf2;

%% Encontrar pior e melhor dias
[~, worst_o] = max(drms_o);
[~, best_o ] = min(drms_o);
[~, worst_b] = max(drms_b);
[~, best_b ] = min(drms_b);

%% Função para shading
function shadeRegion(xp, y_top, y_bot, color)
    fill([xp' fliplr(xp')], [y_top' fliplr(y_bot')], color, ...
         'FaceAlpha',0.3,'EdgeColor','none');
end

%% Plot Ergo-Target – Ombro Pior vs Melhor (com Controle)
figure; hold on;
shadeRegion(xf_plot, y2_plot, zeros(size(y2_plot)), [1 0 0]);
shadeRegion(xf_plot, y2_plot, y1_plot, [0.9 1 0]);
shadeRegion(xf_plot, ones(size(y1_plot)), y1_plot, [0.6 1 0.6]);

% Controle
plot(mean_control_p_ombro, quantiles, '-k^', 'MarkerFaceColor','k','LineWidth',1.5);

% Pior e melhor
plot(p_o(worst_o,:), quantiles, '-ro','MarkerFaceColor','r','MarkerSize',8,'LineWidth',1.5);
plot(p_o(best_o ,:), quantiles, '-bs','MarkerFaceColor','b','MarkerSize',8,'LineWidth',1.5);

yticklabels(0:10:100)
xlim([0 60])
xlabel('MCV (%)'); ylabel('Percentil');
title('Ergo-Target Dinâmico – Ombro (Pior vs Melhor)');
legend({'Risco Alto','Risco Médio','Risco Baixo','Controle','Pior Dia','Melhor Dia'},'Location','east');
grid on; hold off;

%% Plot Ergo-Target – Bíceps Pior vs Melhor (com Controle)
figure; hold on;
shadeRegion(xf_plot, y2_plot, zeros(size(y2_plot)), [1 0 0]);
shadeRegion(xf_plot, y2_plot, y1_plot, [0.9 1 0]);
shadeRegion(xf_plot, ones(size(y1_plot)), y1_plot, [0.6 1 0.6]);

% Controle
plot(mean_control_p_biceps, quantiles, '-k^', 'MarkerFaceColor','k','LineWidth',1.5);

% Pior e melhor
plot(p_b(worst_b,:), quantiles, '-ro','MarkerFaceColor','r','MarkerSize',8,'LineWidth',1.5);
plot(p_b(best_b ,:), quantiles, '-bs','MarkerFaceColor','b','MarkerSize',8,'LineWidth',1.5);

yticklabels(0:10:100)
xlim([0 60])
xlabel('MCV (%)'); ylabel('Percentil');
title('Ergo-Target Dinâmico – Bíceps (Pior vs Melhor)');
legend({'Risco Alto','Risco Médio','Risco Baixo','Controle','Pior Dia','Melhor Dia'},'Location','east');
grid on; hold off;

%% Medianas Semanais (5 quantis)
data_o = reshape(p_o, group_size, [], 5);
med_o = squeeze(median(data_o,1));

data_b = reshape(p_b, group_size, [], 5);
med_b = squeeze(median(data_b,1));

