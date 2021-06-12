clear all;

path_to_shap = 'E:/YandexDisk/Work/dnamvae/data/datasets/combo/GSE40279_GSE87571_EPIC_GSE55763/vt_score_more_0.005_none/interpretation';
path_to_ann = 'E:/YandexDisk/Work/dnamvae/data/annotation';

ann_type = 'full';
fn = sprintf('%s/ann_%s.xlsx', path_to_ann, ann_type);
ann = readtable(fn, 'ReadRowNames', true);

num_features = 100;
num_subjects = 300;

num_features_to_plot = 10;
opacity = 0.65;

fn = sprintf('%s/shap_values_%d_%d.xlsx', path_to_shap, num_subjects, num_features);
opts = detectImportOptions(fn);
tbl = readtable(fn, opts);

features = tbl.Properties.VariableNames';

features_beta = features(contains(features, '_beta'));
features_shap = features(contains(features, '_shap'));

num_cpgs = size(features_beta, 1);

x_positions = linspace(1, 2 * num_features_to_plot - 1, num_features_to_plot)';
fig = figure;
propertyeditor('on');
cpgs = {};
genes = {};
mean_abs_shaps = zeros(num_features_to_plot, 1);
for cpg_id = 1:num_features_to_plot
    
    cg_split = split(features_beta(cpg_id), '_');
    cg = cg_split{1};
    cpgs{cpg_id} = cg;
    gene_raw = string(ann{cg, 'Gene'});
    gene_array = split(gene_raw, ';')';
    gene_unique = unique(gene_array);
    gene = join(gene_unique, ';');
    genes{cpg_id} = gene{1};
    cg_beta = tbl{:, features_beta(cpg_id)};
    cg_beta_scale = rescale(cg_beta);
    cg_shap = tbl{:, features_shap(cpg_id)};
    
    x = ones(size(cg_shap, 1), 1) * x_positions(cpg_id);
    y = cg_shap;
    c = cg_beta_scale;
    
    s = swarmchart(x, y, 20, c, 'filled');
    hold all;
    
    mean_abs_shaps(cpg_id) = mean(abs(cg_shap));
end

labels = [cpgs; genes];
tickLabels = strtrim(sprintf('%s\\newline%s\n', labels{:}));
set(gca,'XTick', x_positions, 'XTickLabel', tickLabels)
xtickangle(90);
axis auto;
xlim([min(x_positions) - 1, max(x_positions) + 1])
box on;
grid on;
set(gca, 'FontSize', 20);
xlabel('', 'Interpreter', 'latex');
ylabel('SHAP value (impact on model output) ', 'Interpreter', 'latex');
ax = gca;
ax.XAxis.FontSize = 14;
h = colorbar;
h.Ticks = [0, 1];
h.TickLabels = {'Min', 'Max'};
title(h, 'Feature', 'FontSize', 20, 'interpreter','latex');
set(h, 'TickLabelInterpreter', 'latex');
fn_fig = sprintf('%s/shap_dataset_swarmchart_%d_%d', path_to_shap, num_subjects, num_features);
oqs_save_fig(fig, fn_fig)

fig = figure;
propertyeditor('on');
labels = [flip(cpgs); flip(genes)];
tickLabels = strtrim(sprintf('%s\\newline%s\n', labels{:}));
set(gca,'XTick', x_positions, 'XTickLabel', tickLabels)
s = barh(flip(mean_abs_shaps), 'FaceColor', 'red');
set(gca, 'FontSize', 30);
box on;
grid on;
yticks(linspace(1, num_features_to_plot, num_features_to_plot))
ylim([0.5, num_features_to_plot + 0.5])
hold all;
set(gca, 'yTickLabel', tickLabels);
ax = gca;
ax.YAxis.FontSize = 20;
xlabel('mean($|$SHAP value$|$)', 'Interpreter', 'latex');
ylabel('', 'Interpreter', 'latex');
fn_fig = sprintf('%s/shap_dataset_barh_%d_%d', path_to_shap, num_subjects, num_features);
oqs_save_fig(fig, fn_fig)


