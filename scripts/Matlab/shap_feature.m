clear all;

path_to_shap = 'E:/YandexDisk/Work/dnamvae/data/datasets/combo/GSE40279_GSE87571_EPIC_GSE55763/vt_score_more_0.005_none/interpretation';
save_path = sprintf('%s/feature', path_to_shap);
if not(isfolder(save_path))
    mkdir(save_path)
end
path_to_ann = 'E:/YandexDisk/Work/dnamvae/data/annotation';

ann_type = 'full';
fn = sprintf('%s/ann_%s.xlsx', path_to_ann, ann_type);
ann = readtable(fn, 'ReadRowNames', true);

num_features = 100;
num_subjects = 1000;
outcome = 'age';

num_features_to_plot = 10;
opacity = 0.65;

fn = sprintf('%s/shap_values_%d_%d.xlsx', path_to_shap, num_subjects, num_features);
opts = detectImportOptions(fn);
tbl = readtable(fn, opts);

features = tbl.Properties.VariableNames';

features_beta = features(contains(features, '_beta'));
features_shap = features(contains(features, '_shap'));

num_cpgs = size(features_beta, 1);

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
    out_real = tbl{:, 'preds'};
    
    x = cg_beta;
    y = cg_shap;
    c = out_real;
    
    fig = figure;
    propertyeditor('on');
    h = scatter(x, y, 100, c, 'filled', 'MarkerEdgeAlpha', opacity, 'MarkerFaceAlpha', opacity);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    box on;
    grid on;
    set(gca, 'FontSize', 30);
    xlabel({cpgs{cpg_id}, genes{cpg_id}}, 'Interpreter', 'latex');
    ylabel('SHAP value', 'Interpreter', 'latex');
    ax = gca;
    set(ax,'TickLabelInterpreter','Latex')
    
    h = colorbar;
    title(h, 'Age', 'FontSize', 20, 'interpreter','latex');
    set(h, 'TickLabelInterpreter', 'latex');
    
    fn_fig = sprintf('%s/shap_feature_scatter_%d_%s_%d_%d', save_path, cpg_id, cpgs{cpg_id}, num_subjects, num_features);
    oqs_save_fig(fig, fn_fig)
    
    
    fig = figure;
    propertyeditor('on');
    h = scatter(out_real, x, 100, y, 'filled', 'MarkerEdgeAlpha', opacity, 'MarkerFaceAlpha', opacity);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    box on;
    grid on;
    set(gca, 'FontSize', 30);
    xlabel(outcome, 'Interpreter', 'latex');
    ylabel({cpgs{cpg_id}, genes{cpg_id}}, 'Interpreter', 'latex');
    ax = gca;
    set(ax,'TickLabelInterpreter','Latex')
    
    h = colorbar;
    title(h, 'SHAP', 'FontSize', 20, 'interpreter','latex');
    set(h, 'TickLabelInterpreter', 'latex');
    
    fn_fig = sprintf('%s/outcome_scatter_%d_%s_%d_%d', save_path, cpg_id, cpgs{cpg_id}, num_subjects, num_features);
    oqs_save_fig(fig, fn_fig)
end
