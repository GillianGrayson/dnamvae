clear all;

path_to_shap = 'E:/YandexDisk/Work/dnamvae/data/datasets/combo/GSE40279_GSE87571_EPIC_GSE55763/vt_score_more_0.005_none/interpretation';
save_path = sprintf('%s/subject', path_to_shap);
if not(isfolder(save_path))
    mkdir(save_path)
end
path_to_ann = 'E:/YandexDisk/Work/dnamvae/data/annotation';

ann_type = 'full';
fn = sprintf('%s/ann_%s.xlsx', path_to_ann, ann_type);
ann = readtable(fn, 'ReadRowNames', true);

num_features = 100;
num_subjects = 1000;

num_features_to_plot = 10;
num_subjects_to_plot = 100;

fn = sprintf('%s/shap_values_%d_%d.xlsx', path_to_shap, num_subjects, num_features);
opts = detectImportOptions(fn);
tbl = readtable(fn, opts);

features = tbl.Properties.VariableNames';

features_beta = features(contains(features, '_beta'));
features_shap = features(contains(features, '_shap'));

num_cpgs = size(features_beta, 1);

x_positions = linspace(1, 2 * num_features_to_plot - 1, num_features_to_plot)';

for s_id = 1:num_subjects_to_plot
    
    cpgs = {};
    genes = {};
    shaps = zeros(num_features_to_plot, 1);
    for cpg_id = 1:num_features_to_plot
        
        cg_split = split(features_beta(cpg_id), '_');
        cg = cg_split{1};
        cpgs{cpg_id} = cg;
        gene_raw = string(ann{cg, 'Gene'});
        gene_array = split(gene_raw, ';')';
        gene_unique = unique(gene_array);
        gene = join(gene_unique, ';');
        genes{cpg_id} = gene{1};
        cg_shap = tbl{s_id, features_shap(cpg_id)};
        
        shaps(cpg_id) = cg_shap;
    end
    
    fig = figure;
    propertyeditor('on');
    labels = [flip(cpgs); flip(genes)];
    tickLabels = strtrim(sprintf('%s\\newline%s\n', labels{:}));
    set(gca,'XTick', x_positions, 'XTickLabel', tickLabels)
    s = barh(flip(shaps), 'FaceColor', 'red');
    set(gca, 'FontSize', 30);
    box on;
    grid on;
    yticks(linspace(1, num_features_to_plot, num_features_to_plot))
    ylim([0.5, num_features_to_plot + 0.5])
    hold all;
    set(gca, 'yTickLabel', tickLabels);
    ax = gca;
    ax.YAxis.FontSize = 20;
    xlabel('SHAP value', 'Interpreter', 'latex');
    ylabel('', 'Interpreter', 'latex');
    title(sprintf('$Age = %0.2f$', tbl{s_id, 'preds'}), 'interpreter', 'latex')
    fn_fig = sprintf('%s/shap_subject_barh_%d_%d_%d', save_path, s_id, num_subjects, num_features);
    oqs_save_fig(fig, fn_fig)

end



