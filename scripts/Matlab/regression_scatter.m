clear all;

path = 'E:/YandexDisk/Work/dnamvae/models/fcmlp_model/logs/runs/2021-06-09/01-45-07';
fn = sprintf('%s/inference.xlsx', path);

color = 'red';
opacity = 0.65;

opts = detectImportOptions(fn);
tbl = readtable(fn, opts);

x_var = 'age_real';
x_label = 'age_real';
xlims = [0; 100];
y_var = 'age_pred';
y_label = 'age_pred';
ylims = [0; 100];

xs = tbl.(x_var);
ys = tbl.(y_var);

fig1 = figure;
propertyeditor('on');
grid on;

h = scatter(xs, ys, 250, 'o', 'LineWidth',  1, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', color, 'MarkerEdgeAlpha', opacity, 'MarkerFaceAlpha', opacity);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold all;

T = table(xs, ys, 'VariableNames', {x_var, y_var});
lm = fitlm(T, sprintf('%s~%s', y_var, x_var));
R2 = lm.Rsquared.Ordinary
RMSE = lm.RMSE
ae = abs(x_var - y_var);
MAE = mean(ae)

x_fit = [xlims(1); xlims(2)];
y_fit = lm.Coefficients{'(Intercept)','Estimate'} + x_fit * lm.Coefficients{x_var,'Estimate'};
h = plot(x_fit, y_fit, 'LineWidth', 2, 'Color', color);
h.Annotation.LegendInformation.IconDisplayStyle = 'off';

set(gca, 'FontSize', 40);
xlabel(strrep(x_label,'_','\_'), 'Interpreter', 'latex');
set(gca, 'FontSize', 40);
ylabel(strrep(y_label,'_','\_'), 'Interpreter', 'latex');
ax = gca;
set(ax,'TickLabelInterpreter','Latex')

%legend(gca,'off');
%legend('Location', legend_location, 'NumColumns', 1, 'Interpreter', 'latex');

box on;
xlim(xlims);
ylim(ylims);
fn_fig = sprintf('%s/inference', path);
oqs_save_fig(fig1, fn_fig)

