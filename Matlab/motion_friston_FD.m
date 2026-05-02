%% --- Paramètres ---
rp_file = '/data/belleville/Arsenii/03-Raw_fMRI/sub-3002498/ses-t04/func/rp_12-fmri-MemoryTask_fmri-MemoryTask_20190527141834_12.txt';
[out_dir, rp_name, ~] = fileparts(rp_file);
FD_threshold = 0.5;             % seuil FD en mm
r = 50;                         % rayon du cerveau en mm pour rotations


%% --- Charger les paramètres de mouvement ---
rp = load(rp_file);  % Nx6 : [Tx Ty Tz Rx Ry Rz]
[nVol, nPar] = size(rp);
if nPar ~= 6
    error('La matrice rp doit être Nx6');
end
%% --- Calculer le déplacement total (Total Displacement) ---
% Isoler les paramètres de translation (x, y, z en mm)
translations = rp(:, 1:3);

% Soustraire la position du premier volume pour avoir un point de référence à 0
translations_rel = translations - translations(1, :);

% Calculer la distance euclidienne pour chaque volume depuis l'origine
distance_euclidienne = sqrt(sum(translations_rel.^2, 2));

% Trouver le déplacement total maximum
total_displacement = max(distance_euclidienne);

fprintf('Déplacement total maximum = %.3f mm\n', total_displacement);

% Vérifier l'exclusion
if total_displacement > 3
    fprintf('Attention : Le déplacement total dépasse 3 mm.\n');
end
%% --- Calculer les différences (delta entre volumes) ---
diff_rp = [zeros(1,6); diff(rp)];  % premier volume = 0

%% --- Convertir rotations en déplacement linéaire ---
diff_rp(:,4:6) = diff_rp(:,4:6) * r;

%% --- Calculer FD (somme des valeurs absolues) ---
FD = sum(abs(diff_rp), 2);  % Nx1

%% --- Identifier volumes à scrubbing ---
scrub_idx = find(FD > FD_threshold);

%% --- Résultats ---
fprintf('FD moyen = %.3f mm\n', mean(FD));
fprintf('Nombre de volumes à scrubbing (FD > %.2f mm) : %d / %d\n', ...
    FD_threshold, length(scrub_idx), nVol);


%% --- Optionnel : visualisation (Server-Safe) ---
% 1. Create the figure invisibly so it doesn't pop up
fig = figure('Visible', 'off');

% 2. Draw the plot as usual
plot(FD, 'b-o');
hold on;
yline(FD_threshold, 'r--', 'Seuil FD 0.5 mm');
xlabel('Volume');
ylabel('FD (mm)');
title('Framewise Displacement');
grid on;

% 3. Save the image to the hard drive (dynamically named)
img_filename = fullfile(out_dir, ['plot_FD_' rp_name '.png']);
saveas(fig, img_filename);

% 4. Close the invisible figure to free up the server's RAM
close(fig);

fprintf('Image FD sauvegardée dans : %s\n', img_filename);


%% 1) Paramètres originaux (6)
P = rp;

%% 2) Dérivées temporelles (6)
dP = [zeros(1,6); diff(P)];

%% 3) Carrés des paramètres (6)
P2 = P.^2;

%% 4) Carrés des dérivées (6)
dP2 = dP.^2;

%% 5) Concaténation -> 24 régressseurs
Friston24 = [P dP P2 dP2];

%% Vérification
disp(size(Friston24))  % doit afficher [N 24]




% 3. Create the dummy regressors (the "spike" columns)
nVol = size(Friston24, 1);
nSpikes = length(scrub_idx); % This should be 5
DummyMatrix = zeros(nVol, nSpikes);

for i = 1:nSpikes
    % Place a 1 at the exact volume that needs to be scrubbed
    DummyMatrix(scrub_idx(i), i) = 1; 
end

% 4. Combine them into one giant 29-column matrix
FinalRegressors = [Friston24, DummyMatrix];

% 5. Save the final file for SPM
output_name = fullfile(out_dir, ['nuisance_combined_' rp_name '.txt']);
save(output_name, 'FinalRegressors', '-ascii');

fprintf('Success! Saved %d columns to %s\n', size(FinalRegressors, 2), output_name);




%% --- NOUVEAU BLOC : Calcul et sauvegarde des métriques ---


% 1. Calculer toutes les statistiques TD (Total Displacement)
TD_mean = mean(distance_euclidienne);
TD_std = std(distance_euclidienne);
TD_median = median(distance_euclidienne);
TD_min = min(distance_euclidienne);
TD_max = max(distance_euclidienne);
TD_sup_1 = sum(distance_euclidienne > 1);
TD_sup_2 = sum(distance_euclidienne > 2);
TD_sup_3 = sum(distance_euclidienne > 3);

% 2. Calculer toutes les statistiques STS / FD
STS_mean = mean(FD);
STS_std = std(FD);
STS_median = median(FD);
STS_min = min(FD);
STS_max = max(FD);
STS_sup_05 = sum(FD > 0.5);
STS_sup_10 = sum(FD > 1.0);
STS_sup_15 = sum(FD > 1.5);

% 3. Créer et sauvegarder un fichier texte avec ces valeurs
% Nom du fichier dynamique (ex: metrics_rp_12-fmri...txt)
metrics_filename = fullfile(out_dir, ['metrics_' rp_name '.txt']);
fid = fopen(metrics_filename, 'w');

% Écrire les en-têtes (pour copier-coller facilement dans Excel)
fprintf(fid, 'TD_mean\tTD_STD\tTD_median\tTD_max\tTD_>1mm\tTD_>2mm\tTD_>3mm\tSTS_mean\tSTS_STD\tSTS_median\tSTS_max\tSTS_>0.5mm\tSTS_>1mm\n');
% Écrire les valeurs
fprintf(fid, '%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\n', ...
    TD_mean, TD_std, TD_median, TD_max, TD_sup_1, TD_sup_2, TD_sup_3, ...
    STS_mean, STS_std, STS_median, STS_max, STS_sup_05, STS_sup_10);

fclose(fid);
fprintf('Métriques sauvegardées dans : %s\n', metrics_filename);