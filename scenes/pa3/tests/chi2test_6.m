obsFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 902, 730, 536, 409, 332, 378, 348, 470, 568, 755, 932, 1163, 1375, 1452, 1577, 1577, 1576, 1474, 1322, 1177; 10013, 7268, 5188, 3992, 3334, 3307, 3497, 4126, 5434, 7513, 10392, 13195, 15190, 16430, 16884, 16792, 17037, 16190, 14875, 12775; 17936, 9611, 5814, 4075, 3238, 2849, 3359, 4118, 5990, 10535, 20037, 34365, 37606, 38944, 39495, 39595, 39505, 38799, 37263, 32656; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2640, 23745, 49033, 68032, 74457, 66451, 45700, 20405, 1262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ];
expFrequencies = [ 32587.3, 32587.3, 32587.3, 32587.3, 32587.3, 32587.3, 32587.3, 32587.3, 32587.3, 32587.3, 32587.3, 32587.3, 32587.3, 32589.6, 32679.7, 32882.8, 32650.9, 32588.5, 32587.3, 32587.3; 28193.7, 28193.7, 28193.7, 28193.7, 28193.7, 28193.7, 28193.7, 28193.7, 28193.7, 28193.7, 28193.7, 28193.7, 28193.7, 28194.2, 28926.8, 35781.8, 28518.2, 28193.8, 28193.7, 28193.7; 25148.3, 25148.3, 25148.3, 25148.3, 25148.3, 25148.3, 25148.3, 25148.3, 25148.3, 25148.3, 25148.3, 25148.3, 25148.3, 25148.3, 25388.6, 63237.7, 25189.1, 25148.3, 25148.3, 25148.3; 22512.8, 22512.8, 22512.8, 22512.8, 22512.8, 22512.8, 22512.8, 22512.8, 22512.8, 22512.8, 22512.8, 22512.8, 22512.8, 22512.8, 22512.8, 69211.4, 22512.8, 22512.8, 22512.8, 22512.8; 20053.6, 20053.6, 20053.6, 20053.6, 20053.6, 20053.6, 20053.6, 20053.6, 20053.6, 20053.6, 20053.6, 20053.6, 20053.6, 20053.6, 20212.4, 35144, 20076.1, 20053.6, 20053.6, 20053.6; 17645.5, 17645.5, 17645.5, 17645.5, 17645.5, 17645.5, 17645.5, 17645.5, 17645.5, 17645.5, 17645.5, 17645.5, 17645.5, 17645.9, 19792.7, 30266.6, 18628, 17645.6, 17645.5, 17645.5; 15186.3, 15186.3, 15186.3, 15186.3, 15186.3, 15186.3, 15186.3, 15186.3, 15186.3, 15186.3, 15186.3, 15186.3, 15186.3, 15243.8, 19675.2, 22560.3, 18354.1, 15210.2, 15186.3, 15186.3; 12550.8, 12550.8, 12550.8, 12550.8, 12550.8, 12550.8, 12550.8, 12550.8, 12550.8, 12550.8, 12550.8, 12550.9, 12558, 13064.3, 17830.6, 19113.8, 16985, 12864.8, 12554.4, 12550.8; 9505.54, 9505.42, 9505.42, 9505.42, 9505.42, 9505.42, 9505.42, 9505.42, 9505.42, 9505.43, 9505.62, 9510.9, 9619.24, 10726, 13991.5, 15315.5, 13538.6, 10435.9, 9582.67, 9508.92; 5134.35, 5121.45, 5117.13, 5115.42, 5114.72, 5114.54, 5114.78, 5115.57, 5117.5, 5122.47, 5137.72, 5194.77, 5417.43, 6064.08, 7029.26, 7455.38, 6907, 5944.02, 5367.77, 5181.53 ];
colormap(jet);
clf; subplot(2,1,1);
imagesc(obsFrequencies);
title('Observed frequencies');
axis equal;
subplot(2,1,2);
imagesc(expFrequencies);
axis equal;
title('Expected frequencies');