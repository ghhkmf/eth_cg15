obsFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1e+06, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ];
expFrequencies = [ 0, 0, 0, 0, 0, 0, 4.77963e-39, 5.7021e-29, 2.29594e-16, 2.02869e-09, 1.02472e-08, 1.19726e-12, 8.0093e-24, 1.66559e-35, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 2.39065e-32, 2.9009e-15, 0.00965547, 0.12169, 2.31612e-09, 3.82822e-25, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 3.0453e-20, 12.2791, 507.385, 1.87958e-10, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 5.54355, 63683.8, 1.30038e-23, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 8.69055e-08, 290430, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 17.3304, 34490.8, 3.92842e-30, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 2.12876e-21, 424.12, 21506.9, 8.01083e-09, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 5.05741e-34, 9.58565e-11, 133.524, 1699.75, 9.26164e-05, 1.62586e-23, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 3.44466e-36, 9.35115e-21, 8.94583e-09, 0.932172, 5.98335, 0.00011277, 8.14941e-16, 1.07948e-29, 0, 0, 0, 0, 0, 0; 2.15628e-32, 2.15628e-32, 2.15628e-32, 2.15628e-32, 2.15628e-32, 2.20321e-32, 7.447e-27, 4.97434e-19, 2.72706e-10, 1.28522e-05, 4.43371e-05, 7.21524e-08, 2.32479e-15, 2.92399e-24, 4.00766e-31, 2.15628e-32, 2.15628e-32, 2.15628e-32, 2.15628e-32, 2.15628e-32 ];
colormap(jet);
clf; subplot(2,1,1);
imagesc(obsFrequencies);
title('Observed frequencies');
axis equal;
subplot(2,1,2);
imagesc(expFrequencies);
axis equal;
title('Expected frequencies');
