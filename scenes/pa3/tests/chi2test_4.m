obsFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1e+06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ];
expFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 7.98349e-36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 2.32026e-23, 6.7964e-21, 2.76909e-29, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 2.43231e-26, 1.64453e-13, 5.79026e-12, 3.09084e-16, 3.50486e-34, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1.08447e-29, 1.75247e-14, 4.11577e-07, 5.64759e-06, 9.07972e-09, 7.11833e-19, 1.44952e-36, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 2.49562e-30, 7.31781e-17, 1.95536e-07, 0.0124238, 0.0908177, 0.000895113, 3.56626e-10, 4.21604e-21, 1.07609e-35, 0, 0, 0, 0, 0, 0; 0, 0, 0, 2.96299e-39, 4.69999e-28, 1.26975e-17, 5.91214e-09, 0.0088094, 17.0187, 77.7443, 2.6839, 0.000135353, 1.08351e-11, 4.57318e-21, 6.9622e-32, 0, 0, 0, 0, 0; 2.80846e-40, 6.26114e-36, 5.65009e-30, 4.30509e-23, 5.51781e-16, 2.4366e-09, 0.000962198, 11.4204, 1825.49, 5458.72, 525.687, 0.710625, 1.59638e-05, 1.54053e-11, 2.04144e-18, 1.58451e-25, 3.62895e-32, 1.20015e-37, 1.43809e-41, 0; 3.07832e-21, 7.26288e-19, 1.5235e-15, 1.16458e-11, 1.30413e-07, 0.000941704, 2.17057, 730.259, 15232.7, 23016, 7686.65, 128.499, 0.179183, 4.83451e-05, 5.30041e-09, 4.88456e-13, 8.99848e-17, 8.25741e-20, 9.08225e-22, 2.72936e-22; 0.0601195, 0.0601198, 0.060129, 0.060534, 0.0821825, 1.10708, 45.1776, 1508.51, 10797.8, 14360.2, 6892.63, 503.868, 12.3383, 0.343259, 0.0657243, 0.0602278, 0.0601224, 0.0601195, 0.0601195, 0.0601195 ];
colormap(jet);
clf; subplot(2,1,1);
imagesc(obsFrequencies);
title('Observed frequencies');
axis equal;
subplot(2,1,2);
imagesc(expFrequencies);
axis equal;
title('Expected frequencies');