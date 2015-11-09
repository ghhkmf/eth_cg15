obsFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 1162, 1266, 1177, 1181, 1243, 1179, 1242, 1194, 1210, 7510, 63770, 3335, 1175, 1122, 1205, 1253, 1158, 1177, 1233, 1210; 3598, 3678, 3493, 3563, 3674, 3623, 3570, 3561, 3699, 18641, 62296, 11914, 3674, 3583, 3587, 3679, 3677, 3587, 3601, 3673; 6081, 5957, 6059, 6088, 5910, 6126, 6146, 5986, 7116, 24867, 49826, 19338, 6590, 5962, 6035, 5963, 5968, 6017, 5807, 6085; 8423, 8241, 8282, 8326, 8499, 8480, 8605, 8568, 11265, 23086, 32400, 20147, 10234, 8496, 8301, 8330, 8331, 8506, 8454, 8368; 10835, 10841, 10780, 10873, 10843, 10797, 10811, 11336, 13045, 15498, 17066, 15137, 12398, 11305, 10871, 10668, 10962, 10908, 10656, 10779 ];
expFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 1200, 1200, 1200, 1200, 1200, 1200, 1200, 1200, 1201.14, 7512.24, 63730.1, 3275.56, 1200.14, 1200, 1200, 1200, 1200, 1200, 1200, 1200; 3600, 3600, 3600, 3600, 3600, 3600, 3600, 3600.07, 3726.35, 18588.1, 62381.5, 11784.8, 3640.48, 3600.02, 3600, 3600, 3600, 3600, 3600, 3600; 6000, 6000, 6000, 6000, 6000, 6000, 6000.06, 6013.73, 7130.37, 24943.5, 49931.6, 19227.2, 6573.57, 6005.57, 6000.01, 6000, 6000, 5999.99, 6000, 6000; 8400, 8400, 8400, 8400, 8400.01, 8400.33, 8409.49, 8614.36, 11188.6, 23004.9, 32541.7, 20206.2, 10276.4, 8527.47, 8405.28, 8400.18, 8400.01, 8399.99, 8400, 8400; 10804.3, 10804.7, 10806.1, 10809.2, 10817.8, 10842.6, 10943.2, 11375.2, 12854.3, 15640.7, 17107.1, 15130.5, 12479.3, 11249.7, 10913.2, 10835.6, 10815.6, 10808.5, 10805.8, 10804.6 ];
colormap(jet);
clf; subplot(2,1,1);
imagesc(obsFrequencies);
title('Observed frequencies');
axis equal;
subplot(2,1,2);
imagesc(expFrequencies);
axis equal;
title('Expected frequencies');
