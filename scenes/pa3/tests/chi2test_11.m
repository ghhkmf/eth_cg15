obsFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 4622, 60156, 36842, 34494, 38680, 55763, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 79181, 208031, 203715, 191269, 47379, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1117, 38751, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ];
expFrequencies = [ 21914, 22241.9, 23442.9, 27590.4, 35679.7, 39832.6, 33892.5, 26330.7, 23038.7, 22131.6, 21883.1, 21802.8, 21771.7, 21758.3, 21752.5, 21751.1, 21753.3, 21760.3, 21776.3, 21813.7; 18795.8, 18797.2, 18878.2, 21195.9, 47024.5, 80586.8, 37577.9, 20049.2, 18832.6, 18796.4, 18795.8, 18795.8, 18795.8, 18795.8, 18795.8, 18795.8, 18795.8, 18795.8, 18795.8, 18795.8; 16765.5, 16765.5, 16765.5, 16774.8, 27679.3, 148554, 20136, 16766.8, 16765.5, 16765.5, 16765.5, 16765.5, 16765.5, 16765.5, 16765.5, 16765.5, 16765.5, 16765.5, 16765.5, 16765.5; 15008.5, 15008.5, 15008.5, 15008.5, 15446.2, 78795.7, 15037.7, 15008.5, 15008.5, 15008.6, 15008.6, 15008.5, 15008.6, 15008.6, 15008.5, 15008.6, 15008.6, 15008.5, 15008.6, 15008.6; 13369.1, 13369.1, 13369.1, 13461.6, 21215.5, 32637.1, 17440.9, 13393.1, 13369.1, 13369.1, 13369.1, 13369.1, 13369.1, 13369.1, 13369.1, 13369.1, 13369.1, 13369.1, 13369.1, 13369.1; 11763.7, 11765, 11854.7, 13938.9, 21852.1, 19574.4, 20517.1, 13005.5, 11803.8, 11764.2, 11763.7, 11763.7, 11763.7, 11763.7, 11763.7, 11763.7, 11763.7, 11763.7, 11763.7, 11763.7; 10136.6, 10273, 11367.7, 15372.5, 17498.1, 13340.6, 17871.5, 14301.2, 10954.9, 10214.1, 10131.4, 10124.7, 10124.3, 10124.2, 10124.2, 10124.2, 10124.2, 10124.2, 10124.3, 10125.1; 8682.67, 9528.32, 11658.4, 14216.6, 13439, 10847.1, 14028, 13810.9, 11093.9, 9262.44, 8603.2, 8427.27, 8384.54, 8373.7, 8370.71, 8370.15, 8371.08, 8374.97, 8389.34, 8446.96; 7842.76, 9118.52, 10658.4, 11419.1, 10531.9, 9643.77, 10831.8, 11387.3, 10353.3, 8806.88, 7644.18, 6999.46, 6690.63, 6551.86, 6494.29, 6480.59, 6502.07, 6572.38, 6737.34, 7100.64; 6437.6, 7072.5, 7557.99, 7701.21, 7536.42, 7422.53, 7582.23, 7703.17, 7479.23, 6943.24, 6304.4, 5737.55, 5315.05, 5039.52, 4890.38, 4849.86, 4912.39, 5086.4, 5391.49, 5846.24 ];
colormap(jet);
clf; subplot(2,1,1);
imagesc(obsFrequencies);
title('Observed frequencies');
axis equal;
subplot(2,1,2);
imagesc(expFrequencies);
axis equal;
title('Expected frequencies');
