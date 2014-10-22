clear all
clf
hold on

% Round-off error results for different orderings of 3D GM rules

% max error
% s, less-than-sorted        unsorted                abs-less-than           abs-greater-than        random                  less-than-interleaved
data = [
  1, 6.9388939039072284e-18, 3.4694469519536142e-18  3.4694469519536142e-18  6.9388939039072284e-18  2.7755575615628914e-17  6.9388939039072284e-18
  2, 5.5511151231257827e-17, 3.4694469519536142e-17  4.8572257327350599e-17  8.3266726846886741e-17  2.7755575615628914e-17  1.3877787807814457e-17
  3, 1.3877787807814457e-16, 5.2041704279304213e-17  3.0531133177191805e-16  1.6653345369377348e-16  1.3877787807814457e-16  1.6653345369377348e-16
  4, 1.4155343563970746e-15, 1.0755285551056204e-16  2.2898349882893854e-16  1.1102230246251565e-16  3.6082248300317588e-16  1.0824674490095276e-15
  5, 3.6359804056473877e-15, 2.2759572004815709e-15  6.8278716014447127e-15  6.8278716014447127e-15  1.8041124150158794e-16  2.7478019859472624e-15
  6, 7.7229889150487452e-15, 1.1657341758564144e-14  9.7422070410857486e-15  1.0463852007092100e-14  6.3837823915946501e-16  3.4972025275692431e-15
  7, 4.8544501751734970e-14, 2.1593837828959295e-14  4.7378767575878555e-14  5.0626169922907138e-14  5.4123372450476381e-15  9.7422070410857486e-15
  8, 2.7641777755604835e-13, 3.0531133177191805e-14  7.0360384185619296e-14  6.8722805224297190e-14  2.4494295480792516e-15  1.6792123247455493e-14
  9, 8.1121220851798626e-13, 1.4993561947562739e-13  1.5404344466674047e-14  1.8374191057546341e-14  8.5764728652293343e-15  8.9026008787129740e-15
 10, 2.8433644327918728e-12, 3.5063618675224006e-13  4.5594084063793616e-13  4.6787573815265660e-13  2.6839641620313159e-14  3.3287261835823756e-13
    ];

s = data(:,1);
std_sort = log10(data(:,2));
unsorted = log10(data(:,3));
abs_sort = log10(data(:,4));
abs_grtr = log10(data(:,5));
rnd_sort = log10(data(:,6));
intrleav = log10(data(:,7));

ph(1) = plot(s, std_sort, 'bo-');
ph(2) = plot(s, unsorted, 'ro-');
ph(3) = plot(s, abs_sort, 'ko-');
ph(4) = plot(s, abs_grtr, 'mo-');
ph(5) = plot(s, rnd_sort, 'co-');
ph(6) = plot(s, intrleav, 'go-');
set(ph, 'linewidth', 6);
legend('less-than', 'unsorted', 'abs-less-than', 'abs-greater-than', 'random', 'less-than-interleaved', 'location', 'northwest');

xlabel('s');
ylabel('log_{10}(max error)');

% Print to PDF
set (gcf, "papersize", [11, 8.5]);
set (gcf, "paperorientation", 'landscape');

% I was using these paper settings in older versions of Octave, but
% they changed in 3.8.0
is_380 = strcmp(version(), '3.8.0');

if (!is_380)
  set (gcf, "paperposition", [0.25, 0.25, 10.75, 8.25]);
else
  % In Octave 3.8.0, the default paperposition is [0.25000, 2.50000, 8.00000, 6.00000],
  % the third number makes the plot taller instead of wider!
  set (gcf, "paperposition", [0.25, 0.25, 8.0, 10.5]);
end

print('-dpdf', 'gm_rule.pdf', '-FHelvetica:20');
