

minin = [1.0347; 1.0826; 0.7275; 0.7351; 2.0817; 0.7417; 0.7556; 0.8980; 1.4480; 0.7196; 0.8013; 1.0838; 0.5539; 1.6154; 1.4343];

maxin = [2.0094; 4.1883; 2.9991; 2.9563; 4.9855; 2.9982; 4.2269; 3.8215; 3.9985; 5.1439; 5.2354; 4.7269; 3.7286; 6.1919; 5.8338];


in = [1.2734; 1.9883; 1.0921; 1.2311; 2.8458; 0.9865; 2.0600; 2.1392; 1.8975; 2.3757; 3.0535; 2.4230; 2.0261; 3.8528; 3.0020]; % minin < in < maxin

in2 = [2.2734; 1.9883; 1.0921; 1.2311; 2.8458; 0.9865; 2.0600; 2.1392; 1.8975; 2.3757; 3.0535; 2.4230; 2.0261; 3.8528; 3.0020];% minin < maxin < in2

in3 = [0.2734; 1.9883; 1.0921; 1.2311; 2.8458; 0.9865; 2.0600; 2.1392; 1.8975; 2.3757; 3.0535; 2.4230; 2.0261; 3.8528; 3.0020];% in3 < minin < maxin

in4=minin+[maxin-minin]./2; % in4 is equidistant from minin and maxin, i.e. minin < in4 < maxin

exp(-1.*[(in-minin).* (maxin-in)  (in2-minin).* (maxin-in2)  (in3-minin).* (maxin-in3)  (in4-minin).* (maxin-in4)])
% [(minin-in).* (in-maxin)  (minin-in2).* (in2-maxin)  (minin-in3).* (in3-maxin)  (minin-in4).* (in4-maxin)]

[(in-minin).* (maxin-in)]
[(in2-minin).* (maxin-in2)]
[(in3-minin).* (maxin-in3)]

sum(exp(-[(in-minin).* (maxin-in)]))
sum(exp(-[(in2-minin).* (maxin-in2)]))
sum(exp(-[(in3-minin).* (maxin-in3)]))

sum(exp(-((in-minin).*(maxin-in))))/size(in,1)
sum(exp(-((in2-minin).*(maxin-in2))))/size(in,1)
sum(exp(-((in3-minin).*(maxin-in3))))/size(in,1)