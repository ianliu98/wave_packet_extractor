clear;
file = csvread('extract_sample.csv');
forf = csvread('a3ug_forf.csv');
T1 = 8000;
T2 = 15000;
size_file = size(file);
time_posi = size_file(1) * T1 / (T2 - T1);
field_forf = forf(time_posi+1:time_posi+size_file(1),:);
% plot
max_value = max(max(file));
colormap(jet)
mesh(file)
view(2)
caxis([max_value-1.5, max_value])
colorbar()
hold on

% select two points
[a, b] = ginput(2);
a = floor(a);
b = floor(b);

% group velocity
dt = 0.004;
tskip = 4;
vg_distance = 280 * (a(2) - a(1)) / size_file(2);
vg_time = (b(2) - b(1)) * (T2 - T1) / size_file(1);
vg = vg_distance / vg_time * (1 / (dt * tskip));

% wavepacket range
k = (b(2) - b(1)) / (a(2) - a(1));
m = b(1) - k*a(1);
ylim([0,size_file(1)])
x = 1:1:size_file(2);
plot(x, k*x+m)
[offset_x, offset_y] = ginput(2);
m_top = m + abs(offset_y(1) - (k*offset_x(1)+m));
m_bottom = m - abs(offset_y(2) - (k*offset_x(2)+m));
m_range = m_top - m_bottom;
m_top = floor(m_top);
m_bottom = floor(m_bottom);
plot(x, k*x+m_top)
plot(x, k*x+m_bottom)

% extract wavepacket
file2 = file;
%file2 = field_forf;
for j = 1:size_file(2)
    for i = 1:size_file(1)
        max_v = k * j + m_top;
        min_v = k * j + m_bottom;
        if (i>max_v) || (i<min_v)
            file2(i,j) = -20;
        end
    end
end
figure,
colormap(jet)
mesh(file2)
view(2)
caxis([max_value-1.5, max_value])
colorbar()

% field
file3 = 10.^(file2);
field = [];
for jj = 1:size_file(2)
    field(jj) = sum(file3(:,jj));
end
field = field / m_range;

% nonlinear growth rate
gn = [];
gn(1) = 0;
for tmp = 2:size_file(2)
    ratio = field(tmp)/field(tmp-1);
    gn(tmp) = log(ratio) * vg;
end

x_bottom_right = floor(-m_bottom / k);
ini = max(x_bottom_right,size_file(2)/2);
x_top_left = floor((size_file(1)-m_top) / k);
eend = min(x_top_left,size_file(2));
x_tick = 0:140/(eend-ini):140;

% linear
gr_linear_all = csvread('a3ug_linear_gr.csv');
gr_linear_extract = gr_linear_all(:,120:200);
gr_linear_extract_2_4 = gr_linear_extract(2:4,:);
gr_linear_extract_av = [];
for kkk = 1: 81
    gr_linear_extract_av(kkk) = sum(gr_linear_extract_2_4(:,kkk));
end
hori_tick = 0.01:0.01:1;
TT1 = T1/1e4;
TT2 = T2/1e4;
vert_tick = TT1+(TT2-TT1)/81:(TT2-TT1)/81:TT2;

% plot
figure,
subplot(2,2,1)
plot(x_tick,field(ini:eend))
ylabel('B_w/B_0')
title('Wave Amplitude')
subplot(2,2,3)
semilogy(x_tick,abs(gn(ini:eend)))
xlabel('h [c\Omega^{-1}_{e0}]')
ylabel('\Gamma_n/\Omega_{e0}')
title('Growth Rate')
subplot(2,2,[2,4])
colormap(jet)
surf(hori_tick,vert_tick,transpose(gr_linear_extract))
view(2)
colorbar()
xlabel('\omega [\Omega_{e0}]')
ylabel('t\times10^4 [\Omega^{-1}_{e0}]')
title('Linear Growth Rate')
shading interp



