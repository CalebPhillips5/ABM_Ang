clear all;
v = VideoWriter('newfile.avi','Uncompressed AVI');
v.FrameRate = 5;
open(v)
saida0_00000;
grad_vegf0_00000;
circle_cell(cells(1,1),cells(1,2),cells(1,3),cells(1,4));
xlim([0 2*cells(1,4)]);
ylim([0 2*cells(1,4)]);
hold on;
quiver(grad(:,1),grad(:,2),grad(:,3),grad(:,4));
s=size(cells);
for c = 2:s(1);
circle_cell(cells(c,1),cells(c,2),cells(c,3),cells(c,4));
end;
print -dpng saida0_00000.png
A = imread('saida0_00000.png');
writeVideo(v,A)
clf;
saida0_00500;
grad_vegf0_00500;
circle_cell(cells(1,1),cells(1,2),cells(1,3),cells(1,4));
xlim([0 2*cells(1,4)]);
ylim([0 2*cells(1,4)]);
hold on;
quiver(grad(:,1),grad(:,2),grad(:,3),grad(:,4));
s=size(cells);
for c = 2:s(1);
circle_cell(cells(c,1),cells(c,2),cells(c,3),cells(c,4));
end;
print -dpng saida0_00500.png
A = imread('saida0_00500.png');
writeVideo(v,A)
clf;
saida0_01000;
grad_vegf0_01000;
circle_cell(cells(1,1),cells(1,2),cells(1,3),cells(1,4));
xlim([0 2*cells(1,4)]);
ylim([0 2*cells(1,4)]);
hold on;
quiver(grad(:,1),grad(:,2),grad(:,3),grad(:,4));
s=size(cells);
for c = 2:s(1);
circle_cell(cells(c,1),cells(c,2),cells(c,3),cells(c,4));
end;
print -dpng saida0_01000.png
A = imread('saida0_01000.png');
writeVideo(v,A)
clf;
saida0_01500;
grad_vegf0_01500;
circle_cell(cells(1,1),cells(1,2),cells(1,3),cells(1,4));
xlim([0 2*cells(1,4)]);
ylim([0 2*cells(1,4)]);
hold on;
quiver(grad(:,1),grad(:,2),grad(:,3),grad(:,4));
s=size(cells);
for c = 2:s(1);
circle_cell(cells(c,1),cells(c,2),cells(c,3),cells(c,4));
end;
print -dpng saida0_01500.png
A = imread('saida0_01500.png');
writeVideo(v,A)
clf;
saida0_02000;
grad_vegf0_02000;
circle_cell(cells(1,1),cells(1,2),cells(1,3),cells(1,4));
xlim([0 2*cells(1,4)]);
ylim([0 2*cells(1,4)]);
hold on;
quiver(grad(:,1),grad(:,2),grad(:,3),grad(:,4));
s=size(cells);
for c = 2:s(1);
circle_cell(cells(c,1),cells(c,2),cells(c,3),cells(c,4));
end;
print -dpng saida0_02000.png
A = imread('saida0_02000.png');
writeVideo(v,A)
clf;
saida0_02500;
grad_vegf0_02500;
circle_cell(cells(1,1),cells(1,2),cells(1,3),cells(1,4));
xlim([0 2*cells(1,4)]);
ylim([0 2*cells(1,4)]);
hold on;
quiver(grad(:,1),grad(:,2),grad(:,3),grad(:,4));
s=size(cells);
for c = 2:s(1);
circle_cell(cells(c,1),cells(c,2),cells(c,3),cells(c,4));
end;
print -dpng saida0_02500.png
A = imread('saida0_02500.png');
writeVideo(v,A)
clf;
saida0_03000;
grad_vegf0_03000;
circle_cell(cells(1,1),cells(1,2),cells(1,3),cells(1,4));
xlim([0 2*cells(1,4)]);
ylim([0 2*cells(1,4)]);
hold on;
quiver(grad(:,1),grad(:,2),grad(:,3),grad(:,4));
s=size(cells);
for c = 2:s(1);
circle_cell(cells(c,1),cells(c,2),cells(c,3),cells(c,4));
end;
print -dpng saida0_03000.png
A = imread('saida0_03000.png');
writeVideo(v,A)
clf;
saida0_03500;
grad_vegf0_03500;
circle_cell(cells(1,1),cells(1,2),cells(1,3),cells(1,4));
xlim([0 2*cells(1,4)]);
ylim([0 2*cells(1,4)]);
hold on;
quiver(grad(:,1),grad(:,2),grad(:,3),grad(:,4));
s=size(cells);
for c = 2:s(1);
circle_cell(cells(c,1),cells(c,2),cells(c,3),cells(c,4));
end;
print -dpng saida0_03500.png
A = imread('saida0_03500.png');
writeVideo(v,A)
clf;
saida0_04000;
grad_vegf0_04000;
circle_cell(cells(1,1),cells(1,2),cells(1,3),cells(1,4));
xlim([0 2*cells(1,4)]);
ylim([0 2*cells(1,4)]);
hold on;
quiver(grad(:,1),grad(:,2),grad(:,3),grad(:,4));
s=size(cells);
for c = 2:s(1);
circle_cell(cells(c,1),cells(c,2),cells(c,3),cells(c,4));
end;
print -dpng saida0_04000.png
A = imread('saida0_04000.png');
writeVideo(v,A)
clf;
close(v)
exit
