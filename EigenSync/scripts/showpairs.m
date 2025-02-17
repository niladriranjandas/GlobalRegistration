function [] = showpairs(pairs,i,j)

hold on;



scatter3(pairs{i,j}(:,1),pairs{i,j}(:,2),pairs{i,j}(:,3),10,'r')
scatter3(pairs{i,j}(:,4),pairs{i,j}(:,5),pairs{i,j}(:,6),10,'g')

hold off;
end

