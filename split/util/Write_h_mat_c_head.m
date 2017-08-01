h=rand(1,5); % Filter coefficients
g = rand(10,1);
fd=fopen('myheader.h','wt');
fprintf(fd,'float h[%d]={%.9g',length(h),h(1));
fprintf(fd,',\n %.9g',h(2:end));
fprintf(fd,'};\n');


fprintf(fd,'float g[%d]={%.9g',length(g),g(1));
fprintf(fd,',\n %.9g',g(2:end));

fprintf(fd,'};\n');
fclose(fd)
