% tic    
% b=0.6;
% D=Diagramm(b);
% D=AllClassify(D);
% Postroit(D);
% save(sprintf('DiagrammFor_%1.1f.mat', b),'D')
% toc


tic
b=0.80;
D=Diagramm(b);
D=AllClassify(D);
Postroit(D);
save(sprintf('DiagrammFor_%1.1f.mat', b),'D')
toc


% tic
% b=1.0;
% D=Diagramm(b);
% D=AllClassify(D);
% Postroit(D);
% save(sprintf('DiagrammFor_%1.1f.mat', b),'D')
% toc


