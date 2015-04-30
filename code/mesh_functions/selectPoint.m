function idx=selectPoint(figure, mesh)

datacursormode on
dcm_obj = datacursormode(figure); 

fprintf(1,'Press any key when you have selected your vertex.\n');
pause;
a = getCursorInfo(dcm_obj);
p = a.Position;

searchResult = (mesh.vertices==repmat(p,size(mesh.vertices,1),1));
searchResult = searchResult(:,1)&searchResult(:,2)&searchResult(:,3);

idx = find(searchResult);