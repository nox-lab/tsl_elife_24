output_parent_dir = "json/";
mkdir json;
files = dir('*.mat');
% disp(files);
for file = files'
    file_name = file.name(1:end-4);
    fid = fopen(output_parent_dir + file_name + ".json", 'w');
    disp(file_name);
    mat = load(file.name);
    json = jsonencode(mat);
    fprintf(fid, json);
end
