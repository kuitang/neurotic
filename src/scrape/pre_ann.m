%% Prepare queries for annotations for slice annotaitons.
%
% [1] http://hssl.cs.jhu.edu/wiki/doku.php?id=randal:hssl:research:brain:annotation_service
% [2] http://hssl.cs.jhu.edu/wiki/doku.php?id=randal:hssl:research:brain:ramon_hdf5_interfaces
% [3] http://hssl.cs.jhu.edu/wiki/doku.php?id=randal:hssl:research:brain:ramon_hdf5_ontology
% [4] http://hssl.cs.jhu.edu/wiki/doku.php?id=randal:hssl:research:brain:data_set_description

%% Configuration
F = fullfile(DATA_PREFIX, 'kat11info.h5');
xyres = 1000;
zres  = 100;
secret;
KBASE = ['http://openconnecto.me/emca/' KAT11_TOKEN '/'];
PUBBASE = 'http://openconnecto.me/emca/kasthuri11/';

%% Download info file
urlwrite([KBASE 'projinfo/'], F);
proj_res      = h5read(F, '/PROJECT/RESOLUTION');
proj_type     = h5read(F, '/PROJECT/TYPE');
proj_res_str  = num2str(proj_res);
proj_cube_dim = h5read(F, ['/DATASET/CUBE_DIMENSION/' proj_res_str]);
proj_xy       = h5read(F, ['/DATASET/IMAGE_SIZE/' proj_res_str]);
slice_range   = h5read(F, '/DATASET/SLICERANGE');

%% Download all ANNOIDS in this extent.
% The documentation says you can do this, but this actually doesn't work.
% Use the code below (brute force specifying ANNOIDS) instead.
% delete(Q);
% h5create(Q, '/RESOLUTION', [1 1], 'Datatype', 'int64');
% h5create(Q, '/XYZOFFSET',  [3 1], 'Datatype', 'int64');
% h5create(Q, '/CUTOUTSIZE', [3 1], 'Datatype', 'int64');
% 
% h5write(Q, '/RESOLUTION', proj_res);
% h5write(Q, '/XYZOFFSET',  [1; 1; slice_range(1)]);
% h5write(Q, '/CUTOUTSIZE', [proj_xy; slice_range(2)]);
% 
% % Since MATLAB can't post a binary file, we system out to curl
% [st, res] = system(['/usr/bin/curl -o kat11_ann.h5 --data-binary @' Q ...
%                   ' ' KBASE 'objects/type/1']);
% disp(res)
% type kat11_ann.h5

%% Brute force

%% Download bounding boxes
% This list is discovered empirically. i.e. watching which requests fail.
topid = 6075;
exclude = [2146 2208 2216 2218];

annoid = 1:topid;
annoid = annoid';
annoid(exclude) = [];

Nids = length(annoid);

tmp = 'bbox_tmp.h5';

res      = zeros(Nids, 1);
xyz_off  = zeros(Nids, 3);
xyz_dim  = zeros(Nids, 3);
conf     = zeros(Nids, 1);
kvpairs  = cell(Nids, 1);

for i = 1:Nids
    a = annoid(i);
        
    urlwrite([KBASE '/' num2str(a) '/boundingbox/'], tmp);
    head = ['/' num2str(a) '/'];
    
    res(i)       = h5read(tmp, [head 'RESOLUTION']);
    xyz_off(i,:) = h5read(tmp, [head 'XYZOFFSET']);
    xyz_dim(i,:) = h5read(tmp, [head 'XYZDIMENSION']);
    conf(i)      = h5read(tmp, [head 'METADATA/CONFIDENCE']);
    kvpairs(i)   = h5read(tmp, [head 'METADATA/KVPAIRS']);
    disp(num2str(a));
end

annotation_bboxen = dataset(annoid, res, xyz_off, xyz_dim, conf, kvpairs);
save('annotation_bboxen.mat', 'annotation_bboxen');

%% Download the cutouts from the server
load('annotation_bboxen.mat');
vox_sz = sum(double(annotation_bboxen(:,'xyz_dim')), 2);
small_idxs = vox_sz < 1000;
small_bboxen = annotation_bboxen(small_idxs,:);

for i = 1:size(small_bboxen, 1)
    vcstr = vol_cutout_string(double(small_bboxen(i,'res')), double(small_bboxen(i,'xyz_off')), double(small_bboxen(i,'xyz_dim')));
    a = num2str(double(small_bboxen(i,'annoid')));
    urlwrite([PUBBASE '/' vcstr], [DATA_PREFIX '/kat11_img_' a]);
    urlwrite([KBASE '/' vcstr], [DATA_PREFIX '/kat11_ann_' a]);
end

%% Chunk our requests
% The server will crash with chunksz = 500. 100 is the highest known good.
chunksz = 10;
parts = partition_rem(annoid', chunksz);

for i = 1:length(parts)
    % We made annoids column vectors for partitioning. But the remote
    % interface can only deal with row vectors, so we tranpose back.
    annoids_part = parts{i}';
       
    first = num2str(annoids_part(1));
    last  = num2str(annoids_part(end));
    suffix = [first '_' last '.hd5'];
    
    query_file = ['queries/kat11_q_' suffix];
    
    h5create(query_file, '/ANNOIDS', size(annoid), 'Datatype', 'int32');
    h5write(query_file, '/ANNOIDS', annoid);
end

% %% Use this code to probe ANNOIDS
% 
% ed = 15000;
% for start = 6119:ed
%     annoids = start:ed;
%     
%     delete(Q);
%     h5create(Q, '/ANNOIDS', [1 length(annoids)], 'Datatype', 'int32');
%     h5write(Q, '/ANNOIDS', annoids);
%     % Since MATLAB can't post a binary file, we system out to curl
%     [st, res] = system(['/usr/bin/curl -o ' OUT ' --data-binary @' Q ...
%                       ' ' KBASE 'objects/voxels/0']);    
%     type(OUT)
% end

%% Cleanup
clear KBASE;
clear KAT11_TOKEN;
