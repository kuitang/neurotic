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

%% Download info file
secret;
KBASE = ['http://openconnecto.me/emca/' KAT11_TOKEN '/'];
urlwrite([KBASE 'projinfo/'], F);


%% Parse the h5 file
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
% This list is discovered empirically. i.e. watching which requests fail.
topid = 6075
exclude = [2146 2208 2216 2218];

annoids = 1:topid;
annoids(exclude) = [];

% Chunk our requests
chunksz = 500;
parts = partition_rem(annoids', chunksz);

for i = 1:length(parts)
    % We made annoids column vectors for partitioning. But the remote
    % interface can only deal with row vectors, so we tranpose back.
    annoids_part = parts{i}';
       
    first = num2str(annoids_part(1));
    last  = num2str(annoids_part(end));
    suffix = [first '_' last '.hd5'];
    
    query_file = ['queries/kat11_q_' suffix];
    
    h5create(query_file, '/ANNOIDS', size(annoids), 'Datatype', 'int32');
    h5write(query_file, '/ANNOIDS', annoids);
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
