%% Get annotations and slices within a bounding box.
%
% [1] http://hssl.cs.jhu.edu/wiki/doku.php?id=randal:hssl:research:brain:annotation_service
% [2] http://hssl.cs.jhu.edu/wiki/doku.php?id=randal:hssl:research:brain:ramon_hdf5_interfaces
% [3] http://hssl.cs.jhu.edu/wiki/doku.php?id=randal:hssl:research:brain:ramon_hdf5_ontology
% [4] http://hssl.cs.jhu.edu/wiki/doku.php?id=randal:hssl:research:brain:data_set_description

%% Configuration
F = fullfile(DATA_PREFIX, 'kat11info.h5');
IMG_F = fullfile(DATA_PREFIX, 'kat11_img.h5');
ANN_F = fullfile(DATA_PREFIX, 'kat11_ann.h5');
Q = 'kat11query.h5';

RES = 1;
% Numbers from email communication
XYZOFF = [4000 7000 1200];
XYZEXT = [1000 1000 100];

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


%% Download the annotation volume
cutout_str = vol_cutout_string(RES, XYZOFF, XYZEXT);

ann_cube_file = fullfile(DATA_PREFIX, 'kat11_ann.hd5');
urlwrite([KBASE cutout_str], ann_cube_file);

img_cube_file = fullfile(DATA_PREFIX, 'kat11_img.hd5');
urlwrite(['http://openconnecto.me/emca/kasthuri11/' cutout_str], img_cube_file);

%% Brute force
% This list is discovered empirically. i.e. watching which requests fail.
topid = 6075;
exclude = [2146 2208 2216 2218];

annoids = 1:topid;
annoids(exclude) = [];

% Download all ANNOIDS in this extent.
% The documentation says you can do this, but this actually doesn't work.
% Use the code below (brute force specifying ANNOIDS) instead.
delete(Q);
h5create(Q, '/ANNOIDS',    size(annoids), 'Datatype', 'int32');
h5create(Q, '/RESOLUTION', [1 1], 'Datatype', 'int64');
h5create(Q, '/XYZOFFSET',  [1 3], 'Datatype', 'int64');
h5create(Q, '/CUTOUTSIZE', [1 3], 'Datatype', 'int64');

h5write(Q, '/ANNOIDS', annoids);
h5write(Q, '/RESOLUTION', RES);
h5write(Q, '/XYZOFFSET',  XYZOFF);
h5write(Q, '/CUTOUTSIZE', XYZEXT);

% Since MATLAB can't post a binary file, we shell out to curl
[st, res] = system(['/usr/bin/curl -o kat11_ann.h5 --data-binary @' Q ...
                  ' ' KBASE 'objects/voxels/' proj_res_str]);
disp(res)
% Make sure we actually got an HDF5 out. Otherwise, this will error
% out.
h5info(outfile, ['/' num2str(annoids(1))]);

%% Cleanup
clear KBASE;
clear KAT11_TOKEN;
