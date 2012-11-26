function dl_annotations(outfile, kbaseurl, annoids, res)
% Download annotation hdf files in one request.
% kbaseurl takes the form 'http://openconnecto.me/emca/TOKEN/'
% annoids is row vector of annotation ids.
% res is the resolution to send to server. Should be the same as the
%   reported resolution from the info hd5 file.

    assert(size(annoids, 1) == 1, 'annoids must be row vector!');

    Q = tempname;    
    h5create(Q, '/ANNOIDS', size(annoids), 'Datatype', 'int32');
    h5write(Q, '/ANNOIDS', annoids);
    
    % Since MATLAB can't post a binary file, we system out to curl
    [st, res] = system(['/usr/bin/curl -o ' outfile ' --data-binary @' Q ...
                        ' ' kbaseurl 'objects/voxels/' num2str(res)]);

    % Make sure we actually got an HDF5 out. Otherwise, this will error
    % out.
    h5info(outfile, ['/' num2str(annoids(1))]);
end
