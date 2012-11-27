function [ s ] = vol_cutout_string( res, xyzoff, xyzext )

    xyzend = xyzoff + xyzext;
    s = ['hdf5/' num2str(res) ];
    for i = 1:3
        s = [s '/' num2str(xyzoff(i)) ',' num2str(xyzend(i))];
    end
    s = [s '/'];

end

