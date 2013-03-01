function points=get_cov_points2(varargin)
% This function returns data points to 'plot an ellipse representing a 
%   covariance matrix'. It differs from the get_cov_points function: this 
%   function allows confidence interval of covariance to be explicitely 
%   specified (currently set to 0.5). to adjust change numeric parameter 
%   after 'conf' in properties struct definition below.
%   Eg: get_cov_point2(C,mu,'conf',0.9);

    default_properties = struct(...
      'C', [], ... % The covaraince matrix (required)
      'mu', [], ... % Center of ellipse (optional)
      'conf', 0.5, ... % Percent confidence/100
      'scale', 1, ... % Scale factor, e.g. 1e-3 to plot m as km
      'style', '', ...  % Plot style
      'clip', inf); % Clipping radius

    if length(varargin) >= 1 & isnumeric(varargin{1})
      default_properties.C = varargin{1};
      varargin(1) = [];
    end

    if length(varargin) >= 1 & isnumeric(varargin{1})
      default_properties.mu = varargin{1};
      varargin(1) = [];
    end

    if length(varargin) >= 1 & isnumeric(varargin{1})
      default_properties.conf = varargin{1};
      varargin(1) = [];
    end

    if length(varargin) >= 1 & isnumeric(varargin{1})
      default_properties.scale = varargin{1};
      varargin(1) = [];
    end

    if length(varargin) >= 1 & ~ischar(varargin{1})
      error('Invalid parameter/value pair arguments.') 
    end

    prop = getopt(default_properties, varargin{:});
    C = prop.C;

    if isempty(prop.mu)
      mu = zeros(length(C),1);
    else
      mu = prop.mu;
    end

    conf = prop.conf;
    scale = prop.scale;
    style = prop.style;

    if conf <= 0 | conf >= 1
      error('conf parameter must be in range 0 to 1, exclusive')
    end

    [r,c] = size(C);
    if r ~= c | (r ~= 2 & r ~= 3)
      error(['Don''t know what to do with ',num2str(r),'x',num2str(c),' matrix'])
    end



    x0=mu(1);
    y0=mu(2);

    % Compute quantile for the desired percentile
    k = sqrt(qchisq(conf,r)); % r is the number of dimensions (degrees of freedom)



    % Make the matrix has positive eigenvalues - else it's not a valid covariance matrix!
    if any(eig(C) <=0)
    error('The covariance matrix must be positive definite (it has non-positive eigenvalues)')
    end

    [x,y,z] = getpoints(C,prop.clip);
    % h1=plot(scale*(x0+k*x),scale*(y0+k*y),prop.style);   %%%%%
    % set(h1,'zdata',z+1)									%%%%%
    % if nargout
    % h=h1;
    % end



    points = [scale*(x0+k*x), scale*(y0+k*y)];






     %---------------------------------------------------------------
    % getpoints - Generate x and y points that define an ellipse, given a 2x2
    %   covariance matrix, C. z, if requested, is all zeros with same shape as
    %   x and y.
    function [x,y,z] = getpoints(C,clipping_radius)

    n=25; % Number of points around ellipse
    p=0:pi/n:2*pi; % angles around a circle

    [eigvec,eigval] = eig(C); % Compute eigen-stuff
    xy = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
    x = xy(:,1);
    y = xy(:,2);
    z = zeros(size(x));

    % Clip data to a bounding radius
    if nargin >= 2
      r = sqrt(sum(xy.^2,2)); % Euclidian distance (distance from center)
      x(r > clipping_radius) = nan;
      y(r > clipping_radius) = nan;
      z(r > clipping_radius) = nan;
    end

    %---------------------------------------------------------------
    function x=qchisq(P,n)
    % QCHISQ(P,N) - quantile of the chi-square distribution.
    if nargin<2
      n=1;
    end

    s0 = P==0;
    s1 = P==1;
    s = P>0 & P<1;
    x = 0.5*ones(size(P));
    x(s0) = -inf;
    x(s1) = inf;
    x(~(s0|s1|s))=nan;

    for ii=1:14
      dx = -(pchisq(x(s),n)-P(s))./dchisq(x(s),n);
      x(s) = x(s)+dx;
      if all(abs(dx) < 1e-6)
        break;
      end
    end

    %---------------------------------------------------------------
    function F=pchisq(x,n)
    % PCHISQ(X,N) - Probability function of the chi-square distribution.
    if nargin<2
      n=1;
    end
    F=zeros(size(x));

    if rem(n,2) == 0
      s = x>0;
      k = 0;
      for jj = 0:n/2-1;
        k = k + (x(s)/2).^jj/factorial(jj);
      end
      F(s) = 1-exp(-x(s)/2).*k;
    else
      for ii=1:numel(x)
        if x(ii) > 0
          F(ii) = quadl(@dchisq,0,x(ii),1e-6,0,n);
        else
          F(ii) = 0;
        end
      end
    end

    %---------------------------------------------------------------
    function f=dchisq(x,n)
    % DCHISQ(X,N) - Density function of the chi-square distribution.
    if nargin<2
      n=1;
    end
    f=zeros(size(x));
    s = x>=0;
    f(s) = x(s).^(n/2-1).*exp(-x(s)/2)./(2^(n/2)*gamma(n/2));

    %---------------------------------------------------------------
    function properties = getopt(properties,varargin)
    %GETOPT - Process paired optional arguments as 'prop1',val1,'prop2',val2,...
    %
    %   getopt(properties,varargin) returns a modified properties structure,
    %   given an initial properties structure, and a list of paired arguments.
    %   Each argumnet pair should be of the form property_name,val where
    %   property_name is the name of one of the field in properties, and val is
    %   the value to be assigned to that structure field.
    %
    %   No validation of the values is performed.
    %
    % EXAMPLE:
    %   properties = struct('zoom',1.0,'aspect',1.0,'gamma',1.0,'file',[],'bg',[]);
    %   properties = getopt(properties,'aspect',0.76,'file','mydata.dat')
    % would return:
    %   properties = 
    %         zoom: 1
    %       aspect: 0.7600
    %        gamma: 1
    %         file: 'mydata.dat'
    %           bg: []
    %
    % Typical usage in a function:
    %   properties = getopt(properties,varargin{:})

    % Process the properties (optional input arguments)
    prop_names = fieldnames(properties);
    TargetField = [];
    for ii=1:length(varargin)
      arg = varargin{ii};
      if isempty(TargetField)
        if ~ischar(arg)
          error('Propery names must be character strings');
        end
        f = find(strcmp(prop_names, arg));
        if length(f) == 0
          error('%s ',['invalid property ''',arg,'''; must be one of:'],prop_names{:});
        end
        TargetField = arg;
      else
        % properties.(TargetField) = arg; % Ver 6.5 and later only
        properties = setfield(properties, TargetField, arg); % Ver 6.1 friendly
        TargetField = '';
      end
    end
    if ~isempty(TargetField)
      error('Property names and values must be specified in pairs.');
    end
