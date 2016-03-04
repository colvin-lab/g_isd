function MDSCluster = MDScluster(varargin)
% function MDSCluster = MDScluster(varargin)
%
% Arguments  :
% 'NClust'   : Max number of clusters to check. If 'index' is left
%            : unspecified then exactly 'NClust' clusters are solved for.
%            : (numeric)
%            : 
% 'index'    : Cluster index type. Options 'dunn' or 'db' optimize.
%            : Anything else skips optimization and makes maxclusters.
%            : Suggest using the optimization on the first run and then 
%            : use the optimal number of clusters in subsequent runs.
%            : 
% 'NSims'    : Display N simulations independently (numeric).
%            : 
% 'NSteps'   : Number of time steps per simulation. Different number of
%            : time steps per simulation is not supported. The function
%            : attempts to figure this out if it is not provided, but
%            : certain options in trjcat (specifically -cat) can make this
%            : less reliable. (numeric)
%            : 
% 'NAvg'     : Runs an averaging window of size 2 * NAvg + 1. (numeric)
%            : 
% 'NSkip'    : In the figure, only display every NSkip'th point. (numeric)
% 'Radius'   : Sphere size. (numeric)
%            : 
% 'Res'      : Sphere resolution. (numeric)
%            : 
% 'CMap'     : Choose a custom color map. Must have same number of colors
%            : as there are clusters. (N-by-3 numeric matrix)
%            : 
% 'bShow'    : Show the clustering output. (boolean)
%            : 
% 'Vis3D'    : Better, but causes error in Octave. (logical)
%            : 
% 'bGrid'    : Show grid lines. (logical)
%            : 

defNClust   = -1;
defIndex    = '';
defNSims    = -1;
defNSteps   = -1;
defNSkip    = -1;
defNAvg     = -1;
defRadius   = -1.0;
defRes      = -1.0;
defCMap     = NaN;
defbShow    = false;
defVis3D    = false;
defbGrid    = false;

p = inputParser;
addRequired(p, 'MDS');
addOptional(p, 'NClust', defNClust, @isnumeric);
addOptional(p, 'index',  defIndex,  @ischar);
addOptional(p, 'NSims',  defNSims,  @isnumeric);
addOptional(p, 'NSteps', defNSteps, @isnumeric);
addOptional(p, 'NSkip',  defNSkip,  @isnumeric);
addOptional(p, 'NAvg',   defNAvg,   @isnumeric);
addOptional(p, 'Radius', defRadius, @isnumeric);
addOptional(p, 'Res',    defRes,    @isnumeric);
addOptional(p, 'CMap',   defCMap,   @isnumeric);
addOptional(p, 'bShow',  defbShow);
addOptional(p, 'Vis3D',  defVis3D);
addOptional(p, 'bGrid',  defbGrid);

% Run parser and store local variables.
parse(p, varargin{:});
MDS    = p.Results.MDS;
NSims  = p.Results.NSims;
NSteps = p.Results.NSteps;
NSkip  = p.Results.NSkip;
NAvg   = p.Results.NAvg;
NClust = p.Results.NClust;


% Test NAvg, NSkip, NSims, and NClust.
if (rem(NSims,  1) ~= 0)
  error('NSims  should have a positive integer value.')
end

if (rem(NSteps, 1) ~= 0)
  error('NSteps should have a positive integer value.')
end

if (rem(NSkip,  1) ~= 0)
  error('NSkip  should have a positive integer value.')
end

if (rem(NAvg,   1) ~= 0)
  error('NAvg   should have a positive integer value.')
end

if (rem(NClust, 1) ~= 0)
  error('NClust should have a positive integer value.')
end

% Rearrange MDS matrix by simulation. MDSmat is MDS reformatted to:
%   MDSmat(time_step, dimension - out of 6, simulation_number)
nframes = size(MDS, 1);
if (NSims < 1)
  NSims  = 1;
  NSteps = nframes;
  MDSmat = MDS;
else
  NSims  = fix(NSims);
  if (NSteps < 1)
    % This should work even if there is one extra frame on the very last
    % simulation (common). But it will probably fail otherwise if the 
    % number of steps varies between simulations.
    NSteps = fix(nframes / NSims);
  end
  MDSmat = zeros(NSteps, 6, NSims);
  for i = 1:NSims
    i1 = (i - 1) * NSteps + 1;
    i2 = i * NSteps;
    MDSmat(:, :, i) = MDS(i1:i2, :);
  end
end
% Apply averaging filter. Saves result to MDSout.
MDSout = MDSmat;
if (NAvg > 2)
  if ((2 * NAvg + 1) > NSteps)
    error('Averaging window size is greater than number of steps.');
  end
  NAvg = fix(NAvg);
  for i = 1:NSims
    for j = 1:NSteps
      j1 = j - NAvg;
      if (j1 < 1)
        j1 = 1;
      end
      j2 = j + NAvg;
      if (j2 > NSteps)
        j2 = NSteps;
      end
      MDSout(j, :, i) = mean(MDSmat(j1:j2, :, i), 1);
    end
  end
  
  % Convert MDSout back to 2D.
  Out2D = true;
  if (Out2D)
    MDSmat = zeros(NSteps * NSims, 6);
    for i = 1:NSims
      i1 = (i - 1) * NSteps + 1;
      i2 = i * NSteps;
      MDSmat(i1:i2, :) = MDSout(:, :, i);
    end
    if (Out2D)
      MDSout = MDSmat;
    end
  end
end

% Setup for clustering.
Y  = pdist(MDSout);
SY = squareform(Y);
Z  = linkage(Y, 'average');

% Clustering.
index = p.Results.index;
if (strcmp(index, ''))
  if (NClust < 1)
    error('Number of clusters not given and not optimized.')
  else
    maxclust = NClust;
  end
else
  if (NClust < 1)
    maxclust = 100;
  else
    maxclust = NClust;
  end
end
if (strcmp(index, 'dunn') || strcmp(index, 'db'))
  % Choose best number of clusters.
  Cindex = zeros(maxclust, 1);
  for n = 2:maxclust
    T = cluster(Z, 'maxclust', n);
    if (strcmp(index, 'dunn'))
      Cindex(n) = dunn_index(SY, T);
    end
    if (strcmp(index, 'db'))
      Cindex(n) = DB_index(SY, T);
    end
  end
  
  % Optimize number of clusters based on index.
  if (strcmp(index, 'dunn'))
    [~, NClust] = max(Cindex);
  elseif (strcmp(index, 'db'))
    Cindex(1)   = 1;
    [~, NClust] = min(Cindex);
  end
else
  
  % Do not optimize number of clusters. Use maxclust instead.
  NClust = maxclust;
end
T = cluster(Z, 'maxclust', NClust);

% Count number in each cluster.
counts = zeros(NClust, 1);
for i = 1:NClust
  counts(i) = sum(T == i);
end

% Reorder by cluster size.
[counts, resort] = sort(counts, 1, 'descend');
OrderedT = zeros(size(T));
for i = 1:NClust
  OrderedT(T == resort(i)) = i;
end
T = OrderedT;

% Calculate center of each cluster.
centers = zeros(NClust, 6);
for i = 1:NClust
  for j = 1:6
    centers(i, j) = mean(MDSout(T == i, j));
  end
end

% Calculate the time step closest to the center of each cluster.
% The result does not use the processed (filtered) MDS.
closest = zeros(NClust, 1);
for i = 1:NClust
  iCenter = centers(i, :);
  iClustr = MDSout;
  iClustr(T ~= i) = NaN;
  iDists  = pdist2(iClustr, iCenter);
  [~, closest(i)] = min(iDists);
end

% For each cluster, calculate the number of structures from each
% simulation. This won't work unless there are the same number of 
% steps per simulation (an extra frame at the end is automatically 
% ignored).
simcount = zeros(NClust, NSims);
for i = 1:NSims
  for j = 1:NSteps
    ij = (i - 1) * NSteps + j;
    simcount(T(ij), i) = simcount(T(ij), i) + 1;
  end
end


if (p.Results.bShow)
  % Display a subset of data.
  if (NSkip < 2)
    NSkip = 1;
  else
    NSkip = fix(NSkip);
  end
  
  % Calculate plot limits.
  bsize = max(max(MDSout(:, 1:3), [], 1) - min(MDSout(:, 1:3), [], 1));
  if (p.Results.Radius < 0.0)
    R = 0.01 * bsize;
  else
    R = p.Results.Radius;
  end
  %bctr  = mean(MDSout(:, 1:3), 1);
  bmin  = min(min(MDSout(:, 1:3), [], 1) - R);
  bmax  = max(max(MDSout(:, 1:3), [], 1) + R);
  
  % Split MDS by dimensions.
  x = MDSout(:, 1);
  y = MDSout(:, 2);
  z = MDSout(:, 3);
  
  % Calculate spheres.
  if (p.Results.Res < 0.0)
    Res = 6;
  else
    Res = p.Results.Res;
  end
  [Sx, Sy, Sz] = sphere(Res);
  Sx = R * Sx; Sy = R * Sy; Sz = R * Sz;
  
  % Assign the color map.
  if (isnan(p.Results.CMap))
    cmap = jet(NClust);
  else
    % Check that the colormap meets the requirements.
    if (size(p.Results.CMap, 1) == NClust)
      cmap = p.Results.CMap;
    else
      warning(['Color map size does not match number of clusters. '...
              'Generating color map automatically.']) %#ok<WNTAG>
      cmap = jet(NClust);
    end
  end
  
  % Choose axis display style.
  figure;
  axis([bmin, bmax, bmin, bmax, bmin, bmax]);
  if (p.Results.Vis3D)
    axis('vis3d');
  else
    axis('equal');
  end
  if (p.Results.bGrid)
    grid on;
  end
  hold on;
  
  % Displays the graphics.
  n = size(MDSout, 1);
  for i = 1:n
    if (mod(i, NSkip) == 0)
      c = cmap(T(i), :);
      h = surf(Sx + x(i), Sy + y(i), Sz + z(i));
      set(h, 'FaceColor', c, 'EdgeColor', 'none');
    end
  end
  
  hold off;
end

MDSCluster.N        = NClust;
MDSCluster.T        = T;
MDSCluster.MDS      = MDS;
MDSCluster.NSims    = NSims;
MDSCluster.NSteps   = NSteps;
MDSCluster.counts   = counts;
MDSCluster.centers  = centers;
MDSCluster.closest  = closest;
MDSCluster.simcount = simcount;

end


function index = dunn_index(SY, T)
% function index = dunn_index(SY, T)
%
% SY: Output of the pdist function (squareform).
% T:  Output of the cluster function.

maxT = max(T);
if (maxT < 2)
    error('Dunn index requires at least two clusters.');
end

min_dist = zeros(maxT, 1);
max_dist = zeros(maxT, 1);
for i = 1:maxT
    Ti = (T == i);
    Tj = (T ~= i);
    if (length(Ti) == 1)
        min_dist(i) = NaN;
        max_dist(i) = NaN;
        continue;
    end
    % Between cluster vs within cluster distances.
    bc_dist = SY(Ti, Tj);
    wc_dist = SY(Ti, Ti);
    min_dist(i) = min(min(bc_dist));
    max_dist(i) = max(max(wc_dist));
end

index = min(min_dist) / max(max_dist);

end


function index = DB_index(SY, T)
% function index = DB_index(SY, T)
%
% SY: Output of pdist function (squareform).
% T:  Output of cluster function.


% Error check.
maxT = max(T);
if (maxT < 2)
    error('DB index requires at least two clusters.');
end

% Calculate the within cluster scatter.
D = zeros(maxT, 1);
S = zeros(maxT, 1);
for i = 1:maxT
    Ti = (T == i);
    if (length(Ti) == 1)
        S(i) = NaN;
        continue;
    end
    S(i) = sqrt(mean(mean(SY(Ti, Ti) .^ 2))) / 2;
end

% Mean distance between clusters.
for i = 1:maxT
    if (isnan(S(i)))
        D(i) = 1;
        continue;
    end
    R  = zeros(maxT, 1);
    Ti = (T == i);
    for j = 1:maxT
        if (i == j)
            R(j) = 0;
            continue;
        end
        Tj   = (T == j);
        M    = sqrt(mean(mean(SY(Ti, Tj) .^ 2)));
        R(j) = (S(i) + S(j)) / M;
    end
    % Ratio of internal scatter to external distance for each cluster.
    D(i) = max(R);
end

index = mean(D);

end
