function pfm_mtx = initPfmMtx(general, varout, metricname)
if nargin < 3
    metricname = '';
end

try
    ismonitor = varout.ismonitor;
catch
    ismonitor = false;
    warning('The variable ''ismonitor'' in ''varout'' is missing.')
end

lj = varout.LJ;
nmethods = varout.nmethods;
nzones = varout.nzones;
numV = varout.numV;
nummu = varout.nummu;

if ismonitor
    npts = varout.nmonpts;
    points = 'monitor';
else
    npts = varout.nctrpts;
    points = 'control';
end

if strcmpi(metricname, 'ac')
    npts = 1;
end

% initialization of the performance metrics
mtx_structure = orderfields(struct('Vmax', lj, 'mu', 0, 'V', 1, ...
    'targetzone', 1, 'metric', metricname, 'scores', zeros(npts,1), ...
    'points', points));

pfm_mtx = cellfun(@(x) mtx_structure, cell(nmethods, nzones), 'UniformOutput', false);

% update the zone index
for zz = 1:nzones
    for nm = 1:nmethods
        pfm_mtx{nm,zz}.targetzone = zz;
    end
end

for nn = 1:nmethods
    if nn < 4
        for zz = 1:nzones
            pfm_mtx{nn,zz} = repmat(pfm_mtx{nn,zz}, numV, nummu);
            for vv = 1:numV
                for mm = 1:nummu
                    pfm_mtx{nn,zz}(vv,mm).V = general.userpar.V(vv);
                    pfm_mtx{nn,zz}(vv,mm).mu = general.userpar.mu(mm);
                end
            end
        end
    else
    end
end


end