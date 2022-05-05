function ES1D_VlasovSim_v2(Plasma, Solver)

%----------------------------------------------------------
%   1-D electrostatic Vlasov simulator entrance
%----------------------------------------------------------

%
%	Description
%
%   This function is entrance of function BSL2D_ES_vf.
%   It will read initial conditions and grid parameters
%   from xml files and construct input of Vlasov
%   simulator
%

%
%   Parameters
%
%   Plasma -> Path of .xml file which contains the initial conditions
%   Solver -> Path .xml file which contains the grid parameters
%   varargin -> Other configurable parameters
%

%
%   Acceptable input parameters
%
%   Plasma and Solver must be validate file paths
%   Now, varargin can only be 'filtration' or path of post-
%   processing code
%

%
%   Author: Guanshan Pu; Last modified: 2021.06.02
%

dom = xmlread(Plasma);
cells = dom.getElementsByTagName('Specie');

%----Find number of particle species----
nspecies = 0;

for i = 0:9999999999
    if ~isempty(cells.item(i))
        nspecies = nspecies + 1;
    else
        break
    end
end

np = zeros(1, nspecies);
qc = zeros(1, nspecies);
qm = zeros(1, nspecies);

xdistr = cell(1, nspecies);
vdistr = cell(1, nspecies);

isSave = false(1, nspecies);

nameList = cell(nspecies, 1);

%----Read the charge, charge-to-mass ratio, density and initial distribution----
for i = 0:nspecies - 1
    temp = strtrim(char(cells.item(i).getAttribute('save')));
    if strcmp(temp, 'true')
        isSave(i + 1) = true;
    end
    
    nameList{i+1} = strtrim(char(cells.item(i).getAttribute('name')));
    
    qc(i+1) = str2double(cells.item(i). ...
        getElementsByTagName('charge').item(0).getTextContent());
    qm(i+1) = str2double(cells.item(i). ...
        getElementsByTagName('cmratio').item(0).getTextContent());
    np(i+1) = str2double(cells.item(i). ...
        getElementsByTagName('density').item(0).getTextContent());
    
    nxdistr = str2double(cells.item(i).getElementsByTagName('xdistr'). ...
        item(0).getAttribute('number'));
    
    xdistrc = cells.item(i).getElementsByTagName('xdistr');
    unitfun_x = @(x)0;
    for j = 0:nxdistr-1
        fun = str2func(strtrim(char(xdistrc.item(0). ...
            getElementsByTagName('distr').item(j).getTextContent())));
        unitfun_x = @(x)unitfun_x(x) + fun(x);
    end
    
    nvdistr = str2double(cells.item(i).getElementsByTagName('vdistr'). ...
        item(0).getAttribute('number'));
    
    vdistrc = cells.item(i).getElementsByTagName('vdistr');
    unitfun_v = @(v)0;
    for j = 0:nvdistr-1
        fun = str2func(strtrim(char(vdistrc.item(0). ...
            getElementsByTagName('distr').item(j).getTextContent())));
        unitfun_v = @(v)unitfun_v(v) + fun(v);
    end
     
    xdistr{i + 1} = unitfun_x;
    vdistr{i + 1} = unitfun_v;
end

dom = xmlread(Solver);

%----Read grid parameters----
device = char(dom.getElementsByTagName('Solver').item(0).getAttribute('device'));
if strcmp(device, 'GPU')
    isGPU = true;
else
    isGPU = false;
end

cells = dom.getElementsByTagName('Grid');

Lx = str2double(cells.item(0).getElementsByTagName('xmax'). ...
    item(0).getTextContent());
ngx = str2double(cells.item(0).getElementsByTagName('xngrids'). ...
    item(0).getTextContent())-1;

Vma = zeros(nspecies, 1);Vmi = zeros(nspecies, 1);nvx = zeros(nspecies, 1);

for i = 1:nspecies
    Vma(i) = str2double(cells.item(0).getElementsByTagName('velGrid'). ...
    item(i-1).getElementsByTagName('vmax').item(0).getTextContent());
    Vmi(i) = str2double(cells.item(0).getElementsByTagName('velGrid'). ...
    item(i-1).getElementsByTagName('vmin').item(0).getTextContent());
    nvx(i) = str2double(cells.item(0).getElementsByTagName('velGrid'). ...
    item(i-1).getElementsByTagName('vngrids').item(0).getTextContent())-1;
end

%----Read temporal parameters----
cells = dom.getElementsByTagName('Temporal');

dt = str2double(cells.item(0).getElementsByTagName('tstep'). ...
    item(0).getTextContent());
nt = str2double(cells.item(0).getElementsByTagName('ntsteps'). ...
    item(0).getTextContent());

%----Read boundary conditions----
cells = dom.getElementsByTagName('Boundary');

type = cells.item(0).getElementsByTagName('type').item(0).getTextContent();
type = strtrim(char(type));

lv = str2double(cells.item(0).getElementsByTagName('left'). ...
    item(0).getTextContent());
ltype = strtrim(char(cells.item(0).getElementsByTagName('left'). ...
    item(0).getAttribute('name')));
rv = str2double(cells.item(0).getElementsByTagName('right'). ...
    item(0).getTextContent());
rtype = strtrim(char(cells.item(0).getElementsByTagName('right'). ...
    item(0).getAttribute('name')));

%----Construct boundary condition of field solver----
if strcmp(ltype, 'Dirichlet') && strcmp(rtype, 'Dirichlet')
    ebc.d1 = lv;ebc.d2 = rv;
elseif strcmp(ltype, 'Dirichlet') && strcmp(rtype, 'Neumann')
    ebc.d1 = lv;ebc.n2 = rv;
elseif strcmp(ltype, 'Neumann') && strcmp(rtype, 'Dirichlet')
    ebc.n1 = lv;ebc.d2 = rv;
else
    error('Boundary type can only be Dirichlet or Neumann')
end

%----Read extra configurations----
cells = dom.getElementsByTagName('Config');

%->Filtration
filtr = cells.item(0).getElementsByTagName('filtration');

filtr_sw = char(filtr.item(0).getAttribute('state'));
filtr_pd = str2double(filtr.item(0).getElementsByTagName('period'). ...
    item(0).getTextContent());
filtr_od = str2double(filtr.item(0).getElementsByTagName('order'). ...
    item(0).getTextContent());
filtr_bd = str2double(filtr.item(0).getElementsByTagName('width'). ...
    item(0).getTextContent());

%->External field
external = cells.item(0).getElementsByTagName('external');

external_sw = char(external.item(0).getAttribute('state'));

nfields = 0;
for i = 0:999999999999
    if ~isempty(external.item(0).getElementsByTagName('field').item(i))
        nfields = nfields + 1;
    else
        break
    end
end

ext_magnetic = cell(nfields, 1);
ext_electric = cell(nfields, 1);

hasGlobal = false;

for k = 1:nfields
    typeName = char(external.item(0). ...
        getElementsByTagName('field').item(k-1).getAttribute('type'));
    if ismember(typeName, nameList)
        ind = find(strcmp(nameList, typeName));
        if length(ind) > 1
            ind = ind(k - 1);
        end
    elseif strcmp(typeName, 'global')
        ind = nfields;
        hasGlobal = true;
    end
    
    ext_magnetic{ind} = ...
        str2func(strtrim(char(external.item(0). ...
        getElementsByTagName('field').item(k-1). ...
        getElementsByTagName('magnetic').item(0).getTextContent())));
    ext_electric{ind} = ...
        str2func(strtrim(char(external.item(0). ...
        getElementsByTagName('field').item(k-1). ...
        getElementsByTagName('electric').item(0).getTextContent())));
end

%->Diagnostics
diag = cells.item(0).getElementsByTagName('diagnostics');

diag_fnc_path = char(diag.item(0).getElementsByTagName('path'). ...
    item(0).getTextContent());
diag_sw = char(diag.item(0).getAttribute('state'));

nvars = 0;
for i = 0:999999999999
    if ~isempty(diag.item(0).getElementsByTagName('var').item(i))
        nvars = nvars + 1;
    else
        break
    end
end

vars = cell(nvars, 1);rate = zeros(nvars, 1);

for k = 0:nvars-1
    vars{k+1} = char(diag.item(0).getElementsByTagName('var'). ...
        item(k).getElementsByTagName('name').item(0).getTextContent());
    rate(k+1) = str2double(diag.item(0).getElementsByTagName('var'). ...
        item(k).getElementsByTagName('rate').item(0).getTextContent());
end

diag = struct;

if strcmp(diag_sw, 'enable')
    diag.sw = true;
elseif strcmp(diag_sw, 'disable')
    diag.sw = false;
else
    error('')
end

diag.path = diag_fnc_path;
diag.vars = vars;
diag.rate = rate;

filter = struct;

if strcmp(filtr_sw, 'enable')
    filter.sw = true;
elseif strcmp(filtr_sw, 'disable')
    filter.sw = false;
else
    error('')
end

filter.period = filtr_pd;
filter.order = filtr_od;
filter.bandwidth = filtr_bd;

external = struct;

if strcmp(external_sw, 'enable')
    external.sw = true;
else
    external.sw = false;
end

external.mag = ext_magnetic;
external.ele = ext_electric;
external.hasGlobal = hasGlobal;

%----Construct input of Vlasov simulator----
Vsim_obj = struct;

%->Particles
Vsim_obj.density = np;
Vsim_obj.charge = qc;
Vsim_obj.qmratio = qm;
Vsim_obj.xdistr = xdistr;
Vsim_obj.vdistr = vdistr;
Vsim_obj.isSave = isSave;
Vsim_obj.isGPU = isGPU;

%->Grids
Vsim_obj.Xmax = Lx;
Vsim_obj.Vmax = Vma;
Vsim_obj.Vmin = Vmi;
Vsim_obj.Xgrids = ngx;
Vsim_obj.Vgrids = nvx;

%->Time
Vsim_obj.tstep = dt;
Vsim_obj.ntsteps = nt;

%->Boundary conditions
Vsim_obj.boundary = type;
Vsim_obj.ebc = ebc;

%->Other configurations
Vsim_obj.filter = filter;
Vsim_obj.diag = diag;
Vsim_obj.external = external;

%----Invoke simulator----
BSL2D_ES_vf_v4(Vsim_obj)
