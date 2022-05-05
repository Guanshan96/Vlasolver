function BSL2D_ES_vf_v4(Vsim_obj)
global ngx nvx

%-------------------------------------------------
%   1-D periodic electrostatic Vlasov simulator
%-------------------------------------------------

%
%	Description
%
%	Using backward semi-Lagrangian (BSL) method to solve 
%	1-D periodic Vlasov-Poisson system (electrostatic Vlasov).
%

%
%   Parameters
%
%   Vsim_obj -> A structure which contains grid parameters
%               and initial conditions
%   varargin -> Other configurable parameters
%

%
%   Generally, the user needs not to set Vsim_obj manually.	
%   It can be generate through function ES1D_VlasovSim
%

%
%   Author: Guanshan Pu; Last modified: 2021.06.01
%

%----Grid parameters in phase space----
%--------------------------------------
ngx = Vsim_obj.Xgrids;
nvx = Vsim_obj.Vgrids;

Lx = Vsim_obj.Xmax; Vmi = Vsim_obj.Vmin; Vma = Vsim_obj.Vmax;

xlims = [0, Lx];
x = linspace(xlims(1), xlims(2), ngx + 1);

vxlims = [Vmi, Vma];vx = cell(length(nvx), 1);

for k = 1:length(nvx)
    vx{k} = linspace(vxlims(k, 1), vxlims(k, 2), nvx(k) + 1);
end

dx = (xlims(2) - xlims(1))/ngx;
dvx = (vxlims(:, 2) - vxlims(:, 1))./nvx;

if Vsim_obj.isGPU
    x = gpuArray(x);
    for k = 1:length(nvx)
        vx{k} = gpuArray(vx{k});
    end
end

%----Intialize distribution function----
%---------------------------------------

%->Parameters for particle system
nspecies = length(Vsim_obj.density);

np = Vsim_obj.density;       %----Number density
qc = Vsim_obj.charge;        %----Electric charge
qm = Vsim_obj.qmratio;       %----Charge to mass ratio

nx = Vsim_obj.xdistr;        %----Distribution in config space
nv = Vsim_obj.vdistr;        %----Distribution in velocity space

%->Generate discretized distribution functions
fs = cell(nspecies, 1);
tic
for k = 1:nspecies        
    if Vsim_obj.isGPU
        fs{k} = gpuArray(zeros(ngx + 1, nvx(k) + 1));
    else
        fs{k} = zeros(ngx + 1, nvx(k) + 1);
    end
    for i = 1:ngx + 1
        fs{k}(i, :) = np(k)*reshape(nx{k}(x(i))*nv{k}(vx{k}), [1 nvx(k) + 1]);
    end
end
t = toc;
disp(['Time for initialization ',num2str(t)]);


%----Temporal parameters----
%---------------------------

dt = Vsim_obj.tstep;
ntstep = Vsim_obj.ntsteps;
t_cal = 0;t = 0;

%----Diagnostics----
%-------------------
diag_names = Vsim_obj.diag.vars;
diag_rates = Vsim_obj.diag.rate;
save_inds = find(Vsim_obj.isSave);
diag_vars = Init_Diag_Vars(diag_names, diag_rates, length(save_inds), ...
    ntstep, nspecies);

if Vsim_obj.isGPU
    for i = 1:length(diag_vars)
        if ~strcmp(diag_names{i}, 'Distribution')
            diag_vars{i} = gpuArray(diag_vars{i});
        else
            for k = 1:length(save_inds)
                diag_vars{i}{k} = gpuArray(diag_vars{i}{k});
            end
        end
    end
end

%->Charge density
if Vsim_obj.isGPU
    rho = gpuArray(zeros(ngx + 1, nspecies));
else
    rho = zeros(ngx + 1, nspecies);
end

%->Select boundary condition
ibc = Vsim_obj.boundary;

switch ibc
    case 'natural'
        PCSIP_x = PCSIP_natural_matrix(ngx + 1, dx);
        
        PCSIP_v = cell(nspecies, 1);
        for k = 1:nspecies
            PCSIP_v{k} = PCSIP_natural_matrix(nvx(k) + 1, dvx(k));
        end
        
        if Vsim_obj.isGPU
            PCSIP_x = gpuArray(PCSIP_x);
            for k = 1:nspecies
                PCSIP_v{k} = gpuArray(PCSIP_v{k});
            end
        end
        
        Spinterp = @(f, dt, dx, A, F, dim)PCSIP_2d_natural(f, dt, dx, A, F, dim);
        Field = @(rho_t, Lx, bc)Natural_Field(rho_t, Lx, bc);
    case 'PCSIP_periodic'
        PCSIP_x = PCSIP_periodic_matrix(ngx + 1, dx);
        
        PCSIP_v = cell(nspecies, 1);
        for k = 1:nspecies
            PCSIP_v{k} = PCSIP_periodic_matrix(nvx(k) + 1, dvx(k));
        end
        
        if Vsim_obj.isGPU
            PCSIP_x = gpuArray(PCSIP_x);
            for k = 1:nspecies
                PCSIP_v{k} = gpuArray(PCSIP_v{k});
            end
        end
        
        Spinterp = @(f, dt, dx, A, F, dim)PCSIP_2d_periodic(f, dt, dx, A, F, dim);
        Field = @(rho_t, Lx, bc)Periodic_Field(rho_t, Lx);
    case 'PCHIP_periodic'
        zdx_fcn = @(zdy, Fx, zdx_p, zdx, dt)...
                  zdy - dt*(bsxfun(@times, Fx([2:end 2]), zdx_p([2:end 2], :)) - ...
                  bsxfun(@times, Fx([end-1 1:end-1]), zdx_p([end-1 1:end-1], :)) + ...
                  bsxfun(@times, Fx([2:end 2]), zdx([2:end 2], :)) - ... 
                  bsxfun(@times, Fx([end-1 1:end-1]), zdx([end-1 1:end-1], :)))/(4*dx);
        zdv_fcn = @(zdx_p, Fy, zdy_p, zdy, dt, ind)...
                  zdx_p - dt*(bsxfun(@times, Fy([2:end 2]), zdy(:, [2:end 2])) - ...
                  bsxfun(@times, Fy([end-1 1:end-1]), zdy(:, [end-1 1:end-1])) + ...
                  bsxfun(@times, Fy([2:end 2]), zdy_p(:, [2:end 2])) - ... 
                  bsxfun(@times, Fy([end-1 1:end-1]), zdy_p(:, [end-1 1:end-1])))/(4*dvx(ind));
        Field = @(rho_t, Lx, bc)Periodic_Field(rho_t, Lx);
        
        zdx = cell(nspecies, 1);
        zdv = cell(nspecies, 1);
        
        for j = 1:nspecies
            zdv{j} = (fs{j}(:, 1:end-4) - ...
                8*fs{j}(:, 2:end-3) + ...
                8*fs{j}(:, 4:end-1) - fs{j}(:, 5:end))/(12*dvx(j));
            
            zdv{j} = [(-3*fs{j}(:, 1)+4*fs{j}(:, 2)-fs{j}(:, 3))/(2*dvx(j)), ...
                      (-3*fs{j}(:, 2)+4*fs{j}(:, 3)-fs{j}(:, 4))/(2*dvx(j)), ...
                      zdv{j},...
                      (3*fs{j}(:, end-1)-4*fs{j}(:, end-2)+fs{j}(:, end-3))/(2*dvx(j)),...
                      (3*fs{j}(:, end-0)-4*fs{j}(:, end-1)+fs{j}(:, end-2))/(2*dvx(j))];
            
            zdx{j} = (fs{j}([end-2 end-1 1:end-2], :) - ...
                8*fs{j}([end-1 end 2:end-1], :) + ...
                8*fs{j}([2:end-1 1 2], :) - fs{j}([3:end 2 3], :))/(12*dx);
        end
        
        zdx_p = cell(nspecies, 1);
        zdv_p = cell(nspecies, 1);
end

ebc = Vsim_obj.ebc;

%----Filamentation filtration----
period = Vsim_obj.filter.period;
order = Vsim_obj.filter.order;
bandwidth = Vsim_obj.filter.bandwidth;

%----External field----
E_ext = Vsim_obj.external.ele;
if Vsim_obj.external.sw
    ext_en = 1;
else
    ext_en = 0;
end

%----First push----
for j = 1:nspecies
    if strcmp(ibc, 'natural') || strcmp(ibc, 'PCSIP_periodic')
        fs{j} = Spinterp(fs{j}, dt/2, dx, PCSIP_x, vx{j}, 1);
    else
        [fs{j}, zdx_p{j}] = PCHIP_2d_periodic(fs{j}, zdx{j}, dt/2, dx, vx{j}, 1, ...
            Vsim_obj.isGPU);
        zdv_p{j} = zdv_fcn(zdv{j}, vx{j}, zdx_p{j}, zdx{j}, dt/2, j);
    end
    rho(:,j) = qc(j)*Moments1D(fs{j}, vx{j}, 'zeroth');
end

indset = ones(length(diag_names), 1);

%----main loop----
for i = 1:ntstep

    tic
    
    %----Calculate electric field at dt/2----
    rho_t = sum(rho, 2);
    E = Field(rho_t, Lx, ebc);
    
    if Vsim_obj.external.hasGlobal
        E = E + ext_en*E_ext{end}(t, x');
    end
    
    %----Phase space advection----
    for j = 1:nspecies
        Fx = qm(j)*(E + ext_en*E_ext{j}(t, x'));
        
        %----Advection in phase space----
        if strcmp(ibc, 'PCHIP_periodic')
            [fs{j}, zdv{j}] = PCHIP_2d_periodic(fs{j}, zdv_p{j}, dt, dvx(j), Fx, 2, ...
                Vsim_obj.isGPU);
            zdx{j} = zdx_fcn(zdx_p{j}, Fx, zdv_p{j}, zdv{j}, dt);
            [fs{j}, zdx_p{j}] = PCHIP_2d_periodic(fs{j}, zdx{j}, dt, dx, vx{j}, 1, ...
                Vsim_obj.isGPU);
            zdv_p{j} = zdv_fcn(zdv{j}, vx{j}, zdx_p{j}, zdx{j}, dt, j);
        else
            fs{j} = Spinterp(fs{j}, dt, dvx(j), PCSIP_v{j}, Fx, 2);
            fs{j} = Spinterp(fs{j}, dt, dx, PCSIP_x, vx{j}, 1);
        end
        
        %----Filamentation filtration----
        if mod(i, period) == 0 && Vsim_obj.filter.sw
            fs{j} = Filtration1D(fs{j}, vxlims(j,:), order, bandwidth);
        end
        
        %----Calculate charge density on grid----
        rho(:,j) = qc(j)*Moments1D(fs{j}, vx{j}, 'zeroth');        
    end
    
    if Vsim_obj.external.hasGlobal
        E = E - ext_en*E_ext{end}(t, x)';
    end
    
    %----Calculate diagnostic parameters----
    for k = 1:length(Vsim_obj.diag.vars)
        
        varname = Vsim_obj.diag.vars{k};
        
        switch varname
            case 'Density'
                momentType = 'zeroth';
            case 'Velocity'
                momentType = 'first';
            case 'Total energy'
                momentType = 'second';
            case 'Temperature'
                momentType = 'secondm';
        end
        
        if strcmp(varname, 'Density')||strcmp(varname, 'Velocity')||...
                strcmp(varname, 'Total energy')||strcmp(varname, 'Temperature')
            if mod(i, diag_rates(k)) == 0
                for j = 1:nspecies
                    diag_vars{k}(indset(k), :, j) = Moments1D(fs{j}, vx{j}, momentType);
                end
                indset(k) = indset(k) + 1;
            end
            continue;
        end
        
        switch varname
            case 'Total charge'
                if mod(i, diag_rates(k)) == 0
                    diag_vars{k}(indset(k)) = trapz(x, rho(:, 1) + rho(:, 2));
                    indset(k) = indset(k) + 1;
                end
            case 'Distribution'
                if mod(i, diag_rates(k)) == 0
                    for j = 1:length(save_inds)
                        diag_vars{k}{save_inds(j)}(:, :, indset(k)) = fs{j};
                    end
                    indset(k) = indset(k) + 1;
                end
            case 'Electric field'
                if mod(i, diag_rates(k)) == 0
                    diag_vars{k}(indset(k), :) = E;
                    indset(k) = indset(k) + 1;
                end
        end
        
    end
    
    t_cal = t_cal + toc;
    
    t = t+dt;
    
    %----Display time elapsed----
    if mod(i, 100) == 0
        disp(['Time elapsed: ', num2str(t_cal), '; ', 'Percentage: ', ...
            num2str(100*i/ntstep), '%'])
        disp(['Total charge:', num2str(trapz(x, sum(rho, 2)))])
    end
end

%----Take results from GPU memory to host memory----
x = gather(x); %#ok<NASGU>
for k = 1:nspecies
    vx{k} = gather(vx{k});
end

for i = 1:length(diag_vars)
    if ~strcmp(diag_names{i}, 'Distribution')
        diag_vars{i} = gather(diag_vars{i});
    else
        for k = 1:length(save_inds)
            diag_vars{i}{k} = gather(diag_vars{i}{k});
        end
    end
end

%----Save diagnostic parameters into .mat files----
[filedir, ~, ~] = fileparts(Vsim_obj.diag.path);

diagdir = fullfile(filedir, 'Diag_ES.mat');
griddir = fullfile(filedir, 'Grid_ES.mat');

save(diagdir, 'diag_vars', 'diag_names', 'diag_rates', '-v7.3')
save(griddir, 'ngx', 'nvx', 'x', 'vx', 'dt', 'ntstep')

% if Vsim_obj.diag.sw
%     run(Vsim_obj.diag.path)
% end

function diag_vars = Init_Diag_Vars(diag_names, diag_rates, saves, ntstep, nspecies)
global ngx nvx
diag_vars = cell(length(diag_names), 1);

for i = 1:length(diag_names)
    points = ceil(ntstep/diag_rates(i));
    switch diag_names{i}
        case 'Density'
            diag_vars{i} = zeros(points, ngx + 1, nspecies);
        case 'Velocity'
            diag_vars{i} = zeros(points, ngx + 1, nspecies);
        case 'Distribution'
            diag_vars{i} = cell(saves, 1);
            for k = 1:saves
                diag_vars{i}{k} = zeros(ngx + 1, nvx(k) + 1, points);
            end
        case 'Total charge'
            diag_vars{i} = zeros(points, 1);
        case 'Total energy'
            diag_vars{i} = zeros(points, ngx + 1, nspecies);
        case 'Temperature'
            diag_vars{i} = zeros(points, ngx + 1, nspecies);
        case 'Electric field'
            diag_vars{i} = zeros(points, ngx + 1);
    end
end

function E = Periodic_Field(rho_t, Lx)
global ngx

    Phi = Poisson1D(-rho_t(1:ngx)', Lx);
    E = -Gradient1D(Phi, Lx);
    E = [E, E(1)]'; 
       
function E = Natural_Field(rho_t, Lx, bc)
global dx

    Phi = Poisson1D_FE(rho_t, Lx, bc);
    E = -(Phi(3:end) - Phi(1:end-2))/(2*dx);
    E = cat(1, -(Phi(2) - Phi(1))/dx, E, -(Phi(end) - Phi(end-1))/dx);
