
% Generate the parameter combinations
v_klk = [5, 7]; 
v_eT = [1e-6 1e-7 1e-8 1e-9 1e-10];
v_s0 = [1e-1 1e-2 1e-3 1e-4 1e-5];
v_cell_v = 0.05/24;
v_n_iT = [0 0.25 0.5 0.75 1.0];

T = 0 : 1 : (25*24);
nT = length(T);

nklk = length(v_klk);
neT = length(v_eT);
ns0 = length(v_s0);
ncell_v = length(v_cell_v);
nn_iT = length(v_n_iT);
nTot = nklk*neT*ns0*ncell_v*nn_iT;

v_klk = repmat(v_klk',nTot/nklk,1);
v_eT = repmat(repelem(v_eT',nklk,1),nTot/(neT*nklk),1);
v_s0 = repmat(repelem(v_s0',nklk*neT,1),nTot/(ns0*neT*nklk),1);
v_cell_v = repmat(repelem(v_cell_v',nklk*neT*ns0,1),nTot/(ncell_v*ns0*neT*nklk),1);
v_n_iT = repelem(v_n_iT',nTot/nn_iT,1);
params = table(v_klk,v_eT,v_s0,v_cell_v,v_n_iT);

% Pre-create the output table
outputarray = cell(nTot,1);

parfor i = 1:nTot
    params_i = params(i,:);
    klk = params_i.v_klk;
    eT = params_i.v_eT;
    n_iT = params_i.v_n_iT;
    s0 = params_i.v_s0;
    v = params_i.v_cell_v;
    iT = n_iT*eT;
    %iT = n_iT * 0.44e-9;

    % Initial conditions
     i0 = max(0,eT*(n_iT-1)/s0);
     ci0 = n_iT-i0*s0/eT;
%     i0 = max(0,(iT-eT)/s0);
%     ci0 = (iT-i0*s0)/eT;
    e0 = 1-ci0;
    ic =[1,e0,0,i0,ci0]; % Note, order of elements: s, e, cs, i, ci
    
    % Solve
    [t,y] = ode15s(@(t,y) migrating_ode_system(t,y,v,eT,iT,s0,klk),T,ic);

    % Generate output data
    z = v*t;
    pH = 6.8482 - 0.3765*z - 5.1663*z.^2 + 3.1792*z.^3;
    output_data = array2table([t, z, pH, y], "VariableNames",{'t_hr', 'z', 'pH', 's','e','cs','i','ci'});
    id_vars = repmat({num2str(klk),num2str(eT),num2str(n_iT),num2str(s0), num2str(v)}, nT, 1);
    id_vars = array2table(id_vars, 'VariableNames', {'klk','eT','n_iT','s0','v'});
    
    outputarray{i} = [id_vars output_data];
end

% Concatenate the results
output = array2table(strings(nT*nTot,width(outputarray{1})));
output.Properties.VariableNames = outputarray{1}.Properties.VariableNames;
for i = 1:length(outputarray)
    output_i = ((i-1)*nT+1):i*nT;
    output(output_i,:) = outputarray{i};
end

% Write output
writetable(output,"data/movingcell_solutions.csv");