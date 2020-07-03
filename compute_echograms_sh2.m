%function [abs_echograms, rec_echograms, echograms] = compute_echograms_sh(room, src, rec, abs_wall, limits, rec_orders)
function [abs_echograms] = compute_echograms_sh2(room, src, rec, abs_wall, limits, rec_orders)
% This version returns only the absorption echogram. It modifies the nested 
% src, rec loop to save ram, in particular it doesn't save or return the 
% echograms structure, as only the absorption echogram is needed.

nRec = size(rec,1);
nSrc = size(src,1);
% limit the RIR by reflection order or by time-limit
type = 'maxTime';
% compute echogram due to pure propagation (frequency-independent)
for ns=1:nSrc
    for nr=1:nRec
        disp('')
        disp(['Compute echogram: Source ' num2str(ns) ' - Receiver ' num2str(nr)])
        % compute echogram
        echograms = ims_coreMtx(room, src(ns,:), rec(nr,:), type, max(limits));
        rec_echograms = rec_moduleSH(echograms, rec_orders);
        abs_echograms(ns, nr, :) = absorption_module(rec_echograms, abs_wall, limits);
    end
end
%disp('Apply SH directivites')
%rec_echograms = rec_moduleSH(echograms, rec_orders);
%clear echograms

% apply boundary absorption
% for ns = 1:nSrc
%     for nr=1:nRec
%         disp('')
%         disp(['Apply absorption: Source ' num2str(ns) ' - Receiver ' num2str(nr)])
%         % compute echogram
%         abs_echograms(ns, nr, :) = absorption_module(rec_echograms(ns,nr), abs_wall, limits);
%     end
% end
