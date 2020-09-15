function [pX_out,xsupport_out] = remove_redundant(pX, xsupport)
%REMOVE_REDUNDANT Delete low probability points and merge close points.
%   Detailed explanation goes here

%SETTINGS HERE
prob_thresh = 1e-10;
sup_thresh = 3e-3; %rel. to sigma

sigma = sqrt(pX.'*xsupport.^2); %std dev
sup_thresh = sup_thresh*sigma;

%discard low probability points
xsupport_temp = xsupport(pX>prob_thresh);
pX_temp = pX(pX>prob_thresh);

%combine points that are too close
pX_out = [];
xsupport_out = [];
merged_prob = 0; %combined prob of points to be merged
weighted_pos = 0; %of points to be merged
prev = xsupport_temp(1);
for i = 1:length(pX_temp)
    current = xsupport_temp(i);
    %if new point no longer close enough to be merged
    if abs(current-prev) > sup_thresh
        %store old cluster of points
        pX_out = [pX_out; merged_prob];
        xsupport_out = [xsupport_out; weighted_pos/merged_prob];
        %reset vars
        merged_prob = 0;
        weighted_pos = 0;
    end
    %accumulate/update
    merged_prob = merged_prob + pX_temp(i);
    weighted_pos = weighted_pos+current*pX_temp(i);
    prev = current;
end
pX_out = [pX_out; merged_prob];
xsupport_out = [xsupport_out; weighted_pos/merged_prob];

%some probability was lost, renormalize
pX_out = pX_out/sum(pX_out);
end

