function [perm, signedCorr, Shat_aligned] = align_components( ...
                                      Strue, Shat, flipSigns)
%ALIGN_COMPONENTS  Align inferred components to ground-truth components.
%
% [perm, signedCorr, Shat_aligned] = ALIGN_COMPONENTS(Strue, Shat, flipSigns)
%
% Uses the Hungarian algorithm ( Statistics & ML Toolbox function
% MATCHPAIRS ) to find the one-to-one assignment that maximises the
% (absolute) Pearson correlation between the *ground-truth* components
% (columns of Strue) and the *inferred* components (columns of Shat).
%
% INPUT
%   Strue      –  (Ntrials × Ktrue)   matrix with the real components
%   Shat       –  (Ntrials × Khat )   matrix with the inferred components
%   flipSigns  –  logical, default true
%                If true, the function flips the sign of every matched
%                inferred component so that its correlation with the
%                corresponding ground-truth component is positive.
%
% OUTPUT
%   perm          (1 × Ktrue)   Index of the matched inferred component
%                               for every ground-truth component. 0 if the
%                               ground-truth component had no match
%                               (possible when Khat < Ktrue).
%   signedCorr    (1 × Ktrue)   Pearson correlation of every matched pair
%                               (after the optional sign flip).
%   Shat_aligned  (Ntrials × Khat)  Copy of Shat with the above sign flips
%                                   applied (empty if flipSigns == false).
%
% Example
%   [perm,corr,Shat2] = align_components(S,Shat);
%   disp([perm; corr])
%
% --------------------------------------------------------------------
%   Fabian-Sim 2025
% --------------------------------------------------------------------
arguments
    Strue   double {mustBeFinite, mustBeReal}
    Shat    double {mustBeFinite, mustBeReal}
    flipSigns logical = true
end

% ---------------------------------------------------------------------
% 1. Z-score columns so correlation is just column dot-product
% ---------------------------------------------------------------------
StrueZ = zscore(Strue ,0,1);      % (N × Ktrue)
ShatZ  = zscore(Shat  ,0,1);      % (N × Khat)

% ---------------------------------------------------------------------
% 2. Similarity (absolute Pearson correlation)
% ---------------------------------------------------------------------
sim    = abs(corr(StrueZ,ShatZ)); % Ktrue × Khat

% ---------------------------------------------------------------------
% 3. Hungarian assignment that MAXIMISES similarity
%    Use min penalty so every possible match is allowed
% ---------------------------------------------------------------------
penalty = min(sim(:)) - eps;

pairs  = matchpairs(sim, penalty, 'max');   % each row: [trueIdx  hatIdx]

Ktrue  = size(Strue,2);
perm         = zeros(1,Ktrue);         % 0 → unmatched
signedCorr   = zeros(1,Ktrue);

% ---------------------------------------------------------------------
% 4. Record permutation and sign information
% ---------------------------------------------------------------------
for r = 1:size(pairs,1)
    k  = pairs(r,1);        % ground-truth index
    j  = pairs(r,2);        % inferred index

    % signed correlation (before abs)
    thisCorr = corr(StrueZ(:,k), ShatZ(:,j));

    % Optional sign flip
    if flipSigns && thisCorr < 0
        Shat(:,j)  = -Shat(:,j);
        ShatZ(:,j) = -ShatZ(:,j);
        thisCorr   = -thisCorr;        % now positive
    end

    perm(k)       = j;
    signedCorr(k) = thisCorr;
end

if flipSigns
    Shat_aligned = Shat;
else
    Shat_aligned = [];
end
end
