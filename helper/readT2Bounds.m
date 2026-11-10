function [bf, bs, info] = readT2Bounds(file, which)
%READT2BOUNDS  fit limits of a per-voxel biexponential fit in a T2estimate.m
%  [bf,bs,info] = readT2Bounds('T2estimate_scan402_418.m')            last per-voxel fit in the file
%  [bf,bs,info] = readT2Bounds('T2estimate_scan402_418.m','first')    first one ('last', 'first' or a number)
%  returns the T2f bounds bf = [lo hi] and the T2s bounds bs = [lo hi] (ms) of the fit that feeds
%  T2f(i,j,k) in a per-voxel loop, for either model used in this study:
%     c+a*((1-d)*exp(-x/b)+d*exp(-x/e))   Lower [0,bf1,c0,0,bs1]  Upper [Inf,bf2,Inf,1,bs2]
%     c+a*exp(-x/b)+d*exp(-x/e)            Lower [0,bf1,c0,0,bs1]  Upper [Inf,bf2,Inf,Inf,bs2]
%  T2estimate_scan402_418.m (2026-09-21) has TWO per-voxel loops: the first writes T2f / T2s / T2 (mask = allave > 12000),
%  the second writes T2allf / T2alls / T2all (mask = ROIall, the whole matrix) -> 'first' and 'last'.
%  Commented-out lines are ignored.  info.model = the fittype string, info.line = its line number, info.n = fits found.
    if nargin < 2 || isempty(which), which = 'last'; end
    txt = fileread(file);
    txt = regexprep(txt, '(?m)^[ \t]*%[^\n]*', '');          % drop commented lines (newlines kept)
    pat = ['''Lower'',\s*\[\s*0\s*,\s*([\d.]+)\s*,\s*[\d.]+\s*,\s*0\s*,\s*([\d.]+)\s*\][^\[]*' ...
           '''Upper'',\s*\[\s*Inf\s*,\s*([\d.]+)\s*,\s*Inf\s*,\s*(?:1|Inf)\s*,\s*([\d.]+)\s*\][^;]*;[^;]*' ...
           'fittype\(''(c\+a\*[^'']*)''[^\n]*\n(?:[^\n]*\n){0,3}\s*T2f\(i,j,k\)'];
    [tk, st] = regexp(txt, pat, 'tokens', 'start');
    if isempty(tk), error('readT2Bounds:notFound', 'no per-voxel biexponential fit found in %s', file); end
    if ischar(which) || isstring(which)
        switch lower(char(which)), case 'last', n = numel(tk); case 'first', n = 1; otherwise, error('which must be first, last or a number'); end
    else
        n = which;
    end
    t = tk{n};  v = str2double(t(1:4));
    bf = [v(1) v(3)];  bs = [v(2) v(4)];
    info.model = t{5};  info.line = 1 + nnz(txt(1:st(n)) == newline);  info.n = numel(tk);
end
