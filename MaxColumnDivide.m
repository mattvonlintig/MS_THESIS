function [scaled_data, M] = MaxColumnDivide(data)
% MAXCOLUMNDIVIDE divides each column of data by the max in that column to
% produce a new dataset with columns of max == 1. For use to compare
% waveforms.

M = max(abs(data(:,:)));
scaled_data = data/(diag(M));
