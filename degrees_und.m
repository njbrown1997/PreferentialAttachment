function [deg] = degrees_und(CIJ)
%DEGREES_UND        Degree
%https://github.com/jblocher/matlab-network-utilities/blob/master/BrainConnectivity/degrees_und.m
%
%   deg = degrees_und(CIJ);
%
%   Node degree is the number of links connected to the node.
%
%   Input:      CIJ,    undirected (binary/weighted) connection matrix
%
%   Output:     deg,    node degree
%
%   Note: Weight information is discarded.
%
%
%   Olaf Sporns, Indiana University, 2002/2006/2008


% ensure CIJ is binary...
CIJ = double(CIJ~=0);
deg = sum(CIJ);