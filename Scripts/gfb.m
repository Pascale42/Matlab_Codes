% function getFileBase

function FileBase=gfb

currentDirectory = pwd;
[~, FileBase, ~] = fileparts(currentDirectory);