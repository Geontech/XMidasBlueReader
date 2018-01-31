# XMidasBlueReader

This MATLAB class is a utility to progressively read through BLUE -format files.

```
bf = XMidasBlueReader('...path/to/bluefile');
% Header control block: bf.hcb

% Read 4096 samples into a cell array
samps = bf.read(4096);
sampsMat = cell2mat(samps); % As a vector or matrix
```