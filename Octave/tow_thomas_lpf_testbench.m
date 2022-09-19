clear; clc; close all;

cfg.g_dc = 1;
cfg.f0 = 5;
cfg.Q = 1;
cfg.fs = 122.88*4;
cfg.rbw = 15/1000;
cfg.opt = "impinvar";

[y,outstruct] = tow_thomas_lpf(cfg,1);