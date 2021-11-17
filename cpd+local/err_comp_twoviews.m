function [errR, errT] = err_comp(grtR, grtT, R, t,s)

errR = 0;
errT = 0;
errR = errR + norm(grtR- R,'fro');
errT = errT + norm(grtT*s-t,2);

