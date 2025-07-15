function posimat = figposi(numrow,numcol)
tempscr = get(0,'ScreenSize'); scrwid = tempscr(3);scrhei = tempscr(4); n_r_f = numcol; n_c_f = numrow; 
wid_f = fix(scrwid/n_r_f);hei_f = fix(scrhei/n_c_f);
fposix = 1:wid_f:scrwid; fposiy = scrhei-hei_f+1:-hei_f:1;
[fposiX,fposiY] = meshgrid(fposix,fposiy);
wid_f_mat = wid_f*ones(size(fposiX))*0.9; hei_f_mat = hei_f*ones(size(fposiX))*0.75; 
posimat = [fposiX(:) fposiY(:) wid_f_mat(:) hei_f_mat(:)];
end