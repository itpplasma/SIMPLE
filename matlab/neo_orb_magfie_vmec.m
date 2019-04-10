function [bmod, sqrtg, bder, hcovar, hctrvar, hcurl] = neo_orb_magfie_vmec(x)
% Fortran routine magfie_vmec in magfie.f90
    p_x = libpointer('doublePtr', x);

    bmod = 0.0;
    sqrtg = 0.0;
    bder = zeros(3);
    hcovar = zeros(3);
    hctrvar = zeros(3);
    hcurl = zeros(3);

    p_bmod = libpointer('doublePtr', bmod);
    p_sqrtg = libpointer('doublePtr', sqrtg);
    p_bder = libpointer('doublePtr', bder);
    p_hcovar = libpointer('doublePtr', hcovar);
    p_hctrvar = libpointer('doublePtr', hctrvar);
    p_hcurl = libpointer('doublePtr', hcurl);

    calllib('libneo_orb', 'magfie_vmec_', p_x, p_bmod, ...
      p_sqrtg, p_bder, p_hcovar, p_hctrvar, p_hcurl);
  
    bmod = p_bmod.Value;
    sqrtg = p_sqrtg.Value;
    bder  = p_bder.Value;
    hcovar = p_hcovar.Value;
    hctrvar = p_hctrvar.Value;
    hcurl = p_hcurl.Value;
  
end
