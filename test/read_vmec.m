function f=read_vmec(filename)
% READ_VMEC(filename) This function reads the VMEC wout file.
% This funciton reads the wout file and returns the data from the file
% in a structure.  Based on the 'readw_only_priv.f' subroutine.
% In addtion to the raw variables, the strcture also contains the
% reconstituted Fourier arrays in their full form.
% Currently this function can read VMEC files up to version 8+ (netCDF and
% text).  It returns fourier harmonics in the NESCOIL (nu+nv) format.
%
% Example usage
%      data=read_vmec('wout.test');     % Reads VMEC wout file
%      mdata=read_vmec('mercier.test'); % Reads VMEC mercier file
%
% Maintained by: Samuel Lazerson (lazerson@pppl.gov)
% Version:       1.96

% NOTES:
%   1/5/11      Modified output variables to use cos sine (c/s) notation
%               instead of the e notation.  Also added modifications to
%               support non-axisymmetric VMEC runs.
%   1/13/11     Overloaded to read vmec mercier file.
%   1/31/11     All quantities now mapped to full mesh.
%   2/01/11     Updated for version 8.47
%   2/28/11     Properly Handles half grid quantities (see wrfcn in pgplot)
%   5/31/11     Added support for +8.0 text files
%   3/21/12     Modified so opening netCDF via path is possible.
%               Fixed issue with mu constant when reading netCDF files.
%   1/30/13     Added calculation of chip
%   3/31/14     Fixed calculation of non-stellarator symmetric J terms.
%   1/10/14     Modified to read 8.51 files
%               Uses new methods for calculating J taking into acount the
%               odd modes properly.
%   3/1/16      Corrected variable names so text files use the new method
%               calculation of current densities.

% Defaults
netcdffile=0;
version=999999.0;
% Number of Arguments
if (nargin == 0)
    disp('read_vmec requires a filename.');
    f=1;
    return
end
% Check to see the file exists
fid=fopen(filename,'r+');
if (fid < 0)
    disp(['ERROR: Could not find file ' filename]);
    f=fid;
    return
else
    fclose(fid);
end
% Handle File type
nfilename=size(filename,2);
if strcmp(filename(nfilename-2:nfilename),'.nc') && ...
        ~isempty(strfind(filename,'wout'))
    netcdffile=1; % This is a netCDF file
elseif strcmp(filename(1:7),'mercier')
    disp('mercier file')
    f=read_vmec_mercier(filename);
    return
elseif strcmp(filename(1:6),'jxbout')
    disp('jxbout file')
    f=read_vmec_jxbout(filename);
    return
else
    fid=fopen(filename,'r'); % Open File
    % Read First Line and extract version information
    line=fgetl(fid);
    index=find(line == '=')+1;
    length=size(line,2);
    version=str2num(strtrim(line(index:length)));
    % Handle unknown versions
    if (version < 0) || (version > 8.52)
        disp('Unknown file type or version.');
        f=2;
        return
    end
    % VMEC files with comma delimited values do exist, in an attempt to
    % handle them we dynamically create the format specifier for fscanf.
    start=ftell(fid);  % Get the current poisition to rewind to
    line=fgetl(fid);  % Get the next line
    if isempty(strfind(line,','))
        fmt='%g';
    else
        fmt='%g,';
        disp('Comma delimited file detected!');
        disp('     Only 6.54+ version supported!');
        if version <6.54
            return
        end
    end
    fseek(fid,start,'bof'); % Go back to just after the version.
end
% Handle versions
%   pre 5.1, 6.05, 6.20, 6.50, 6.54
if (version <= 5.10) && ~netcdffile
    f=read_vmec_orig(fid,fmt);
elseif (version <= 6.05) && ~netcdffile
    f=read_vmec_605(fid,fmt);
elseif (version <= 6.20) && ~netcdffile
    f=read_vmec_620(fid,fmt);
elseif (version <= 6.50) && ~netcdffile
    f=read_vmec_650(fid,fmt);
elseif (version <= 6.95) && ~netcdffile
    f=read_vmec_695(fid,fmt);
elseif (version <= 8.00) && ~netcdffile
    f=read_vmec_800(fid,fmt);
elseif (version <= 8.52) && ~netcdffile
    f=read_vmec_847(fid,fmt);
elseif netcdffile
    f=read_vmec_netcdf(filename);
    version=f.version;
end
if ~netcdffile
    fclose(fid); % Close the File
    f.version=version; % Add Version info to structure
end
% If VMEC through an error then return what was read
if (f.ierr_vmec && (f.ierr_vmec ~= 4))
    disp('VMEC runtime error detected!');
    return
end
if (~isfield(f,'lrfplogical'))
    f.lrfplogical=0;
end
% Need to convert various quantities to full mesh
f=half2fullmesh(f);
% Now recompose the Fourier arrays
msize=max(f.xm);
nsize=max(f.xn/f.nfp)-min(f.xn/f.nfp)+1;
f.rbc=zeros(msize+1,nsize,f.ns);
f.zbs=zeros(msize+1,nsize,f.ns);
f.lbs=zeros(msize+1,nsize,f.ns);
% Create derivative terms
% B_R   = B^u*(dR/du) + B^v (dR/dv)
% B_phi = R*B^v
% B_Z   = B^u*(dZ/du) + B^v (dZ/dv)
f.rsc=zeros(msize+1,nsize,f.ns); % dRmn/ds*cos
f.rus=zeros(msize+1,nsize,f.ns); % dRmn/du*sin
f.rvs=zeros(msize+1,nsize,f.ns); % dRmn/dv*sin
f.zss=zeros(msize+1,nsize,f.ns); % dZmn/ds*sin
f.zuc=zeros(msize+1,nsize,f.ns); % dZmn/du*cos
f.zvc=zeros(msize+1,nsize,f.ns); % dZmn/dv*cos
offset=min(f.xn./f.nfp)-1;
f.xn=-f.xn;
for i=1:f.ns
    for j=1:f.mnmax
        m=f.xm(j)+1;
        n=-offset+f.xn(j)./f.nfp;
        f.rbc(m,n,i)    =   f.rmnc(j,i);
        f.rus(m,n,i)    =  -f.rmnc(j,i)     .*  f.xm(j);
        f.rvs(m,n,i)    =  -f.rmnc(j,i)     .*  f.xn(j);
        f.zbs(m,n,i)    =   f.zmns(j,i);
        f.zuc(m,n,i)    =   f.zmns(j,i)     .*  f.xm(j);
        f.zvc(m,n,i)    =   f.zmns(j,i)     .*  f.xn(j);
        f.lbs(m,n,i)    =   f.lmns(j,i);
    end
end
f.rumns=-f.rmnc.*repmat(f.xm',[1 f.ns]);
f.rvmns=-f.rmnc.*repmat(f.xn',[1 f.ns]);
f.zumnc=f.zmns.*repmat(f.xm',[1 f.ns]);
f.zvmnc=f.zmns.*repmat(f.xn',[1 f.ns]);
% Handle Radial Derivatives
f.rsc=f.rbc;
f.zss=f.zbs;
f.rsmnc=f.rmnc;
f.zsmns=f.zmns;
for i=2:f.ns-1
    f.rsmnc(:,i)=f.rmnc(:,i+1)-f.rmnc(:,i-1);
    f.zsmns(:,i)=f.zmns(:,i+1)-f.zmns(:,i-1);
    f.rsc(:,:,i)=f.rbc(:,:,i+1)-f.rbc(:,:,i-1);
    f.zss(:,:,i)=f.zbs(:,:,i+1)-f.zbs(:,:,i-1);
end
f.rsc=0.5.*f.rsc;
f.zss=0.5.*f.zss;
f.rsmnc=0.5*f.rsmnc;
f.zsmns=0.5*f.zsmns;
f.rsc(:,:,1)=f.rbc(:,:,2)-2.*f.rbc(:,:,1);
f.zss(:,:,1)=f.zbs(:,:,2)-2.*f.zbs(:,:,1);
f.rsmnc(:,1)=f.rsmnc(:,2)-2.*f.rsmnc(:,1);
f.zsmns(:,1)=f.zsmns(:,2)-2.*f.zsmns(:,1);
f.rsc(:,:,f.ns)=2.*f.rbc(:,:,f.ns)-f.rbc(:,:,f.ns-1);
f.zss(:,:,f.ns)=2.*f.zbs(:,:,f.ns)-f.zbs(:,:,f.ns-1);
f.rsmnc(:,f.ns)=2.*f.rsmnc(:,f.ns)-f.rsmnc(:,f.ns-1);
f.zsmns(:,f.ns)=2.*f.zsmns(:,f.ns)-f.zsmns(:,f.ns-1);
% Handle Vector values seperately to deal with nyqyist
if isfield(f,'xm_nyq')
    msize=max(f.xm_nyq);
    nsize=max(f.xn_nyq/f.nfp)-min(f.xn_nyq/f.nfp)+1;
    mnmax=f.mnmax_nyq;
    offset=min(f.xn_nyq./f.nfp)-1;
    f.xn_nyq=-f.xn_nyq;
    xn=f.xn_nyq;
    xm=f.xm_nyq;
    f.mpol=max(f.xm_nyq);  % Do this for VMECplot
    f.ntor=max(f.xn_nyq/f.nfp);  % Do this for VMECplot
else
    msize=max(f.xm);
    nsize=max(f.xn/f.nfp)-min(f.xn/f.nfp)+1;
    mnmax=f.mnmax;
    offset=min(f.xn./f.nfp)-1;
    xn=f.xn;
    xm=f.xm;
end
f.bc=zeros(msize+1,nsize,f.ns);
f.gc=zeros(msize+1,nsize,f.ns);
f.b_ss=zeros(msize+1,nsize,f.ns);
f.b_uc=zeros(msize+1,nsize,f.ns);
f.b_vc=zeros(msize+1,nsize,f.ns);
f.buc=zeros(msize+1,nsize,f.ns);
f.bvc=zeros(msize+1,nsize,f.ns);
f.currvc=zeros(msize+1,nsize,f.ns);
if f.version > 8.0, f.curruc=zeros(msize,nsize,f.ns); end
for i=1:f.ns
    for j=1:mnmax
        m=xm(j)+1;
        n=-offset+xn(j)./f.nfp;
        f.bc(m,n,i)     =   f.bmnc(j,i);
        f.gc(m,n,i)     =   f.gmnc(j,i);
        f.b_ss(m,n,i)   =   f.bsubsmns(j,i);
        f.b_uc(m,n,i)   =   f.bsubumnc(j,i);
        f.b_vc(m,n,i)   =   f.bsubvmnc(j,i);
        f.buc(m,n,i)    =   f.bsupumnc(j,i);
        f.bvc(m,n,i)    =   f.bsupvmnc(j,i);
        f.currvc(m,n,i) =   f.currvmnc(j,i);
        if f.version > 8.0
            f.curruc(m,n,i) =   f.currumnc(j,i);
        end
    end
end
% Note derivative terms not implemented for non-axisymmetric runs
% Handle non-axisymmetric runs
if f.iasym==1
    msize=max(f.xm);
    nsize=max(f.xn/f.nfp)-min(f.xn/f.nfp)+1;
    f.rbs=zeros(msize+1,nsize,f.ns);
    f.zbc=zeros(msize+1,nsize,f.ns);
    f.lbc=zeros(msize+1,nsize,f.ns);
    % Create derivative terms
    % B_R   = B^u*(dR/du) + B^v (dR/dv)
    % B_phi = R*B^v
    % B_Z   = B^u*(dZ/du) + B^v (dZ/dv)
    f.rss=zeros(msize+1,nsize,f.ns); % dRmn/ds*cos
    f.ruc=zeros(msize+1,nsize,f.ns); % dRmn/du*sin
    f.rvc=zeros(msize+1,nsize,f.ns); % dRmn/dv*sin
    f.zsc=zeros(msize+1,nsize,f.ns); % dZmn/ds*sin
    f.zus=zeros(msize+1,nsize,f.ns); % dZmn/du*cos
    f.zvs=zeros(msize+1,nsize,f.ns); % dZmn/dv*cos
    offset=min(f.xn./f.nfp)-1;
    for i=1:f.ns
        for j=1:f.mnmax
            m=f.xm(j)+1;
            n=-offset+f.xn(j)./f.nfp;
            f.rbs(m,n,i)    =   f.rmns(j,i);
            f.ruc(m,n,i)    =  -f.rmns(j,i)     .*  f.xm(j);
            f.rvc(m,n,i)    =  -f.rmns(j,i)     .*  f.xn(j);
            f.zbc(m,n,i)    =   f.zmnc(j,i);
            f.zus(m,n,i)    =   f.zmnc(j,i)     .*  f.xm(j);
            f.zvs(m,n,i)    =   f.zmnc(j,i)     .*  f.xn(j);
            f.lbc(m,n,i)    =   f.lmnc(j,i);
        end
    end
    f.rumnc=-f.rmns.*repmat(f.xm',[1 f.ns]);
    f.rvmnc=-f.rmns.*repmat(f.xn',[1 f.ns]);
    f.zumns=f.zmnc.*repmat(f.xm',[1 f.ns]);
    f.zvmns=f.zmnc.*repmat(f.xn',[1 f.ns]);
    % Handle Radial Derivatives
    f.rss=f.rbs;
    f.zsc=f.zbc;
    f.rsmns=f.rmns;
    f.zsmnc=f.zmnc;
    for i=2:f.ns-1
        f.rsmns(:,i)=f.rmns(:,i+1)-f.rmns(:,i-1);
        f.zsmnc(:,i)=f.zmnc(:,i+1)-f.zmnc(:,i-1);
        f.rss(:,:,i)=f.rbs(:,:,i+1)-f.rbs(:,:,i-1);
        f.zsc(:,:,i)=f.zbc(:,:,i+1)-f.zbc(:,:,i-1);
    end
    f.rss=0.5.*f.rss;
    f.zsc=0.5.*f.zsc;
    f.rsmns=0.5*f.rsmns;
    f.zsmnc=0.5*f.zsmnc;
    f.rss(:,:,1)=f.rbs(:,:,2)-2.*f.rbs(:,:,1);
    f.zsc(:,:,1)=f.zbc(:,:,2)-2.*f.zbc(:,:,1);
    f.rsmns(:,1)=f.rsmns(:,2)-2.*f.rsmns(:,1);
    f.zsmnc(:,1)=f.zsmnc(:,2)-2.*f.zsmnc(:,1);
    f.rss(:,:,f.ns)=2.*f.rbs(:,:,f.ns)-f.rbs(:,:,f.ns-1);
    f.zsc(:,:,f.ns)=2.*f.zbc(:,:,f.ns)-f.zbc(:,:,f.ns-1);
    f.rsmns(:,f.ns)=2.*f.rsmns(:,f.ns)-f.rsmns(:,f.ns-1);
    f.zsmnc(:,f.ns)=2.*f.zsmnc(:,f.ns)-f.zsmnc(:,f.ns-1);
    % Handle Vector values seperately to deal with nyqyist
    if isfield(f,'xm_nyq')
        msize=max(f.xm_nyq);
        nsize=max(f.xn_nyq/f.nfp)-min(f.xn_nyq/f.nfp)+1;
        mnmax=f.mnmax_nyq;
        offset=min(f.xn_nyq./f.nfp)-1;
        xn=f.xn_nyq;
        xm=f.xm_nyq;
    else
        msize=max(f.xm);
        nsize=max(f.xn/f.nfp)-min(f.xn/f.nfp)+1;
        mnmax=f.mnmax;
        offset=min(f.xn./f.nfp)-1;
        xn=f.xn;
        xm=f.xm;
    end
    f.bs=zeros(msize+1,nsize,f.ns);
    f.gs=zeros(msize+1,nsize,f.ns);
    f.b_sc=zeros(msize+1,nsize,f.ns);
    f.b_us=zeros(msize+1,nsize,f.ns);
    f.b_vs=zeros(msize+1,nsize,f.ns);
    f.bus=zeros(msize+1,nsize,f.ns);
    f.bvs=zeros(msize+1,nsize,f.ns);
    f.currvs=zeros(msize+1,nsize,f.ns);
    if f.version > 8.0, f.currus=zeros(msize,nsize,f.ns); end
    for i=1:f.ns
        for j=1:mnmax
            m=xm(j)+1;
            n=-offset+xn(j)./f.nfp;
            f.bs(m,n,i)     =   f.bmns(j,i);
            f.gs(m,n,i)     =   f.gmns(j,i);
            f.b_sc(m,n,i)   =   f.bsubsmnc(j,i);
            f.b_us(m,n,i)   =   f.bsubumns(j,i);
            f.b_vs(m,n,i)   =   f.bsubvmns(j,i);
            f.bus(m,n,i)    =   f.bsupumns(j,i);
            f.bvs(m,n,i)    =   f.bsupvmns(j,i);
            f.currvs(m,n,i) =   f.currvmns(j,i);
            if f.version > 8.0
                f.currus(m,n,i) =   f.currumnc(j,i);
            end
        end
    end
end
% Calculate chi and chip
if (f.lrfplogical)
    %f.chipf=ones(1,f.ns);
else
    f.chipf=f.iotaf.*f.phipf;
    f.chif=cumsum(f.chipf);
end
% Create the Resonance Array
f.M=1:max(f.xm);
f.N=min(f.xn(f.xn > 0)):min(f.xn(f.xn > 0)):max(f.xn);
Minv=1./f.M;
f.iota_res=f.N'*Minv;
%iota_spline=pchip(1:f.ns,abs(f.iotaf));
%f.s_res=0.0.*f.iota_res;
%for i=1:size(f.s_res,1)
%    for j=1:size(f.s_res,2)
%        if (f.iota_res(i,j) > min(abs(f.iotaf))) && (f.iota_res(i,j) < max(abs(f.iotaf)))
%            f_iota=@(x)ppval(iota_spline,x)-f.iota_res(i,j);
%            f.s_res(i,j)=fzero(f_iota,50);
%        end
%    end
%end
% Fix the flow variables
if isfield(f,'protmnc'), f.protmnc = f.protmnc ./ (pi*4.0E-7); end;
if isfield(f,'protrsqmnc'), f.protrsqmnc = f.protrsqmnc ./ (pi*4.0E-7); end;
if isfield(f,'prprmnc'), f.prprmnc = f.prprmnc ./ (pi*4.0E-7); end;
if isfield(f,'protmns'), f.protmns = f.protmns ./ (pi*4.0E-7); end;
if isfield(f,'protrsqmns'), f.protrsqmns = f.protrsqmns ./ (pi*4.0E-7); end;
if isfield(f,'prprmns'), f.prprmns = f.prprmns ./ (pi*4.0E-7); end;
% Calculate the stored energy
f.eplasma=1.5*pi*pi*sum(f.vp.*f.presf)./f.ns;
% Set some defaults
f.nu=2.*(f.mpol+1)+6;
f.datatype='wout';
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=read_vmec_orig(fid,fmt)
% For 5.10 and earlier

% Unit Conversions
dmu0=(2.0e-7)/pi;
% Read Data
data=fscanf(fid,'%g',13);
f.wb=data(1);
f.wp=data(2);
f.gamma=data(3);
f.pfac=data(4);
f.nfp=data(5);
f.ns=data(6);
f.mpol=data(7);
f.ntor=data(8);
f.mnmax=data(9);
f.itfsq=data(10);
f.niter=data(11);
f.iasym=data(12);
f.ireconstruct(13);
data=fscanf(fid,'%d',5);
f.imse=data(1);
f.itse=data(2);
f.nbsets=data(3);
f.nobd=data(4);
f.nextcur=data(5);
f.nstore_seq=100;
% Error Check
if (f.ierr_vmec && (f.ierr_vmec ~= 4))
    disp('ierr_vmec >0');
    return
end
% Read nbfld
if (f.nbsets > 0)
    f.nbfld=fscanf(fid,'%g',f.nbsets);
end
% Read mgrid filename
if isempty(strfind(fmt,','))
    f.mgrid_file=fscanf(fid,'%s',1);
    fmt2='%g%g';
    fmt3='%g%g%g';
    fmt6='%g%g%g%g%g%g';
    fmt7='%g%g%g%g%g%g%g';
    fmt11='%g%g%g%g%g%g%g%g%g%g%g';
    fmt12='%g%g%g%g%g%g%g%g%g%g%g%g';
    fmt14='%g%g%g%g%g%g%g%g%g%g%g%g%g%g';
    fmt2_11='%d%d%g%g%g%g%g%g%g%g%g%g%g';
    fmt2_14='%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g';
else
    f.mgrid_file=fscanf(fid,'%s,',1);
    fmt2='%g,%g';
    fmt3='%g,%g,%g';
    fmt6='%g,%g,%g,%g,%g,%g';
    fmt7='%g,%g,%g,%g,%g,%g,%g';
    fmt11='%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt12='%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt14='%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt2_11='%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt2_14='%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
end
% Read Arrays
if f.iasym > 0
    data1=fscanf(fid,fmt2_14,[16 f.mnmax]);
    data=fscanf(fid,fmt14,[14 f.mnmax*(f.ns-1)]);
else
    data1=fscanf(fid,fmt2_11,[13 f.mnmax]);
    data=fscanf(fid,fmt11,[11 f.mnmax*(f.ns-1)]);
end
% Extract Data from Arrays
f.xm=data1(1,:);
f.xn=data1(2,:);
% First reshape data
f.rmnc=[data1(3,:)' reshape(data(1,:),f.mnmax,f.ns-1)];
f.zmns=[data1(4,:)' reshape(data(2,:),f.mnmax,f.ns-1)];
f.lmns=[data1(5,:)' reshape(data(3,:),f.mnmax,f.ns-1)];
f.bmn=[data1(6,:)' reshape(data(4,:),f.mnmax,f.ns-1)];
f.gmn=[data1(7,:)' reshape(data(5,:),f.mnmax,f.ns-1)];
f.bsubumn=[data1(8,:)' reshape(data(6,:),f.mnmax,f.ns-1)];
f.bsubvmn=[data1(9,:)' reshape(data(7,:),f.mnmax,f.ns-1)];
f.bsubsmn=[data1(10,:)' reshape(data(8,:),f.mnmax,f.ns-1)];
f.bsupumn=[data1(11,:)' reshape(data(9,:),f.mnmax,f.ns-1)];
f.bsupvmn=[data1(12,:)' reshape(data(10,:),f.mnmax,f.ns-1)];
f.currvmn=[data1(13,:)' reshape(data(11,:),f.mnmax,f.ns-1)];
% Read the half-mesh quantities
data=fscanf(fid,fmt12,[12 f.ns/2]);
f.iotas=data(1,:);
f.mass=data(2,:);
f.pres=data(3,:);
f.phip=data(4,:);
f.buco=data(5,:);
f.bvco=data(6,:);
f.phi=data(7,:);
f.vp=data(8,:);
f.overr=data(9,:);
f.jcuru=data(10,:);
f.jcurv=data(11,:);
f.specw=data(12,:);
data=fscanf(fid,fmt6,6);
f.aspect=data(1);
f.betatot=data(2);
f.betapol=data(3);
f.betator=data(4);
f.betaxis=data(5);
f.b0=data(6);
% Mercier Criterion
data=fscanf(fid,fmt6,[6 f.ns-2]);
f.Dmerc=data(1,:);
f.Dshear=data(2,:);
f.Dwell=data(3,:);
f.Dcurr=data(4,:);
f.Dgeod=data(5,:);
f.equif=data(6,:);
if (f.nextcur > 0)
    f.extcur=fscanf(fid,fmt,f.nextcur);
    f.curlabel=fscanf(fid,fmt,f.nextcur);
end
data=fscanf(fid,fmt2,[2 f.nstore_seq]);
f.sqt=data(1,:);
f.wdot=data(1,:);
% Convert from Internal Units to Physical Units
f.mass=f.mass./dmu0;
f.pres=f.pres./dmu0;
f.jcuru=f.jcuru./dmu0;
f.jcurv=f.jcurv./dmu0;
f.jdotb=f.jdotb./dmu0;
f.phi=-f.phi;% Data and MSE Fits
if (f.ireconstruct > 0)
    if ((f.imse >= 2) || (f.itse >0))
        f.twsgt=fscanf(fid,fmt,1);
        f.msewgt=fscanf(fid,fmt,1);
        f.isnodes=fscanf(fid,fmt,1);
        data=fscanf(fid,fmt3,[3 f.isnodes]);
        f.sknots=data(1,:);
        f.ystark=data(2,:);
        f.y2stark=data(3,:);
        f.ipnodes=fscanf(fid,fmt,1);
        data=fscanf(fid,fmt3,[3 f.ipnodes]);
        f.pknots=data(1,:);
        f.ythom=data(2,:);
        f.y2thom=data(3,:);
        data=fscanf(fid,fmt7,[7 (2*f.ns)-1]);
        f.anglemse=data(1,:);
        f.rmid=data(2,:);
        f.qmid=data(3,:);
        f.shear=data(4,:);
        f.presmid=data(5,:);
        f.alfa=data(6,:);
        f.curmid=data(7,:);
        data=fscanf(fid,fmt3,[3 f.imse]);
        f.rstark=data(1,:);
        f.datastark=data(2,:);
        f.qmeas=data(3,:);
        data=fscanf(fid,fmt2,[2 f.itse]);
        f.rthom=data(1,:);
        f.datathom=data(2,:);
    end
    if (f.nobd > 0)
        data=fscanf(fid,fmt3,[3,f.nobd]);
        f.dsiext=data(1,:);
        f.plflux=data(2,:);
        f.dsiobt=data(3,:);
        f.flmwgt=fscanf(fid,fmt,1);
    end
    nbfldn=sum(nbldf(1:f.nbsets));
    if (nbfldn > 0)
        for n=1:nbsets
            data=fscanf(fid,fmt3,[3 f.nbfld(n)]);
            f.bcoil(:,n)=data(1,:);
            f.plbfld(:,n)=data(2,:);
            f.bbc(:,n)=data(3,:);
        end
        f.bcwgt=fscanf(fid,fmt,1);
    end
    f.phidiam=fscanf(fid,fmt,1);
    f.delphid=fscanf(fid,fmt,1);
    % Read Limiter and Prout Plotting Specs
    f.nsets=fscanf(fid,fmt,1);
    f.nparts=fscanf(fid,fmt,1);
    f.nlim=fscanf(fid,fmt,1);
    f.nsetsn=fscanf(fid,fmt,f.nsets);
    f.pfcspec=zeros(f.nparts,max(f.nxsetsn),f.nsets);
    for k=1:f.nsets
        for j=1:f.nsetsn(k)
            for i=1:f.nparts
                f.pfcspec(i,j,k)=fscanf(fid,fmt,1);
            end
        end
    end
    f.limitr=fscanf(fid,fmt,f.nlim);
    f.rlim=zeros(max(f.limitr),f.nlim);
    f.zlim=zeros(max(f.limitr),f.nlim);
    for j=1:f.nlim
        for i=1:f.limitr(j)
            data=fscanf(fid,fmt2,2);
            f.rlim(i,j)=data(1);
            f.zlim(i,j)=data(2);
        end
    end
    f.nrgrid=fscanf(fid,fmt,1);
    f.nzgrid=fscanf(fid,fmt,1);
    f.tokid=fscanf(fid,fmt,1);
    f.rx1=fscanf(fid,fmt,1);
    f.rx2=fscanf(fid,fmt,1);
    f.zy1=fscanf(fid,fmt,1);
    f.zy2=fscanf(fid,fmt,1);
    f.conif=fscanf(fid,fmt,1);
    f.imatch_phiedge=fscanf(fid,fmt,1);
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=read_vmec_605(fid,fmt)
% For 6.05

% Unit Conversions
dmu0=(2.0e-7)/pi;
% Read Data
data=fscanf(fid,'%g',6);
f.wb=data(1);
f.wp=data(2);
f.gamma=data(3);
f.pfac=data(4);
f.rmax_surf=data(5);
f.rmin_surf=data(6);
data=fscanf(fid,'%d',10);
f.nfp=data(1);
f.ns=data(2);
f.mpol=data(3);
f.ntor=data(4);
f.mnmax=data(5);
f.itfsq=data(6);
f.niter=data(7);
f.iasym=data(8);
f.ireconstruct=data(9);
f.ierr_vmec=data(10);
data=fscanf(fid,'%d',5);
f.imse=data(1);
f.itse=data(2);
f.nbsets=data(3);
f.nobd=data(4);
f.nextcur=data(5);
f.nstore_seq=100;
% Error Check
if (f.ierr_vmec && (f.ierr_vmec ~= 4))
    disp('ierr_vmec >0');
    return
end
% Read nbfld
if (f.nbsets > 0)
    f.nbfld=fscanf(fid,'%g',f.nbsets);
end
% Read mgrid filename
f.mgrid_file=fscanf(fid,'%s',1);
% Read Arrays
if f.iasym > 0
    data1=fscanf(fid,'%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g',[16 f.mnmax]);
    data=fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%g%g%g%g',[14 f.mnmax*(f.ns-1)]);
else
    data1=fscanf(fid,'%d%d%g%g%g%g%g%g%g%g%g%g%g',[13 f.mnmax]);
    data=fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%g',[11 f.mnmax*(f.ns-1)]);
end
% Extract Data from Arrays
f.xm=data1(1,:);
f.xn=data1(2,:);
% First reshape data
f.rmnc=[data1(3,:)' reshape(data(1,:),f.mnmax,f.ns-1)];
f.zmns=[data1(4,:)' reshape(data(2,:),f.mnmax,f.ns-1)];
f.lmns=[data1(5,:)' reshape(data(3,:),f.mnmax,f.ns-1)];
f.bmnc=[data1(6,:)' reshape(data(4,:),f.mnmax,f.ns-1)];
f.gmnc=[data1(7,:)' reshape(data(5,:),f.mnmax,f.ns-1)];
f.bsubumnc=[data1(8,:)' reshape(data(6,:),f.mnmax,f.ns-1)];
f.bsubvmnc=[data1(9,:)' reshape(data(7,:),f.mnmax,f.ns-1)];
f.bsubsmns=[data1(10,:)' reshape(data(8,:),f.mnmax,f.ns-1)];
f.bsupumnc=[data1(11,:)' reshape(data(9,:),f.mnmax,f.ns-1)];
f.bsupvmnc=[data1(12,:)' reshape(data(10,:),f.mnmax,f.ns-1)];
f.currvmnc=[data1(13,:)' reshape(data(11,:),f.mnmax,f.ns-1)];
% Read the half-mesh quantities
data=fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%g%g',[12 f.ns/2]);
f.iotas=data(1,:);
f.mass=data(2,:);
f.pres=data(3,:);
f.phip=data(4,:);
f.buco=data(5,:);
f.bvco=data(6,:);
f.phi=data(7,:);
f.vp=data(8,:);
f.overr=data(9,:);
f.jcuru=data(10,:);
f.jcurv=data(11,:);
f.specw=data(12,:);
data=fscanf(fid,'%g%g%g%g%g%g',6);
f.aspect=data(1);
f.betatot=data(2);
f.betapol=data(3);
f.betator=data(4);
f.betaxis=data(5);
f.b0=data(6);
% Mercier Criterion
data=fscanf(fid,'%g%g%g%g%g%g',[6 f.ns-2]);
f.Dmerc=data(1,:);
f.Dshear=data(2,:);
f.Dwell=data(3,:);
f.Dcurr=data(4,:);
f.Dgeod=data(5,:);
f.equif=data(6,:);
if (f.nextcur > 0)
    f.extcur=fscanf(fid,'%g',f.nextcur);
    f.curlabel=fscanf(fid,'%g',f.nextcur);
end
data=fscanf(fid,'%g%g',[2 f.nstore_seq]);
f.sqt=data(1,:);
f.wdot=data(1,:);
% Convert from Internal Units to Physical Units
f.mass=f.mass./dmu0;
f.pres=f.pres./dmu0;
f.jcuru=f.jcuru./dmu0;
f.jcurv=f.jcurv./dmu0;
f.jdotb=f.jdotb./dmu0;
f.phi=-f.phi;% Data and MSE Fits
if (f.ireconstruct > 0)
    if ((f.imse >= 2) || (f.itse >0))
        f.twsgt=fscanf(fid,'%g',1);
        f.msewgt=fscanf(fid,'%g',1);
        f.isnodes=fscanf(fid,'%d',1);
        data=fscanf(fid,'%g%g%g',[3 f.isnodes]);
        f.sknots=data(1,:);
        f.ystark=data(2,:);
        f.y2stark=data(3,:);
        f.ipnodes=fscanf(fid,'%d',1);
        data=fscanf(fid,'%g%g%g',[3 f.ipnodes]);
        f.pknots=data(1,:);
        f.ythom=data(2,:);
        f.y2thom=data(3,:);
        data=fscanf(fid,'%g%g%g%g%g%g%g',[7 (2*f.ns)-1]);
        f.anglemse=data(1,:);
        f.rmid=data(2,:);
        f.qmid=data(3,:);
        f.shear=data(4,:);
        f.presmid=data(5,:);
        f.alfa=data(6,:);
        f.curmid=data(7,:);
        data=fscanf(fid,'%g%g%g',[3 f.imse]);
        f.rstark=data(1,:);
        f.datastark=data(2,:);
        f.qmeas=data(3,:);
        data=fscanf(fid,'%g%g',[2 f.itse]);
        f.rthom=data(1,:);
        f.datathom=data(2,:);
    end
    if (f.nobd > 0)
        data=fscanf(fid,'%g%g%g',[3,f.nobd]);
        f.dsiext=data(1,:);
        f.plflux=data(2,:);
        f.dsiobt=data(3,:);
        f.flmwgt=fscanf(fid,'%g',1);
    end
    nbfldn=sum(nbldf(1:f.nbsets));
    if (nbfldn > 0)
        for n=1:nbsets
            data=fscanf(fid,'%g%g%g',[3 f.nbfld(n)]);
            f.bcoil(:,n)=data(1,:);
            f.plbfld(:,n)=data(2,:);
            f.bbc(:,n)=data(3,:);
        end
        f.bcwgt=fscanf(fid,'%g',1);
    end
    f.phidiam=fscanf(fid,'%g',1);
    f.delphid=fscanf(fid,'%g',1);
    % Read Limiter and Prout Plotting Specs
    f.nsets=fscanf(fid,'%g',1);
    f.nparts=fscanf(fid,'%g',1);
    f.nlim=fscanf(fid,'%g',1);
    f.nsetsn=fscanf(fid,'%g',f.nsets);
    f.pfcspec=zeros(f.nparts,max(f.nxsetsn),f.nsets);
    for k=1:f.nsets
        for j=1:f.nsetsn(k)
            for i=1:f.nparts
                f.pfcspec(i,j,k)=fscanf(fid,'%g',1);
            end
        end
    end
    f.limitr=fscanf(fid,'%g',f.nlim);
    f.rlim=zeros(max(f.limitr),f.nlim);
    f.zlim=zeros(max(f.limitr),f.nlim);
    for j=1:f.nlim
        for i=1:f.limitr(j)
            data=fscanf(fid,'%g%g',2);
            f.rlim(i,j)=data(1);
            f.zlim(i,j)=data(2);
        end
    end
    f.nrgrid=fscanf(fid,'%g',1);
    f.nzgrid=fscanf(fid,'%g',1);
    f.tokid=fscanf(fid,'%g',1);
    f.rx1=fscanf(fid,'%g',1);
    f.rx2=fscanf(fid,'%g',1);
    f.zy1=fscanf(fid,'%g',1);
    f.zy2=fscanf(fid,'%g',1);
    f.conif=fscanf(fid,'%g',1);
    f.imatch_phiedge=fscanf(fid,'%g',1);
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=read_vmec_620(fid,fmt)
% For 6.20

% Unit Conversions
dmu0=(2.0e-7)/pi;
% Read Data
data=fscanf(fid,'%g',6);
f.wb=data(1);
f.wp=data(2);
f.gamma=data(3);
f.pfac=data(4);
f.rmax_surf=data(5);
f.rmin_surf=data(6);
data=fscanf(fid,'%d',10);
f.nfp=data(1);
f.ns=data(2);
f.mpol=data(3);
f.ntor=data(4);
f.mnmax=data(5);
f.itfsq=data(6);
f.niter=data(7);
f.iasym=data(8);
f.ireconstruct=data(9);
f.ierr_vmec=data(10);
data=fscanf(fid,'%d',5);
f.imse=data(1);
f.itse=data(2);
f.nbsets=data(3);
f.nobd=data(4);
f.nextcur=data(5);
f.nstore_seq=100;
% Error Check
if (f.ierr_vmec && (f.ierr_vmec ~= 4))
    disp('ierr_vmec >0');
    return
end
% Read nbfld
if (f.nbsets > 0)
    f.nbfld=fscanf(fid,'%g',f.nbsets);
end
% Read mgrid filename
f.mgrid_file=fscanf(fid,'%s',1);
% Read Arrays
if f.iasym > 0
    data1=fscanf(fid,'%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g',[16 f.mnmax]);
    data=fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%g%g%g%g',[14 f.mnmax*(f.ns-1)]);
else
    data1=fscanf(fid,'%d%d%g%g%g%g%g%g%g%g%g%g%g',[13 f.mnmax]);
    data=fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%g',[11 f.mnmax*(f.ns-1)]);
end
% Extract Data from Arrays
f.xm=data1(1,:);
f.xn=data1(2,:);
% First reshape data
f.rmnc=[data1(3,:)' reshape(data(1,:),f.mnmax,f.ns-1)];
f.zmns=[data1(4,:)' reshape(data(2,:),f.mnmax,f.ns-1)];
f.lmns=[data1(5,:)' reshape(data(3,:),f.mnmax,f.ns-1)];
f.bmn=[data1(6,:)' reshape(data(4,:),f.mnmax,f.ns-1)];
f.gmn=[data1(7,:)' reshape(data(5,:),f.mnmax,f.ns-1)];
f.bsubumnc=[data1(8,:)' reshape(data(6,:),f.mnmax,f.ns-1)];
f.bsubvmnc=[data1(9,:)' reshape(data(7,:),f.mnmax,f.ns-1)];
f.bsubsmns=[data1(10,:)' reshape(data(8,:),f.mnmax,f.ns-1)];
f.bsupumnc=[data1(11,:)' reshape(data(9,:),f.mnmax,f.ns-1)];
f.bsupvmnc=[data1(12,:)' reshape(data(10,:),f.mnmax,f.ns-1)];
f.currvmnc=[data1(13,:)' reshape(data(11,:),f.mnmax,f.ns-1)];
% Read the half-mesh quantities
data=fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%g%g%g',[13 f.ns/2]);
f.iotas=data(1,:);
f.mass=data(2,:);
f.pres=data(3,:);
f.beta_vol(4,:);
f.phip=data(5,:);
f.buco=data(6,:);
f.bvco=data(7,:);
f.phi=data(8,:);
f.vp=data(9,:);
f.overr=data(10,:);
f.jcuru=data(11,:);
f.jcurv=data(12,:);
f.specw=data(13,:);
data=fscanf(fid,'%g%g%g%g%g%g',6);
f.aspect=data(1);
f.betatot=data(2);
f.betapol=data(3);
f.betator=data(4);
f.betaxis=data(5);
f.b0=data(6);
f.isigna=fscanf(fid,'%d\n',1);
f.input_extension=strtrim(fgetl(fid));
data=fscanf(fid,'%g',8);
f.IonLarmor=data(1);
f.VolAvgB=data(2);
f.RBtor0=data(3);
f.RBtor=data(4);
f.Itor=data(5);
f.Aminor=data(6);
f.Rmajor=data(7);
f.Volume=data(8);
% Mercier Criterion
data=fscanf(fid,'%g%g%g%g%g%g',[6 f.ns-2]);
f.Dmerc=data(1,:);
f.Dshear=data(2,:);
f.Dwell=data(3,:);
f.Dcurr=data(4,:);
f.Dgeod=data(5,:);
f.equif=data(6,:);
if (f.nextcur > 0)
    f.extcur=fscanf(fid,'%g',f.nextcur);
    f.curlabel=fscanf(fid,'%g',f.nextcur);
end
data=fscanf(fid,'%g%g',[2 f.nstore_seq]);
f.sqt=data(1,:);
f.wdot=data(1,:);
data=fscanf(fid,'%g%g',[2 f.nstore_seq]);
f.jdotb=data(1,:);
f.bdotgradv=data(2,:);
% Data and MSE Fits
if (f.ireconstruct > 0)
    if ((f.imse >= 2) || (f.itse >0))
        f.twsgt=fscanf(fid,'%g',1);
        f.msewgt=fscanf(fid,'%g',1);
        f.isnodes=fscanf(fid,'%d',1);
        data=fscanf(fid,'%g%g%g',[3 f.isnodes]);
        f.sknots=data(1,:);
        f.ystark=data(2,:);
        f.y2stark=data(3,:);
        f.ipnodes=fscanf(fid,'%d',1);
        data=fscanf(fid,'%g%g%g',[3 f.ipnodes]);
        f.pknots=data(1,:);
        f.ythom=data(2,:);
        f.y2thom=data(3,:);
        data=fscanf(fid,'%g%g%g%g%g%g%g',[7 (2*f.ns)-1]);
        f.anglemse=data(1,:);
        f.rmid=data(2,:);
        f.qmid=data(3,:);
        f.shear=data(4,:);
        f.presmid=data(5,:);
        f.alfa=data(6,:);
        f.curmid=data(7,:);
        data=fscanf(fid,'%g%g%g',[3 f.imse]);
        f.rstark=data(1,:);
        f.datastark=data(2,:);
        f.qmeas=data(3,:);
        data=fscanf(fid,'%g%g',[2 f.itse]);
        f.rthom=data(1,:);
        f.datathom=data(2,:);
    end
    if (f.nobd > 0)
        data=fscanf(fid,'%g%g%g',[3,f.nobd]);
        f.dsiext=data(1,:);
        f.plflux=data(2,:);
        f.dsiobt=data(3,:);
        f.flmwgt=fscanf(fid,'%g',1);
    end
    nbfldn=sum(nbldf(1:f.nbsets));
    if (nbfldn > 0)
        for n=1:nbsets
            data=fscanf(fid,'%g%g%g',[3 f.nbfld(n)]);
            f.bcoil(:,n)=data(1,:);
            f.plbfld(:,n)=data(2,:);
            f.bbc(:,n)=data(3,:);
        end
        f.bcwgt=fscanf(fid,'%g',1);
    end
    f.phidiam=fscanf(fid,'%g',1);
    f.delphid=fscanf(fid,'%g',1);
    % Read Limiter and Prout Plotting Specs
    f.nsets=fscanf(fid,'%g',1);
    f.nparts=fscanf(fid,'%g',1);
    f.nlim=fscanf(fid,'%g',1);
    f.nsetsn=fscanf(fid,'%g',f.nsets);
    f.pfcspec=zeros(f.nparts,max(f.nxsetsn),f.nsets);
    for k=1:f.nsets
        for j=1:f.nsetsn(k)
            for i=1:f.nparts
                f.pfcspec(i,j,k)=fscanf(fid,'%g',1);
            end
        end
    end
    f.limitr=fscanf(fid,'%g',f.nlim);
    f.rlim=zeros(max(f.limitr),f.nlim);
    f.zlim=zeros(max(f.limitr),f.nlim);
    for j=1:f.nlim
        for i=1:f.limitr(j)
            data=fscanf(fid,'%g%g',2);
            f.rlim(i,j)=data(1);
            f.zlim(i,j)=data(2);
        end
    end
    f.nrgrid=fscanf(fid,'%g',1);
    f.nzgrid=fscanf(fid,'%g',1);
    f.tokid=fscanf(fid,'%g',1);
    f.rx1=fscanf(fid,'%g',1);
    f.rx2=fscanf(fid,'%g',1);
    f.zy1=fscanf(fid,'%g',1);
    f.zy2=fscanf(fid,'%g',1);
    f.conif=fscanf(fid,'%g',1);
    f.imatch_phiedge=fscanf(fid,'%g',1);
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=read_vmec_650(fid,fmt)
% For 6.50

% Unit Conversions
dmu0=(2.0e-7)/pi;
% Read Data
data=fscanf(fid,'%g',6);
f.wb=data(1);
f.wp=data(2);
f.gamma=data(3);
f.pfac=data(4);
f.rmax_surf=data(5);
f.rmin_surf=data(6);
data=fscanf(fid,'%d',10);
f.nfp=data(1);
f.ns=data(2);
f.mpol=data(3);
f.ntor=data(4);
f.mnmax=data(5);
f.itfsq=data(6);
f.niter=data(7);
f.iasym=data(8);
f.ireconstruct=data(9);
f.ierr_vmec=data(10);
data=fscanf(fid,'%d',6);
f.imse=data(1);
f.itse=data(2);
f.nbsets=data(3);
f.nobd=data(4);
f.nextcur=data(5);
f.nstore_seq=data(6);
% Error Check
if (f.ierr_vmec && (f.ierr_vmec ~= 4))
    disp('ierr_vmec >0');
    return
end
% Read nbfld
if (f.nbsets > 0)
    f.nbfld=fscanf(fid,'%g',f.nbsets);
end
% Read mgrid filename
f.mgrid_file=fscanf(fid,'%s',1);
% Read Arrays
if f.iasym > 0
    data1=fscanf(fid,'%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g',[16 f.mnmax]);
    data=fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%g%g%g%g',[14 f.mnmax*(f.ns-1)]);
else
    data1=fscanf(fid,'%d%d%g%g%g%g%g%g%g%g%g%g%g',[13 f.mnmax]);
    data=fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%g',[11 f.mnmax*(f.ns-1)]);
end
% Extract Data from Arrays
f.xm=data1(1,:);
f.xn=data1(2,:);
% First reshape data
f.rmnc=[data1(3,:)' reshape(data(1,:),f.mnmax,f.ns-1)];
f.zmns=[data1(4,:)' reshape(data(2,:),f.mnmax,f.ns-1)];
f.lmns=[data1(5,:)' reshape(data(3,:),f.mnmax,f.ns-1)];
f.bmnc=[data1(6,:)' reshape(data(4,:),f.mnmax,f.ns-1)];
f.gmnc=[data1(7,:)' reshape(data(5,:),f.mnmax,f.ns-1)];
f.bsubumnc=[data1(8,:)' reshape(data(6,:),f.mnmax,f.ns-1)];
f.bsubvmnc=[data1(9,:)' reshape(data(7,:),f.mnmax,f.ns-1)];
f.bsubsmns=[data1(10,:)' reshape(data(8,:),f.mnmax,f.ns-1)];
f.bsupumnc=[data1(11,:)' reshape(data(9,:),f.mnmax,f.ns-1)];
f.bsupvmnc=[data1(12,:)' reshape(data(10,:),f.mnmax,f.ns-1)];
f.currvmnc=[data1(13,:)' reshape(data(11,:),f.mnmax,f.ns-1)];
% Read the half-mesh quantities
data=fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%g%g%g',[13 f.ns/2]);
f.iotas=data(1,:);
f.mass=data(2,:);
f.pres=data(3,:);
f.beta_vol(4,:);
f.phip=data(5,:);
f.buco=data(6,:);
f.bvco=data(7,:);
f.phi=data(8,:);
f.vp=data(9,:);
f.overr=data(10,:);
f.jcuru=data(11,:);
f.jcurv=data(12,:);
f.specw=data(13,:);
data=fscanf(fid,'%g%g%g%g%g%g',6);
f.aspect=data(1);
f.betatot=data(2);
f.betapol=data(3);
f.betator=data(4);
f.betaxis=data(5);
f.b0=data(6);
f.isigna=fscanf(fid,'%d\n',1);
f.input_extension=strtrim(fgetl(fid));
data=fscanf(fid,'%g',8);
f.IonLarmor=data(1);
f.VolAvgB=data(2);
f.RBtor0=data(3);
f.RBtor=data(4);
f.Itor=data(5);
f.Aminor=data(6);
f.Rmajor=data(7);
f.Volume=data(8);
% Mercier Criterion
data=fscanf(fid,'%g%g%g%g%g%g',[6 f.ns-2]);
f.Dmerc=data(1,:);
f.Dshear=data(2,:);
f.Dwell=data(3,:);
f.Dcurr=data(4,:);
f.Dgeod=data(5,:);
f.equif=data(6,:);
if (f.nextcur > 0)
    f.extcur=fscanf(fid,'%g',f.nextcur);
    f.curlabel=fscanf(fid,'%g',f.nextcur);
end
data=fscanf(fid,'%g%g',[2 f.nstore_seq]);
f.sqt=data(1,:);
f.wdot=data(1,:);
data=fscanf(fid,'%g%g',[2 f.nstore_seq]);
f.jdotb=data(1,:);
f.bdotgradv=data(2,:);
% Data and MSE Fits
if (f.ireconstruct > 0)
    if ((f.imse >= 2) || (f.itse >0))
        f.twsgt=fscanf(fid,'%g',1);
        f.msewgt=fscanf(fid,'%g',1);
        f.isnodes=fscanf(fid,'%d',1);
        data=fscanf(fid,'%g%g%g',[3 f.isnodes]);
        f.sknots=data(1,:);
        f.ystark=data(2,:);
        f.y2stark=data(3,:);
        f.ipnodes=fscanf(fid,'%d',1);
        data=fscanf(fid,'%g%g%g',[3 f.ipnodes]);
        f.pknots=data(1,:);
        f.ythom=data(2,:);
        f.y2thom=data(3,:);
        data=fscanf(fid,'%g%g%g%g%g%g%g',[7 (2*f.ns)-1]);
        f.anglemse=data(1,:);
        f.rmid=data(2,:);
        f.qmid=data(3,:);
        f.shear=data(4,:);
        f.presmid=data(5,:);
        f.alfa=data(6,:);
        f.curmid=data(7,:);
        data=fscanf(fid,'%g%g%g',[3 f.imse]);
        f.rstark=data(1,:);
        f.datastark=data(2,:);
        f.qmeas=data(3,:);
        data=fscanf(fid,'%g%g',[2 f.itse]);
        f.rthom=data(1,:);
        f.datathom=data(2,:);
    end
    if (f.nobd > 0)
        data=fscanf(fid,'%g%g%g',[3,f.nobd]);
        f.dsiext=data(1,:);
        f.plflux=data(2,:);
        f.dsiobt=data(3,:);
        f.flmwgt=fscanf(fid,'%g',1);
    end
    nbfldn=sum(nbldf(1:f.nbsets));
    if (nbfldn > 0)
        for n=1:nbsets
            data=fscanf(fid,'%g%g%g',[3 f.nbfld(n)]);
            f.bcoil(:,n)=data(1,:);
            f.plbfld(:,n)=data(2,:);
            f.bbc(:,n)=data(3,:);
        end
        f.bcwgt=fscanf(fid,'%g',1);
    end
    f.phidiam=fscanf(fid,'%g',1);
    f.delphid=fscanf(fid,'%g',1);
    % Read Limiter and Prout Plotting Specs
    f.nsets=fscanf(fid,'%g',1);
    f.nparts=fscanf(fid,'%g',1);
    f.nlim=fscanf(fid,'%g',1);
    f.nsetsn=fscanf(fid,'%g',f.nsets);
    f.pfcspec=zeros(f.nparts,max(f.nxsetsn),f.nsets);
    for k=1:f.nsets
        for j=1:f.nsetsn(k)
            for i=1:f.nparts
                f.pfcspec(i,j,k)=fscanf(fid,'%g',1);
            end
        end
    end
    f.limitr=fscanf(fid,'%g',f.nlim);
    f.rlim=zeros(max(f.limitr),f.nlim);
    f.zlim=zeros(max(f.limitr),f.nlim);
    for j=1:f.nlim
        for i=1:f.limitr(j)
            data=fscanf(fid,'%g%g',2);
            f.rlim(i,j)=data(1);
            f.zlim(i,j)=data(2);
        end
    end
    f.nrgrid=fscanf(fid,'%g',1);
    f.nzgrid=fscanf(fid,'%g',1);
    f.tokid=fscanf(fid,'%g',1);
    f.rx1=fscanf(fid,'%g',1);
    f.rx2=fscanf(fid,'%g',1);
    f.zy1=fscanf(fid,'%g',1);
    f.zy2=fscanf(fid,'%g',1);
    f.conif=fscanf(fid,'%g',1);
    f.imatch_phiedge=fscanf(fid,'%g',1);
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=read_vmec_695(fid,fmt)
f.lfreeb=0;
data=fscanf(fid,fmt,7);
f.wb=data(1);
f.wp=data(2);
f.gamma=data(3);
f.pfac=data(4);
f.rmax_surf=data(5);
f.rmin_surf=data(6);
f.zmax_surf=data(7);
data=fscanf(fid,fmt,10);
f.nfp=data(1);
f.ns=data(2);
f.mpol=data(3);
f.ntor=data(4);
f.mnmax=data(5);
f.itfsq=data(6);
f.niter=data(7);
f.iasym=data(8);
f.ireconstruct=data(9);
f.ierr_vmec=data(10);
data=fscanf(fid,fmt,6);
f.imse=data(1);
f.itse=data(2);
f.nbsets=data(3);
f.nobd=data(4);
f.nextcur=data(5);
f.nstore_seq=data(6);
% Error Check
if (f.ierr_vmec && (f.ierr_vmec ~= 4))
    disp(strcat('ierr_vmec:',num2str(f.ierr_vmec)));
    return
end
% Read nbfld
if (f.nbsets > 0)
    f.nbfld=fscanf(fid,fmt,f.nbsets);
end
% Read mgrid filename and setup other format statements
if isempty(strfind(fmt,','))
%   f.mgrid_file=fscanf(fid,'%s',1);
    fmt2='%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g';
    fmt3='%g%g%g%g%g%g%g%g%g%g%g%g%g%g';
    fmt4='%d%d%g%g%g%g%g%g%g%g%g%g%g';
    fmt5='%g%g%g%g%g%g%g%g%g%g%g';
    fmt6='%g%g%g%g%g%g%g%g%g%g%g%g%g';
    fmt7='%g%g%g%g%g%g';
    fmt8='%g%g';
    fmt9='%g%g%g';
    fmt10='%g%g%g%g%g%g%g';
else
%    f.mgrid_file=fscanf(fid,'%s,',1);
    fmt2='%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt3='%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt4='%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt5='%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt6='%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt7='%g,%g,%g,%g,%g,%g';
    fmt8='%g,%g';
    fmt9='%g,%g,%g';
    fmt10='%g,%g,%g,%g,%g,%g,%g';
end
if isempty(strfind(fmt,','))
    f.mgrid_file=fscanf(fid,'%s',1);
    fmt2='%g%g';
    fmt3='%g%g%g';
    fmt6='%g%g%g%g%g%g';
    fmt7='%g%g%g%g%g%g%g';
    fmt11='%g%g%g%g%g%g%g%g%g%g%g';
    fmt12='%g%g%g%g%g%g%g%g%g%g%g%g';
    fmt13='%g%g%g%g%g%g%g%g%g%g%g%g%g';
    fmt14='%g%g%g%g%g%g%g%g%g%g%g%g%g%g';
    fmt2_11='%d%d%g%g%g%g%g%g%g%g%g%g%g';
    fmt2_14='%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g';
else
    f.mgrid_file=fscanf(fid,'%s,',1);
    fmt2='%g,%g';
    fmt3='%g,%g,%g';
    fmt6='%g,%g,%g,%g,%g,%g';
    fmt7='%g,%g,%g,%g,%g,%g,%g';
    fmt11='%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt12='%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt13='%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt14='%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt2_11='%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt2_14='%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
end
% Read Arrays
if f.iasym > 0
    data1=fscanf(fid,fmt2_14,[16 f.mnmax]);
    data=fscanf(fid,fmt14,[14 f.mnmax*(f.ns-1)]);
else
    data1=fscanf(fid,fmt2_11,[13 f.mnmax]);
    data=fscanf(fid,fmt11,[11 f.mnmax*(f.ns-1)]);
end
% Extract Data from Arrays
f.xm=data1(1,:);
f.xn=data1(2,:);
% First reshape data
f.rmnc=[data1(3,:)' reshape(data(1,:),f.mnmax,f.ns-1)];
f.zmns=[data1(4,:)' reshape(data(2,:),f.mnmax,f.ns-1)];
f.lmns=[data1(5,:)' reshape(data(3,:),f.mnmax,f.ns-1)]; % On half grid
f.bmnc=[data1(6,:)' reshape(data(4,:),f.mnmax,f.ns-1)];
f.gmnc=[data1(7,:)' reshape(data(5,:),f.mnmax,f.ns-1)];
f.bsubumnc=[data1(8,:)' reshape(data(6,:),f.mnmax,f.ns-1)];
f.bsubvmnc=[data1(9,:)' reshape(data(7,:),f.mnmax,f.ns-1)];
f.bsubsmns=[data1(10,:)' reshape(data(8,:),f.mnmax,f.ns-1)];
f.bsupumnc=[data1(11,:)' reshape(data(9,:),f.mnmax,f.ns-1)];
f.bsupvmnc=[data1(12,:)' reshape(data(10,:),f.mnmax,f.ns-1)];
f.currvmnc=[data1(13,:)' reshape(data(11,:),f.mnmax,f.ns-1)];
if f.iasym >0
    f.rmns=[data1(14,:)' reshape(data(12,:),f.mnmax,f.ns-1)];
    f.zmnc=[data1(15,:)' reshape(data(13,:),f.mnmax,f.ns-1)];
    f.lmnc=[data1(16,:)' reshape(data(14,:),f.mnmax,f.ns-1)]; % On half grid
end
% Read the half-mesh quantities
data=fscanf(fid,fmt13,[13 f.ns-1]);
f.iotas=data(1,:);
f.mass=data(2,:);
f.pres=data(3,:);
f.beta_vol=data(4,:);
f.phip=data(5,:);
f.buco=data(6,:);
f.bvco=data(7,:);
f.phi=data(8,:);
f.vp=data(9,:);
f.overr=data(10,:);
f.jcuru=data(11,:);
f.jcurv=data(12,:);
f.specw=data(13,:);
data=fscanf(fid,fmt,6);
f.aspect=data(1);
f.betatot=data(2);
f.betapol=data(3);
f.betator=data(4);
f.betaxis=data(5);
f.b0=data(6);
f.isigna=fscanf(fid,strcat(fmt,'\n'),1);
f.input_extension=strtrim(fgetl(fid));
data=fscanf(fid,fmt,8);
f.IonLarmor=data(1);
f.VolAvgB=data(2);
f.RBtor0=data(3);
f.RBtor=data(4);
f.Itor=data(5);
f.Aminor=data(6);
f.Rmajor=data(7);
f.Volume=data(8);
% Mercier Criterion
data=fscanf(fid,fmt6,[6 f.ns-2]);
f.Dmerc=data(1,:);
f.Dshear=data(2,:);
f.Dwell=data(3,:);
f.Dcurr=data(4,:);
f.Dgeod=data(5,:);
f.equif=data(6,:);
f.curlabel=cell(f.nextcur,1);
if (f.nextcur > 0)
    f.lfreeb=1;
    f.extcur=fscanf(fid,fmt,f.nextcur);
    fscanf(fid,'\n');
    rem=f.nextcur;
    j=0;
    while rem > 0
        line=fgetl(fid);
        fscanf(fid,'\n');
        test=line(1);
        index=findstr(line,test);
        for i=1:size(index,2)/2;
            f.curlabel{i+j}=strtrim(line(index(2*i-1)+1:index(2*i)-1));
        end
        j=j+size(index,2)/2;
        rem=rem-size(index,2)/2;
    end
end
data=fscanf(fid,fmt2,[2 f.nstore_seq]);
f.sqt=data(1,:);
f.wdot=data(2,:);
data=fscanf(fid,fmt2,[2 f.nstore_seq]);
f.jdotb=data(1,:);
f.bdotgradv=data(2,:);
% No Unit Conversion Necessary
% Data and MSE Fits
if (f.ireconstruct > 0)
    if ((f.imse >= 2) || (f.itse >0))
        f.twsgt=fscanf(fid,fmt,1);
        f.msewgt=fscanf(fid,fmt,1);
        f.isnodes=fscanf(fid,fmt,1);
        data=fscanf(fid,fmt3,[3 f.isnodes]);
        f.sknots=data(1,:);
        f.ystark=data(2,:);
        f.y2stark=data(3,:);
        f.ipnodes=fscanf(fid,fmt,1);
        data=fscanf(fid,fmt3,[3 f.ipnodes]);
        f.pknots=data(1,:);
        f.ythom=data(2,:);
        f.y2thom=data(3,:);
        data=fscanf(fid,fmt7,[7 (2*f.ns)-1]);
        f.anglemse=data(1,:);
        f.rmid=data(2,:);
        f.qmid=data(3,:);
        f.shear=data(4,:);
        f.presmid=data(5,:);
        f.alfa=data(6,:);
        f.curmid=data(7,:);
        data=fscanf(fid,fmt3,[3 f.imse]);
        f.rstark=data(1,:);
        f.datastark=data(2,:);
        f.qmeas=data(3,:);
        data=fscanf(fid,fmt2,[2 f.itse]);
        f.rthom=data(1,:);
        f.datathom=data(2,:);
    end
    if (f.nobd > 0)
        data=fscanf(fid,fmt3,[3,f.nobd]);
        f.dsiext=data(1,:);
        f.plflux=data(2,:);
        f.dsiobt=data(3,:);
        f.flmwgt=fscanf(fid,fmt,1);
    end
    nbfldn=sum(nbldf(1:f.nbsets));
    if (nbfldn > 0)
        for n=1:nbsets
            data=fscanf(fid,fmt3,[3 f.nbfld(n)]);
            f.bcoil(:,n)=data(1,:);
            f.plbfld(:,n)=data(2,:);
            f.bbc(:,n)=data(3,:);
        end
        f.bcwgt=fscanf(fid,fmt,1);
    end
    f.phidiam=fscanf(fid,fmt,1);
    f.delphid=fscanf(fid,fmt,1);
    % Read Limiter and Prout Plotting Specs
    f.nsets=fscanf(fid,fmt,1);
    f.nparts=fscanf(fid,fmt,1);
    f.nlim=fscanf(fid,fmt,1);
    f.nsetsn=fscanf(fid,fmt,f.nsets);
    f.pfcspec=zeros(f.nparts,max(f.nxsetsn),f.nsets);
    for k=1:f.nsets
        for j=1:f.nsetsn(k)
            for i=1:f.nparts
                f.pfcspec(i,j,k)=fscanf(fid,fmt,1);
            end
        end
    end
    f.limitr=fscanf(fid,fmt,f.nlim);
    f.rlim=zeros(max(f.limitr),f.nlim);
    f.zlim=zeros(max(f.limitr),f.nlim);
    for j=1:f.nlim
        for i=1:f.limitr(j)
            data=fscanf(fid,fmt2,2);
            f.rlim(i,j)=data(1);
            f.zlim(i,j)=data(2);
        end
    end
    f.nrgrid=fscanf(fid,fmt,1);
    f.nzgrid=fscanf(fid,fmt,1);
    f.tokid=fscanf(fid,fmt,1);
    f.rx1=fscanf(fid,fmt,1);
    f.rx2=fscanf(fid,fmt,1);
    f.zy1=fscanf(fid,fmt,1);
    f.zy2=fscanf(fid,fmt,1);
    f.conif=fscanf(fid,fmt,1);
    f.imatch_phiedge=fscanf(fid,fmt,1);
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=read_vmec_800(fid,fmt)
f.lfreeb=0;
data=fscanf(fid,fmt,7);
f.wb=data(1);
f.wp=data(2);
f.gamma=data(3);
f.pfac=data(4);
f.rmax_surf=data(5);
f.rmin_surf=data(6);
f.zmax_surf=data(7);
data=fscanf(fid,fmt,10);
f.nfp=data(1);
f.ns=data(2);
f.mpol=data(3);
f.ntor=data(4);
f.mnmax=data(5);
f.itfsq=data(6);
f.niter=data(7);
f.iasym=data(8);
f.ireconstruct=data(9);
f.ierr_vmec=data(10);
f.mnmax_nyq=f.mnmax;
data=fscanf(fid,fmt,6);
f.imse=data(1);
f.itse=data(2);
f.nbsets=data(3);
f.nobd=data(4);
f.nextcur=data(5);
f.nstore_seq=data(6);
% Error Check
if (f.ierr_vmec && (f.ierr_vmec ~= 4))
    disp(strcat('ierr_vmec:',num2str(f.ierr_vmec)));
    return
end
% Read nbfld
if (f.nbsets > 0)
    f.nbfld=fscanf(fid,fmt,f.nbsets);
end
% Read mgrid filename and setup other format statements
if isempty(strfind(fmt,','))
    f.mgrid_file=fscanf(fid,'%s',1);
    fmt2='%g%g';
    fmt3='%g%g%g';
    fmt6='%g%g%g%g%g%g';
    fmt7='%g%g%g%g%g%g%g';
    fmt10='%g%g%g%g%g%g%g%g%g%g';
    fmt11='%g%g%g%g%g%g%g%g%g%g%g';
    fmt12='%g%g%g%g%g%g%g%g%g%g%g%g';
    fmt13='%g%g%g%g%g%g%g%g%g%g%g%g%g';
    fmt14='%g%g%g%g%g%g%g%g%g%g%g%g%g%g';
    fmt2_11='%d%d%g%g%g%g%g%g%g%g%g%g%g';
    fmt2_14='%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g';
else
    f.mgrid_file=fscanf(fid,'%s,',1);
    fmt2='%g,%g';
    fmt3='%g,%g,%g';
    fmt6='%g,%g,%g,%g,%g,%g';
    fmt7='%g,%g,%g,%g,%g,%g,%g';
    fmt11='%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt12='%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt13='%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt14='%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt2_11='%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt2_14='%d,%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
end
% Read Arrays
if f.iasym > 0
    data1=fscanf(fid,fmt2_14,[16 f.mnmax]);
    data=fscanf(fid,fmt14,[14 f.mnmax*(f.ns-1)]);
else
    data1=fscanf(fid,fmt2_11,[13 f.mnmax]);
    data=fscanf(fid,fmt11,[11 f.mnmax*(f.ns-1)]);
end
% Extract Data from Arrays
f.xm=data1(1,:);
f.xn=data1(2,:);
% First reshape data
f.rmnc=[data1(3,:)' reshape(data(1,:),f.mnmax,f.ns-1)];
f.zmns=[data1(4,:)' reshape(data(2,:),f.mnmax,f.ns-1)];
f.lmns=[data1(5,:)' reshape(data(3,:),f.mnmax,f.ns-1)]; % On half grid
f.bmnc=[data1(6,:)' reshape(data(4,:),f.mnmax,f.ns-1)];
f.gmnc=[data1(7,:)' reshape(data(5,:),f.mnmax,f.ns-1)];
f.bsubumnc=[data1(8,:)' reshape(data(6,:),f.mnmax,f.ns-1)];
f.bsubvmnc=[data1(9,:)' reshape(data(7,:),f.mnmax,f.ns-1)];
f.bsubsmns=[data1(10,:)' reshape(data(8,:),f.mnmax,f.ns-1)];
f.bsupumnc=[data1(11,:)' reshape(data(9,:),f.mnmax,f.ns-1)];
f.bsupvmnc=[data1(12,:)' reshape(data(10,:),f.mnmax,f.ns-1)];
f.currvmnc=[data1(13,:)' reshape(data(11,:),f.mnmax,f.ns-1)];
if f.iasym >0
    f.rmns=[data1(14,:)' reshape(data(12,:),f.mnmax,f.ns-1)];
    f.zmnc=[data1(15,:)' reshape(data(13,:),f.mnmax,f.ns-1)];
    f.lmnc=[data1(16,:)' reshape(data(14,:),f.mnmax,f.ns-1)]; % On half grid
end
% Read the full-mesh quantities
data=fscanf(fid,fmt6,[6 f.ns]);
f.iotaf=data(1,:);
f.presf=data(2,:);
f.phipf=data(3,:);
f.phi=data(4,:);
f.jcuru=data(5,:);
f.jcurv=data(6,:);
% Read the half-mesh quantities
data=fscanf(fid,fmt10,[10 f.ns-1]);
f.iotas=data(1,:);
f.mass=data(2,:);
f.pres=data(3,:);
f.beta_vol=data(4,:);
f.phip=data(5,:);
f.buco=data(6,:);
f.bvco=data(7,:);
f.vp=data(8,:);
f.overr=data(9,:);
f.specw=data(10,:);
data=fscanf(fid,fmt,6);
f.aspect=data(1);
f.betatot=data(2);
f.betapol=data(3);
f.betator=data(4);
f.betaxis=data(5);
f.b0=data(6);
f.isigna=fscanf(fid,strcat(fmt,'\n'),1);
f.input_extension=strtrim(fgetl(fid));
data=fscanf(fid,fmt,8);
f.IonLarmor=data(1);
f.VolAvgB=data(2);
f.RBtor0=data(3);
f.RBtor=data(4);
f.Itor=data(5);
f.Aminor=data(6);
f.Rmajor=data(7);
f.Volume=data(8);
% Mercier Criterion
data=fscanf(fid,fmt6,[6 f.ns-2]);
f.Dmerc=data(1,:);
f.Dshear=data(2,:);
f.Dwell=data(3,:);
f.Dcurr=data(4,:);
f.Dgeod=data(5,:);
f.equif=data(6,:);
f.curlabel=cell(f.nextcur,1);
if (f.nextcur > 0)
    f.lfreeb=1;
    f.extcur=fscanf(fid,fmt,f.nextcur);
    fscanf(fid,'\n');
    rem=f.nextcur;
    j=0;
    while rem > 0
        line=fgetl(fid);
        fscanf(fid,'\n');
        test=line(1);
        index=findstr(line,test);
        for i=1:size(index,2)/2;
            f.curlabel{i+j}=strtrim(line(index(2*i-1)+1:index(2*i)-1));
        end
        j=j+size(index,2)/2;
        rem=rem-size(index,2)/2;
    end
end
data=fscanf(fid,fmt2,[2 f.nstore_seq]);
f.sqt=data(1,:);
f.wdot=data(2,:);
data=fscanf(fid,fmt2,[2 f.nstore_seq]);
f.jdotb=data(1,:);
f.bdotgradv=data(2,:);
% No Unit Conversion Necessary
% Data and MSE Fits
if (f.ireconstruct > 0)
    if ((f.imse >= 2) || (f.itse >0))
        f.twsgt=fscanf(fid,fmt,1);
        f.msewgt=fscanf(fid,fmt,1);
        f.isnodes=fscanf(fid,fmt,1);
        data=fscanf(fid,fmt3,[3 f.isnodes]);
        f.sknots=data(1,:);
        f.ystark=data(2,:);
        f.y2stark=data(3,:);
        f.ipnodes=fscanf(fid,fmt,1);
        data=fscanf(fid,fmt3,[3 f.ipnodes]);
        f.pknots=data(1,:);
        f.ythom=data(2,:);
        f.y2thom=data(3,:);
        data=fscanf(fid,fmt7,[7 (2*f.ns)-1]);
        f.anglemse=data(1,:);
        f.rmid=data(2,:);
        f.qmid=data(3,:);
        f.shear=data(4,:);
        f.presmid=data(5,:);
        f.alfa=data(6,:);
        f.curmid=data(7,:);
        data=fscanf(fid,fmt3,[3 f.imse]);
        f.rstark=data(1,:);
        f.datastark=data(2,:);
        f.qmeas=data(3,:);
        data=fscanf(fid,fmt2,[2 f.itse]);
        f.rthom=data(1,:);
        f.datathom=data(2,:);
    end
    if (f.nobd > 0)
        data=fscanf(fid,fmt3,[3,f.nobd]);
        f.dsiext=data(1,:);
        f.plflux=data(2,:);
        f.dsiobt=data(3,:);
        f.flmwgt=fscanf(fid,fmt,1);
    end
    nbfldn=sum(nbldf(1:f.nbsets));
    if (nbfldn > 0)
        for n=1:nbsets
            data=fscanf(fid,fmt3,[3 f.nbfld(n)]);
            f.bcoil(:,n)=data(1,:);
            f.plbfld(:,n)=data(2,:);
            f.bbc(:,n)=data(3,:);
        end
        f.bcwgt=fscanf(fid,fmt,1);
    end
    f.phidiam=fscanf(fid,fmt,1);
    f.delphid=fscanf(fid,fmt,1);
    % Read Limiter and Prout Plotting Specs
    f.nsets=fscanf(fid,fmt,1);
    f.nparts=fscanf(fid,fmt,1);
    f.nlim=fscanf(fid,fmt,1);
    f.nsetsn=fscanf(fid,fmt,f.nsets);
    f.pfcspec=zeros(f.nparts,max(f.nxsetsn),f.nsets);
    for k=1:f.nsets
        for j=1:f.nsetsn(k)
            for i=1:f.nparts
                f.pfcspec(i,j,k)=fscanf(fid,fmt,1);
            end
        end
    end
    f.limitr=fscanf(fid,fmt,f.nlim);
    f.rlim=zeros(max(f.limitr),f.nlim);
    f.zlim=zeros(max(f.limitr),f.nlim);
    for j=1:f.nlim
        for i=1:f.limitr(j)
            data=fscanf(fid,fmt2,2);
            f.rlim(i,j)=data(1);
            f.zlim(i,j)=data(2);
        end
    end
    f.nrgrid=fscanf(fid,fmt,1);
    f.nzgrid=fscanf(fid,fmt,1);
    f.tokid=fscanf(fid,fmt,1);
    f.rx1=fscanf(fid,fmt,1);
    f.rx2=fscanf(fid,fmt,1);
    f.zy1=fscanf(fid,fmt,1);
    f.zy2=fscanf(fid,fmt,1);
    f.conif=fscanf(fid,fmt,1);
    f.imatch_phiedge=fscanf(fid,fmt,1);
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=read_vmec_847(fid,fmt)
mu0=4.*pi.*1e-7;
f.lfreeb=0;
data=fscanf(fid,fmt,7);
f.wb=data(1);
f.wp=data(2);
f.gamma=data(3);
f.pfac=data(4);
f.rmax_surf=data(5);
f.rmin_surf=data(6);
f.zmax_surf=data(7);
data=fscanf(fid,fmt,11);
f.nfp=data(1);
f.ns=data(2);
f.mpol=data(3);
f.ntor=data(4);
f.mnmax=data(5);
f.mnmax_nyq=data(6);
f.itfsq=data(7);
f.niter=data(8);
f.iasym=data(9);
f.ireconstruct=data(10);
f.ierr_vmec=data(11);
data=fscanf(fid,fmt,6);
f.imse=data(1);
f.itse=data(2);
f.nbsets=data(3);
f.nobd=data(4);
f.nextcur=data(5);
f.nstore_seq=data(6);
% Error Check
if (f.ierr_vmec && (f.ierr_vmec ~= 4))
    disp(strcat('ierr_vmec:',num2str(f.ierr_vmec)));
    return
end
% Read nbfld
if (f.nbsets > 0)
    f.nbfld=fscanf(fid,fmt,f.nbsets);
end
% Read mgrid filename and setup other format statements
if isempty(strfind(fmt,','))
    f.mgrid_file=fscanf(fid,'%s',1);
    fmt2='%g%g';
    fmt3='%g%g%g';
    fmt6='%g%g%g%g%g%g';
    fmt7='%g%g%g%g%g%g%g';
    fmt10='%g%g%g%g%g%g%g%g%g%g';
    fmt13='%g%g%g%g%g%g%g%g%g%g%g%g%g';
    fmt14='%g%g%g%g%g%g%g%g%g%g%g%g%g%g';
    fmt20='%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g%g';
    fmt2_3_2_7='%d%d%g%g%g%d%d%g%g%g%g%g%g%g';
    fmt2_6_2_14='%d%d%g%g%g%g%g%g%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g';
    fmt2_7='%d%d%g%g%g%g%g%g%g';
    fmt2_14='%d%d%g%g%g%g%g%g%g%g%g%g%g%g%g%g';
else
    f.mgrid_file=fscanf(fid,'%s,',1);
    fmt2='%g,%g';
    fmt3='%g,%g,%g';
    fmt6='%g,%g,%g,%g,%g,%g';
    fmt7='%g,%g,%g,%g,%g,%g,%g';
    fmt13='%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g';
    fmt2_3='%d,%d,%g,%g,%g';
    fmt2_6='%d,%d,%g,%g,%g,%g,%g,%g';
end
% Read Arrays
f.xm=zeros(1,f.mnmax);
f.xn=zeros(1,f.mnmax);
f.rmnc=zeros(f.mnmax,f.ns);
f.zmns=zeros(f.mnmax,f.ns);
f.lmns=zeros(f.mnmax,f.ns);
f.xm_nyq=zeros(1,f.mnmax_nyq);
f.xn_nyq=zeros(1,f.mnmax_nyq);
f.bmnc=zeros(f.mnmax_nyq,f.ns);
f.gmnc=zeros(f.mnmax_nyq,f.ns);
f.bsubumnc=zeros(f.mnmax_nyq,f.ns);
f.bsubvmnc=zeros(f.mnmax_nyq,f.ns);
f.bsubsmns=zeros(f.mnmax_nyq,f.ns);
f.bsupumnc=zeros(f.mnmax_nyq,f.ns);
f.bsupvmnc=zeros(f.mnmax_nyq,f.ns);
if f.iasym >0
    f.rmns=zeros(f.mnmax,f.ns);
    f.zmnc=zeros(f.mnmax,f.ns);
    f.lmnc=zeros(f.mnmax,f.ns);
    f.bmns=zeros(f.mnmax_nyq,f.ns);
    f.gmns=zeros(f.mnmax_nyq,f.ns);
    f.bsubumns=zeros(f.mnmax_nyq,f.ns);
    f.bsubvmns=zeros(f.mnmax_nyq,f.ns);
    f.bsubsmnc=zeros(f.mnmax_nyq,f.ns);
    f.bsupumns=zeros(f.mnmax_nyq,f.ns);
    f.bsupvmns=zeros(f.mnmax_nyq,f.ns);
end
for i=1:f.ns
    for j=1:f.mnmax
        if i==1
            f.xm(j)=fscanf(fid,'%d',1);
            f.xn(j)=fscanf(fid,'%d',1);
        end
        f.rmnc(j,i)=fscanf(fid,'%g',1);
        f.zmns(j,i)=fscanf(fid,'%g',1);
        f.lmns(j,i)=fscanf(fid,'%g',1);
        if f.iasym > 0
            f.rmns(j,i)=fscanf(fid,'%g',1);
            f.zmnc(j,i)=fscanf(fid,'%g',1);
            f.lmnc(j,i)=fscanf(fid,'%g',1);
        end
    end
    for j=1:f.mnmax_nyq
        if i==1
            f.xm_nyq(j)=fscanf(fid,'%d',1);
            f.xn_nyq(j)=fscanf(fid,'%d',1);
        end
        f.bmnc(j,i)=fscanf(fid,'%g',1);
        f.gmnc(j,i)=fscanf(fid,'%g',1);
        f.bsubumnc(j,i)=fscanf(fid,'%g',1);
        f.bsubvmnc(j,i)=fscanf(fid,'%g',1);
        f.bsubsmns(j,i)=fscanf(fid,'%g',1);
        f.bsupumnc(j,i)=fscanf(fid,'%g',1);
        f.bsupvmnc(j,i)=fscanf(fid,'%g',1);
        if f.iasym > 0
            f.bmns(j,i)=fscanf(fid,'%g',1);
            f.gmns(j,i)=fscanf(fid,'%g',1);
            f.bsubumns(j,i)=fscanf(fid,'%g',1);
            f.bsubvmns(j,i)=fscanf(fid,'%g',1);
            f.bsubsmnc(j,i)=fscanf(fid,'%g',1);
            f.bsupumns(j,i)=fscanf(fid,'%g',1);
            f.bsupvmns(j,i)=fscanf(fid,'%g',1);
        end
    end
end
f.mnyq=max(f.xm_nyq);
f.nnyq=max(f.xn_nyq)./f.nfp;
% Calculate the Currentsf.currumnc=zeros(f.mnmaxnyq,f.ns);
f.currvmnc=zeros(f.mnmax_nyq,f.ns);
ohs = f.ns-1;
hs  = 1.0/ohs;
ns = f.ns;
for i=2:ns
    shalf(i) = sqrt(hs*(i-1.5));
    sfull(i) = sqrt(hs*(i-1));
end
js1 = 3:ns;
js  = 2:(ns-1);
for mn = 1:f.mnmax_nyq
    if (mod(f.xm_nyq,2) == 1)
        t1  = 0.5.*(shalf(js1).*f.bsubsmns(mn,js1)+...
            shalf(js).*f.bsubsmns(mn,js))./sfull(js);
        bu0 = f.bsubumnc(mn,js)./shalf(js);
        bu1 = f.bsubumnc(mn,js1)./shalf(js1);
        t2  = ohs.*(bu1-bu0).*sfull(js)+0.25.*(bu0+bu1)./sfull(js);
        bv0 = f.bsubvmnc(mn,js)./shalf(js);
        bv1 = f.bsubvmnc(mn,js1)./shalf(js1);
        t3  = ohs.*(bv1-bv0).*sfull(js)+0.25.*(bv0+bv1)./sfull(js);
    else
        t1  = 0.5.*(f.bsubsmns(mn,js1)+f.bsubsmns(mn,js));
        t2  = ohs.*(f.bsubumnc(mn,js1)+f.bsubumnc(mn,js));
        t3  = ohs.*(f.bsubvmnc(mn,js1)+f.bsubvmnc(mn,js));
    end
    f.currumnc(mn,js) = -double(f.xn_nyq(mn)).*t1 - t3;
    f.currvmnc(mn,js) = -double(f.xn_nyq(mn)).*t1 + t2;
end
% OLD Way
%for i=2:f.ns-1
%    f.currumnc(:,i)=-double(f.xnnyq)'.*f.bsubsmns(:,i)-(f.ns-1).*(f.bsubvmnc(:,i+1)-f.bsubvmnc(:,i));
%    f.currvmnc(:,i)=-double(f.xmnyq)'.*f.bsubsmns(:,i)+(f.ns-1).*(f.bsubumnc(:,i+1)-f.bsubumnc(:,i));
%end
f.currumnc(:,1)=0.0;
f.currvmnc(:,1)=0.0;
for i=1:f.mnmax_nyq
    if (f.xm_nyq(i)==0)
        f.currumnc(i,1)=2.*f.currumnc(i,2)-f.currumnc(i,3);
        f.currvmnc(i,1)=2.*(f.ns-1).*f.bsubumnc(i,2);
    end
end
f.currumnc(:,f.ns)=2.*f.currumnc(:,f.ns-1)-f.currumnc(:,f.ns-2);
f.currvmnc(:,f.ns)=2.*f.currvmnc(:,f.ns-1)-f.currvmnc(:,f.ns-2);
f.currumnc=f.currumnc./mu0;
f.currvmnc=f.currvmnc./mu0;
if f.iasym
    f.currumns=zeros(f.mnmaxnyq,f.ns);
    f.currvmns=zeros(f.mnmaxnyq,f.ns);
    for mn = 1:f.mnmax_nyq
        if (mod(f.xm_nyq,2) == 1)
            t1  = 0.5.*(shalf(js1).*f.bsubsmnc(mn,js1)+...
                shalf(js).*f.bsubsmnc(mn,js))./sfull(js);
            bu0 = f.bsubumns(mn,js)./shalf(js);
            bu1 = f.bsubumns(mn,js1)./shalf(js1);
            t2  = ohs.*(bu1-bu0).*sfull(js)+0.25.*(bu0+bu1)./sfull(js);
            bv0 = f.bsubvmns(mn,js)./shalf(js);
            bv1 = f.bsubvmns(mn,js1)./shalf(js1);
            t3  = ohs.*(bv1-bv0).*sfull(js)+0.25.*(bv0+bv1)./sfull(js);
        else
            t1  = 0.5.*(f.bsubsmnc(mn,js1)+f.bsubsmnc(mn,js));
            t2  = ohs.*(f.bsubumns(mn,js1)+f.bsubumns(mn,js));
            t3  = ohs.*(f.bsubvmns(mn,js1)+f.bsubvmns(mn,js));
        end
        f.currumns(mn,js) = -double(f.xn_nyq(mn)).*t1 - t3;
        f.currvmns(mn,js) = -double(f.xn_nyq(mn)).*t1 + t2;
    end
    % OLD WAY
    %for i=2:f.ns-1
    %    f.currumns(:,i)= double(f.xnnyq)'.*f.bsubsmnc(:,i)-(f.ns-1).*(f.bsubvmns(:,i+1)-f.bsubvmns(:,i));
    %    f.currvmns(:,i)= double(f.xmnyq)'.*f.bsubsmnc(:,i)+(f.ns-1).*(f.bsubumns(:,i+1)-f.bsubumns(:,i));
    %end
    f.currumns(:,1)=0.0;
    f.currvmns(:,1)=0.0;
    for i=1:f.mnmax_nyq
        if (f.xm_nyq(i)==0)
            f.currumns(i,1)=2.*f.currumns(i,2)-f.currumns(i,3);
            f.currvmns(i,1)=2.*(f.ns-1).*f.bsubumns(i,2);
        end
    end
    f.currumns(:,f.ns)=2.*f.currumns(:,f.ns-1)-f.currumns(:,f.ns-2);
    f.currvmns(:,f.ns)=2.*f.currvmns(:,f.ns-1)-f.currvmns(:,f.ns-2);
    f.currumns=f.currumns./mu0;
    f.currvmns=f.currvmns./mu0;
end
% Read the full-mesh quantities
data=fscanf(fid,fmt6,[6 f.ns]);
f.iotaf=data(1,:);
f.presf=data(2,:);
f.phipf=data(3,:);
f.phi=data(4,:);
f.jcuru=data(5,:);
f.jcurv=data(6,:);
% Read the half-mesh quantities
data=fscanf(fid,fmt10,[10 f.ns-1]);
f.iotas=data(1,:);
f.mass=data(2,:);
f.pres=data(3,:);
f.beta_vol=data(4,:);
f.phip=data(5,:);
f.buco=data(6,:);
f.bvco=data(7,:);
f.vp=data(8,:);
f.overr=data(9,:);
f.specw=data(10,:);
data=fscanf(fid,fmt,6);
f.aspect=data(1);
f.betatot=data(2);
f.betapol=data(3);
f.betator=data(4);
f.betaxis=data(5);
f.b0=data(6);
f.isigna=fscanf(fid,strcat(fmt,'\n'),1);
f.input_extension=strtrim(fgetl(fid));
data=fscanf(fid,fmt,8);
f.IonLarmor=data(1);
f.VolAvgB=data(2);
f.RBtor0=data(3);
f.RBtor=data(4);
f.Itor=data(5);
f.Aminor=data(6);
f.Rmajor=data(7);
f.Volume=data(8);
% Mercier Criterion
data=fscanf(fid,fmt6,[6 f.ns-2]);
f.Dmerc=data(1,:);
f.Dshear=data(2,:);
f.Dwell=data(3,:);
f.Dcurr=data(4,:);
f.Dgeod=data(5,:);
f.equif=data(6,:);
f.curlabel=cell(f.nextcur,1);
if (f.nextcur > 0)
    f.lfreeb=1;
    f.extcur=fscanf(fid,fmt,f.nextcur);
    lcurr=strtrim(fscanf(fid,'%s',1));
    if strcmpi(lcurr,'T')
        fscanf(fid,'\n');
        rem=f.nextcur;
        j=0;
        while rem > 0
            line=fgetl(fid);
            fscanf(fid,'\n');
            test=line(1);
            index=findstr(line,test);
            for i=1:size(index,2)/2;
                f.curlabel{i+j}=strtrim(line(index(2*i-1)+1:index(2*i)-1));
            end
            j=j+size(index,2)/2;
            rem=rem-size(index,2)/2;
        end
    end
end
data=fscanf(fid,fmt2,[2 f.nstore_seq]);
f.sqt=data(1,:);
f.wdot=data(2,:);
data=fscanf(fid,fmt2,[2 f.nstore_seq]);
f.jdotb=data(1,:);
f.bdotgradv=data(2,:);
% No Unit Conversion Necessary
% Data and MSE Fits
if (f.ireconstruct > 0)
    if ((f.imse >= 2) || (f.itse >0))
        f.twsgt=fscanf(fid,fmt,1);
        f.msewgt=fscanf(fid,fmt,1);
        f.isnodes=fscanf(fid,fmt,1);
        data=fscanf(fid,fmt3,[3 f.isnodes]);
        f.sknots=data(1,:);
        f.ystark=data(2,:);
        f.y2stark=data(3,:);
        f.ipnodes=fscanf(fid,fmt,1);
        data=fscanf(fid,fmt3,[3 f.ipnodes]);
        f.pknots=data(1,:);
        f.ythom=data(2,:);
        f.y2thom=data(3,:);
        data=fscanf(fid,fmt7,[7 (2*f.ns)-1]);
        f.anglemse=data(1,:);
        f.rmid=data(2,:);
        f.qmid=data(3,:);
        f.shear=data(4,:);
        f.presmid=data(5,:);
        f.alfa=data(6,:);
        f.curmid=data(7,:);
        data=fscanf(fid,fmt3,[3 f.imse]);
        f.rstark=data(1,:);
        f.datastark=data(2,:);
        f.qmeas=data(3,:);
        data=fscanf(fid,fmt2,[2 f.itse]);
        f.rthom=data(1,:);
        f.datathom=data(2,:);
    end
    if (f.nobd > 0)
        data=fscanf(fid,fmt3,[3,f.nobd]);
        f.dsiext=data(1,:);
        f.plflux=data(2,:);
        f.dsiobt=data(3,:);
        f.flmwgt=fscanf(fid,fmt,1);
    end
    nbfldn=sum(nbldf(1:f.nbsets));
    if (nbfldn > 0)
        for n=1:nbsets
            data=fscanf(fid,fmt3,[3 f.nbfld(n)]);
            f.bcoil(:,n)=data(1,:);
            f.plbfld(:,n)=data(2,:);
            f.bbc(:,n)=data(3,:);
        end
        f.bcwgt=fscanf(fid,fmt,1);
    end
    f.phidiam=fscanf(fid,fmt,1);
    f.delphid=fscanf(fid,fmt,1);
    % Read Limiter and Prout Plotting Specs
    f.nsets=fscanf(fid,fmt,1);
    f.nparts=fscanf(fid,fmt,1);
    f.nlim=fscanf(fid,fmt,1);
    f.nsetsn=fscanf(fid,fmt,f.nsets);
    f.pfcspec=zeros(f.nparts,max(f.nxsetsn),f.nsets);
    for k=1:f.nsets
        for j=1:f.nsetsn(k)
            for i=1:f.nparts
                f.pfcspec(i,j,k)=fscanf(fid,fmt,1);
            end
        end
    end
    f.limitr=fscanf(fid,fmt,f.nlim);
    f.rlim=zeros(max(f.limitr),f.nlim);
    f.zlim=zeros(max(f.limitr),f.nlim);
    for j=1:f.nlim
        for i=1:f.limitr(j)
            data=fscanf(fid,fmt2,2);
            f.rlim(i,j)=data(1);
            f.zlim(i,j)=data(2);
        end
    end
    f.nrgrid=fscanf(fid,fmt,1);
    f.nzgrid=fscanf(fid,fmt,1);
    f.tokid=fscanf(fid,fmt,1);
    f.rx1=fscanf(fid,fmt,1);
    f.rx2=fscanf(fid,fmt,1);
    f.zy1=fscanf(fid,fmt,1);
    f.zy2=fscanf(fid,fmt,1);
    f.conif=fscanf(fid,fmt,1);
    f.imatch_phiedge=fscanf(fid,fmt,1);
end
f.mgrid_mode=strtrim(fscanf(fid,'%s'));
return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=half2fullmesh(f)
% Interpolate various quantities to the full mesh
% The following quantities are on 1D half mesh
% iotas, mass, pres, beta_vol, buco, bvco, vp, specw, phip, jdotb,
% bdotgradv
% The following quantities are on 1D full mesh
% iotaf, presf, phi, phipf, jcuru, jcurv
% The following quantities are on 1D full mesh from 2:ns-1
% Dmerc, Dshear, Dwell, Dcurr, Dgeod

% Now we create full mesh quantities for those that have full mesh names
if ~isfield(f,'iotaf')
    f.iotaf=h2f(f.iotas,f.ns);
end
if ~isfield(f,'presf')
    f.presf=h2f(f.pres,f.ns);
end
if ~isfield(f,'phipf')
    f.phipf=h2f(f.phip,f.ns);
end
% Now the rest of the quantities are remapped to full mesh
f.beta_vol=h2f(f.beta_vol,f.ns);
f.buco=h2f(f.buco,f.ns);
f.bvco=h2f(f.bvco,f.ns);
f.vp=h2f(f.vp,f.ns);
f.overr=h2f(f.overr,f.ns);
f.specw=h2f(f.specw,f.ns);
if(length(f.jdotb) == f.ns)
    f.jdotb=h2f(f.jdotb,f.ns);
    f.bdotgradv=h2f(f.bdotgradv,f.ns);
end
% Check to make sure phi,jcuru and jcurv are the same size as f.ns
if length(f.phi) ~= f.ns
    temp=zeros(1,f.ns);
    temp(2:f.ns)=f.phi;
    temp(1)=2.*temp(2)-temp(3);
    temp(f.ns)=2.*temp(f.ns-1)-temp(f.ns-2);
    f.phi=temp;
end
if length(f.jcuru) ~= f.ns
    temp=zeros(1,f.ns);
    temp(2:f.ns)=f.jcuru;
    temp(1)=2.*temp(2)-temp(3);
    temp(f.ns)=2.*temp(f.ns-1)-temp(f.ns-2);
    f.jcuru=temp;
end
if length(f.jcurv) ~= f.ns
    temp=zeros(1,f.ns);
    temp(2:f.ns)=f.jcurv;
    temp(1)=2.*temp(2)-temp(3);
    temp(f.ns)=2.*temp(f.ns-1)-temp(f.ns-2);
    f.jcurv=temp;
end
% Put the FLOW arrays on full mesh
if isfield(f,'pmap')
    temp = zeros(1,f.ns);
    temp(2:f.ns) = f.pmap(1:f.ns-1);
    temp(1)=2.*temp(2)-temp(3);
    temp(f.ns)=2.*temp(f.ns-1)-temp(f.ns-2);
    f.pmap = temp;
end
if isfield(f,'omega')  % Note this is zero on axis so extrapolate first
    temp = zeros(1,f.ns);
    temp(2:f.ns) = f.omega(1:f.ns-1);
    temp(2) = 2*temp(3)-temp(4);
    temp(1)=2.*temp(2)-temp(3);
    temp(f.ns)=2.*temp(f.ns-1)-temp(f.ns-2);
    f.omega = temp;
end
if isfield(f,'tpotb') % Note this is zero on axis so extrapolate first
    temp = zeros(1,f.ns);
    temp(2:f.ns) = f.tpotb(1:f.ns-1);
    temp(2) = 2*temp(3)-temp(4);
    temp(1)=2.*temp(2)-temp(3);
    temp(f.ns)=2.*temp(f.ns-1)-temp(f.ns-2);
    f.tpotb = temp;
end
% Fix the Stability values by making into ns arrays
if length(f.Dmerc) ~= f.ns
    temp=zeros(1,f.ns);
    temp(2:f.ns-1)=f.Dmerc;
    f.Dmerc=temp;
    f.Dmerc(1)=2 * f.Dmerc(2) - f.Dmerc(3);
    f.Dmerc(f.ns)=2 * f.Dmerc(f.ns-1) - f.Dmerc(f.ns-2);
    temp(2:f.ns-1)=f.Dshear;
    f.Dshear=temp;
    f.Dshear(1)=2 * f.Dshear(2) - f.Dshear(3);
    f.Dshear(f.ns)=2 * f.Dshear(f.ns-1) - f.Dshear(f.ns-2);
    temp(2:f.ns-1)=f.Dwell;
    f.Dwell=temp;
    f.Dwell(1)=2 * f.Dwell(2) - f.Dwell(3);
    f.Dwell(f.ns)=2 * f.Dwell(f.ns-1) - f.Dwell(f.ns-2);
    temp(2:f.ns-1)=f.Dcurr;
    f.Dcurr=temp;
    f.Dcurr(1)=2 * f.Dcurr(2) - f.Dcurr(3);
    f.Dcurr(f.ns)=2 * f.Dcurr(f.ns-1) - f.Dcurr(f.ns-2);
    temp(2:f.ns-1)=f.Dgeod;
    f.Dgeod=temp;
    f.Dgeod(1)=2 * f.Dgeod(2) - f.Dgeod(3);
    f.Dgeod(f.ns)=2 * f.Dgeod(f.ns-1) - f.Dgeod(f.ns-2);
end
% Now do the matrix values
% First Index (note indexing on vectors is 2:ns when VMEC outputs)
f.lmns(:,1)=     1.5 *     f.lmns(:,2) - 0.5 *     f.lmns(:,3);
f.bsupumnc(:,1)= 1.5 * f.bsupumnc(:,2) - 0.5 * f.bsupumnc(:,3);
f.bsupvmnc(:,1)= 1.5 * f.bsupvmnc(:,2) - 0.5 * f.bsupvmnc(:,3);
f.bsubsmns(:,1)= 1.5 * f.bsubsmns(:,2) - 0.5 * f.bsubsmns(:,3);
f.bsubumnc(:,1)= 1.5 * f.bsubumnc(:,2) - 0.5 * f.bsubumnc(:,3);
f.bsubvmnc(:,1)= 1.5 * f.bsubvmnc(:,2) - 0.5 * f.bsubvmnc(:,3);
f.gmnc(:,1)=     1.5 *     f.gmnc(:,2) - 0.5 *     f.gmnc(:,3);
f.bmnc(:,1)=     1.5 *     f.bmnc(:,2) - 0.5 *     f.bmnc(:,3);
% Average
for i=2:f.ns-1
    f.lmns(:,i)=     0.5 * (     f.lmns(:,i) +     f.lmns(:,i+1) );
    f.bsupumnc(:,i)= 0.5 * ( f.bsupumnc(:,i) + f.bsupumnc(:,i+1) );
    f.bsupvmnc(:,i)= 0.5 * ( f.bsupvmnc(:,i) + f.bsupvmnc(:,i+1) );
    f.bsubsmns(:,i)= 0.5 * ( f.bsubsmns(:,i) + f.bsubsmns(:,i+1) );
    f.bsubumnc(:,i)= 0.5 * ( f.bsubumnc(:,i) + f.bsubumnc(:,i+1) );
    f.bsubvmnc(:,i)= 0.5 * ( f.bsubvmnc(:,i) + f.bsubvmnc(:,i+1) );
    f.gmnc(:,i)    = 0.5 * (     f.gmnc(:,i) +     f.gmnc(:,i+1) );
    f.bmnc(:,i)    = 0.5 * (     f.bmnc(:,i) +     f.bmnc(:,i+1) );
end
% Last Index (note indexing on vectors is 2:ns when VMEC outputs)
f.lmns(:,f.ns)=     2.0 *     f.lmns(:,f.ns-1) -     f.lmns(:,f.ns-2);
f.bsupumnc(:,f.ns)= 2.0 * f.bsupumnc(:,f.ns-1) - f.bsupumnc(:,f.ns-2);
f.bsupvmnc(:,f.ns)= 2.0 * f.bsupvmnc(:,f.ns-1) - f.bsupvmnc(:,f.ns-2);
f.bsubsmns(:,f.ns)= 2.0 * f.bsubsmns(:,f.ns-1) - f.bsubsmns(:,f.ns-2);
f.bsubumnc(:,f.ns)= 2.0 * f.bsubumnc(:,f.ns-1) - f.bsubumnc(:,f.ns-2);
f.bsubvmnc(:,f.ns)= 2.0 * f.bsubvmnc(:,f.ns-1) - f.bsubvmnc(:,f.ns-2);
f.gmnc(:,f.ns)=     2.0 *     f.gmnc(:,f.ns-1) -     f.gmnc(:,f.ns-2);
f.bmnc(:,f.ns)=     2.0 *     f.bmnc(:,f.ns-1) -     f.bmnc(:,f.ns-2);
% Handle ANI/FLOW Values
if isfield(f,'prprmnc')
    f.prprmnc(:,1) =     1.5 *     f.prprmnc(:,2) - 0.5 *     f.prprmnc(:,3);
    for i=2:f.ns-1
        f.prprmnc(:,i) =     0.5 * (     f.prprmnc(:,i) +     f.prprmnc(:,i+1) );
    end
    f.prprmnc(:,f.ns)=     2.0 *     f.prprmnc(:,f.ns-1) -     f.prprmnc(:,f.ns-2);
end
if isfield(f,'protmnc')
    f.protmnc(:,1) =     1.5 *     f.protmnc(:,2) - 0.5 *     f.protmnc(:,3);
    for i=2:f.ns-1
        f.protmnc(:,i) =     0.5 * (     f.protmnc(:,i) +     f.protmnc(:,i+1) );
    end
    f.protmnc(:,f.ns)=     2.0 *     f.protmnc(:,f.ns-1) -     f.protmnc(:,f.ns-2);
end
if isfield(f,'protrsqmnc')
    f.protrsqmnc(:,1) =     1.5 *     f.protrsqmnc(:,2) - 0.5 *     f.protrsqmnc(:,3);
    for i=2:f.ns-1
        f.protrsqmnc(:,i) =     0.5 * (     f.protrsqmnc(:,i) +     f.protrsqmnc(:,i+1) );
    end
    f.protrsqmnc(:,f.ns)=     2.0 *     f.protrsqmnc(:,f.ns-1) -     f.protrsqmnc(:,f.ns-2);
end
if f.iasym >0 % Handle existance of lmnc on half mesh
    f.lmnc(:,1)=     1.5 *     f.lmnc(:,2) - 0.5 *     f.lmnc(:,3);
    f.bsupumns(:,1)= 1.5 * f.bsupumns(:,2) - 0.5 * f.bsupumns(:,3);
    f.bsupvmns(:,1)= 1.5 * f.bsupvmns(:,2) - 0.5 * f.bsupvmns(:,3);
    f.bsubsmnc(:,1)= 1.5 * f.bsubsmnc(:,2) - 0.5 * f.bsubsmnc(:,3);
    f.bsubumns(:,1)= 1.5 * f.bsubumns(:,2) - 0.5 * f.bsubumns(:,3);
    f.bsubvmns(:,1)= 1.5 * f.bsubvmns(:,2) - 0.5 * f.bsubvmns(:,3);
    f.gmns(:,1)=     1.5 *     f.gmns(:,2) - 0.5 *     f.gmns(:,3);
    f.bmns(:,1)=     1.5 *     f.bmns(:,2) - 0.5 *     f.bmns(:,3);
    for i=2:f.ns-1
        f.lmnc(:,i)=     0.5 * (     f.lmnc(:,i) +     f.lmnc(:,i+1) );
        f.bsupumns(:,i)= 0.5 * ( f.bsupumns(:,i) + f.bsupumns(:,i+1) );
        f.bsupvmns(:,i)= 0.5 * ( f.bsupvmns(:,i) + f.bsupvmns(:,i+1) );
        f.bsubsmnc(:,i)= 0.5 * ( f.bsubsmnc(:,i) + f.bsubsmnc(:,i+1) );
        f.bsubumns(:,i)= 0.5 * ( f.bsubumns(:,i) + f.bsubumns(:,i+1) );
        f.bsubvmns(:,i)= 0.5 * ( f.bsubvmns(:,i) + f.bsubvmns(:,i+1) );
        f.gmns(:,i)    = 0.5 * (     f.gmns(:,i) +     f.gmns(:,i+1) );
        f.bmns(:,i)    = 0.5 * (     f.bmns(:,i) +     f.bmns(:,i+1) );
    end
    f.lmnc(:,f.ns)=     2.0 *     f.lmnc(:,f.ns-1) -     f.lmnc(:,f.ns-2);
    f.bsupumns(:,f.ns)= 2.0 * f.bsupumns(:,f.ns-1) - f.bsupumns(:,f.ns-2);
    f.bsupvmns(:,f.ns)= 2.0 * f.bsupvmns(:,f.ns-1) - f.bsupvmns(:,f.ns-2);
    f.bsubsmnc(:,f.ns)= 2.0 * f.bsubsmnc(:,f.ns-1) - f.bsubsmnc(:,f.ns-2);
    f.bsubumns(:,f.ns)= 2.0 * f.bsubumns(:,f.ns-1) - f.bsubumns(:,f.ns-2);
    f.bsubvmns(:,f.ns)= 2.0 * f.bsubvmns(:,f.ns-1) - f.bsubvmns(:,f.ns-2);
    f.gmns(:,f.ns)=     2.0 *     f.gmns(:,f.ns-1) -     f.gmns(:,f.ns-2);
    f.bmns(:,f.ns)=     2.0 *     f.bmns(:,f.ns-1) -     f.bmns(:,f.ns-2);
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = h2f(var,ns)
% Map quantitiy from half to full grid
temp=zeros(1,ns);
temp(1)  =  1.5 *   var( 1) - 0.5 *    var(2);
temp(2:ns-1)=0.5 * ( var(1:ns-2) + var(2:ns-1));
%for i=2:ns-1
%    temp(i)= 0.5 * ( var(i) +    var(i+1) );
%end
temp(ns) =  1.5 *   var(ns-1) - 0.5 * var(ns-2);
f=temp;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=read_vmec_netcdf(filename)
mu0=4.*pi.*1e-7;
f=read_netcdf(filename,'strip','flipdim');
% Now fix named fields so they match the old way of doing things
f.input_extension=f.inputextension;
f.mgrid_file=char(f.mgridfile);
f.ierr_vmec=f.ierflag;
f.rmax_surf=f.rmaxsurf;
f.rmin_surf=f.rminsurf;
f.zmax_surf=f.zmaxsurf;
f.ireconstruct=f.lreconlogical;
f.imse=-1;
f.itse=-1;
f.RBtor=f.rbtor;
f.Rmajor=f.Rmajorp;
f.Aminor=f.Aminorp;
f.betatot=f.betatotal;
f.Volume=f.volumep;
f.VolAvgB=f.volavgB;
if isfield(f,'betavol')
    f.beta_vol=f.betavol';
end
if isfield(f,'specw')
    f.specw=f.specw';
end
if ~isfield(f,'iasym')
    f.iasym=0;
end
f.iasym=f.lasymlogical;
f.freeb=f.lfreeblogical;
f.lfreeb=f.freeb;
f.Itor=f.ctor;
f.Dmerc=f.DMerc;
f.Dwell=f.DWell;
f.Dshear=f.DShear;
f.Dcurr=f.DCurr;
f.Dgeod=f.DGeod;
% Cast some values
f.ntor=double(f.ntor);
f.mpol=double(f.mpol);
f.nfp=double(f.nfp);
f.ns=double(f.ns);
% Calculate Currents
f.currumnc=zeros(f.mnmaxnyq,f.ns);
f.currvmnc=zeros(f.mnmaxnyq,f.ns);
ohs = f.ns-1;
hs  = 1.0/double(ohs);
ns = f.ns;
for i=2:ns
    shalf(i) = sqrt(hs*(i-1.5));
    sfull(i) = sqrt(hs*(i-1));
end
js1 = 3:ns;
js  = 2:(ns-1);
for mn = 1:f.mnmaxnyq
    if (mod(f.xmnyq,2) == 1)
        t1  = 0.5.*(shalf(js1).*f.bsubsmns(mn,js1)+...
            shalf(js).*f.bsubsmns(mn,js))./sfull(js);
        bu0 = f.bsubumnc(mn,js)./shalf(js);
        bu1 = f.bsubumnc(mn,js1)./shalf(js1);
        t2  = ohs.*(bu1-bu0).*sfull(js)+0.25.*(bu0+bu1)./sfull(js);
        bv0 = f.bsubvmnc(mn,js)./shalf(js);
        bv1 = f.bsubvmnc(mn,js1)./shalf(js1);
        t3  = ohs.*(bv1-bv0).*sfull(js)+0.25.*(bv0+bv1)./sfull(js);
    else
        t1  = 0.5.*(f.bsubsmns(mn,js1)+f.bsubsmns(mn,js));
        t2  = ohs.*(f.bsubumnc(mn,js1)-f.bsubumnc(mn,js));
        t3  = ohs.*(f.bsubvmnc(mn,js1)-f.bsubvmnc(mn,js));
    end
    f.currumnc(mn,js) = -double(f.xnnyq(mn)).*t1 - t3;
    f.currvmnc(mn,js) = -double(f.xmnyq(mn)).*t1 + t2;
end
% OLD Way
%for i=2:f.ns-1
%    f.currumnc(:,i)=-double(f.xnnyq)'.*f.bsubsmns(:,i)-(f.ns-1).*(f.bsubvmnc(:,i+1)-f.bsubvmnc(:,i));
%    f.currvmnc(:,i)=-double(f.xmnyq)'.*f.bsubsmns(:,i)+(f.ns-1).*(f.bsubumnc(:,i+1)-f.bsubumnc(:,i));
%end
f.currumnc(:,1)=0.0;
f.currvmnc(:,1)=0.0;
for i=1:f.mnmaxnyq
    if (f.xmnyq(i)<=1)
        f.currumnc(i,1)=2.*f.currumnc(i,2)-f.currumnc(i,3);
        f.currvmnc(i,1)=2.*f.currvmnc(i,2)-f.currvmnc(i,3);
    end
end
f.currumnc(:,f.ns)=2.*f.currumnc(:,f.ns-1)-f.currumnc(:,f.ns-2);
f.currvmnc(:,f.ns)=2.*f.currvmnc(:,f.ns-1)-f.currvmnc(:,f.ns-2);
f.currumnc=f.currumnc./mu0;
f.currvmnc=f.currvmnc./mu0;
if f.iasym
    f.currumns=zeros(f.mnmaxnyq,f.ns);
    f.currvmns=zeros(f.mnmaxnyq,f.ns);
    for mn = 1:f.mnmaxnyq
        if (mod(f.xmnyq,2) == 1)
            t1  = 0.5.*(shalf(js1).*f.bsubsmnc(mn,js1)+...
                shalf(js).*f.bsubsmnc(mn,js))./sfull(js);
            bu0 = f.bsubumns(mn,js)./shalf(js);
            bu1 = f.bsubumns(mn,js1)./shalf(js1);
            t2  = ohs.*(bu1-bu0).*sfull(js)+0.25.*(bu0+bu1)./sfull(js);
            bv0 = f.bsubvmns(mn,js)./shalf(js);
            bv1 = f.bsubvmns(mn,js1)./shalf(js1);
            t3  = ohs.*(bv1-bv0).*sfull(js)+0.25.*(bv0+bv1)./sfull(js);
        else
            t1  = 0.5.*(f.bsubsmnc(mn,js1)+f.bsubsmnc(mn,js));
            t2  = ohs.*(f.bsubumns(mn,js1)-f.bsubumns(mn,js));
            t3  = ohs.*(f.bsubvmns(mn,js1)-f.bsubvmns(mn,js));
        end
        f.currumns(mn,js) = double(f.xnnyq(mn)).*t1 - t3;
        f.currvmns(mn,js) = double(f.xmnyq(mn)).*t1 + t2;
    end
    % OLD WAY
    %for i=2:f.ns-1
    %    f.currumns(:,i)= double(f.xnnyq)'.*f.bsubsmnc(:,i)-(f.ns-1).*(f.bsubvmns(:,i+1)-f.bsubvmns(:,i));
    %    f.currvmns(:,i)= double(f.xmnyq)'.*f.bsubsmnc(:,i)+(f.ns-1).*(f.bsubumns(:,i+1)-f.bsubumns(:,i));
    %end
    f.currumns(:,1)=0.0;
    f.currvmns(:,1)=0.0;
    for i=1:f.mnmaxnyq
        if (f.xmnyq(i)<=1)
            f.currumns(i,1)=2.*f.currumns(i,2)-f.currumns(i,3);
            f.currvmns(i,1)=2.*f.currvmns(i,2)-f.currvmns(i,3);
        end
    end
    f.currumns(:,f.ns)=2.*f.currumns(:,f.ns-1)-f.currumns(:,f.ns-2);
    f.currvmns(:,f.ns)=2.*f.currvmns(:,f.ns-1)-f.currvmns(:,f.ns-2);
    f.currumns=f.currumns./mu0;
    f.currvmns=f.currvmns./mu0;
end
% Remove Renamed Fields
f=rmfield(f,'inputextension');
f=rmfield(f,'mgridfile');
f=rmfield(f,'ierflag');
f=rmfield(f,'Rmajorp');
f=rmfield(f,'Aminorp');
f=rmfield(f,'betatotal');
f=rmfield(f,'volumep');
f=rmfield(f,'volavgB');
f=rmfield(f,'betavol');
f=rmfield(f,'lasymlogical');
f=rmfield(f,'lfreeblogical');
f=rmfield(f,'lreconlogical');
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=read_vmec_mercier(filename)
fid=fopen(filename,'r');
fgetl(fid); % First Header
fgetl(fid); % Second Header ----
% Read first line
line=fgetl(fid);
val=sscanf(line,'%e');
f.s=val(1);
f.phi=val(2);
f.iota=val(3);
f.shear=val(4);
f.vp=val(5);
f.well=val(6);
f.itor=val(7);
f.ditor=val(8);
f.pres=val(9);
f.dpres=val(10);
line=fgetl(fid);
while ~strcmp(line,'');
    val=sscanf(line,'%e');
    f.s=[f.s; val(1)];
    f.phi=[f.phi; val(2)];
    f.iota=[f.iota; val(3)];
    f.shear=[f.shear; val(4)];
    f.vp=[f.vp; val(5)];
    f.well=[f.well; val(6)];
    f.itor=[f.itor; val(7)];
    f.ditor=[f.ditor; val(8)];
    f.pres=[f.pres; val(9)];
    f.dpres=[f.dpres; val(10)];
    line=fgetl(fid);
end
fgetl(fid);
fgetl(fid);
% Read first line
line=fgetl(fid);
val=sscanf(line,'%e');
f.dmerc=val(2);
f.dshear=val(3);
f.dcurr=val(4);
f.dwell=val(5);
f.dgeod=val(6);
while ~feof(fid);
    line=fgetl(fid);
    val=sscanf(line,'%e');
    f.dmerc=[f.dmerc; val(2)];
    f.dshear=[f.dshear; val(3)];
    f.dcurr=[f.dcurr; val(4)];
    f.dwell=[f.dwell; val(5)];
    f.dgeod=[f.dgeod; val(6)];
end
fclose(fid);
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=read_vmec_jxbout(filename)
fid=fopen(filename,'r');
fgetl(fid); % Blank Line
line=fgetl(fid); % Radial Poloidal Toroidal points
val=sscanf(line,'%*20c %d %*24c %d %*24c %d');
f.nrad=val(1);
f.ntheta=val(2);
f.nzeta=val(3);
line=fgetl(fid) % mpol ntor
val=sscanf(line,'%*18c %d %*18c %d');
f.mpol=val(1);
f.ntor=val(2);
for i=1:13      %Header stuff
    fgetl(fid);
end
line=fgetl(fid);    % Values
val=sscanf(line,'%*18c %e %*17c %e %*11c %e %*17c %e');
f.tor_flux=val(1);
f.fnorm=val(2);
f.jdotb=val(3);
f.bdotv=val(4);
line=fgetl(fid);    %percentages
fgetl(fid);
fgetl(fid);
fgetl(fid);
f.data=zeros(f.nrad,f.nzeta,f.ntheta,13);
for i=1:nrad
    for j=1:f.nzeta
        line=fgetl(fid);    % Angle Information
        for k=1:f.ntheta
            line=fgetl(fid);    % Data
            f.data(i,j,k,:)=sscanf(line,'%d %12e');
        end
    end
end
fclose(fid);
return
end


