function data = read_netcdf(filename,varargin)
%data=READ_NETCDF(filename[,'flipdim','strip']) Extracts netcdf file to structure
%   This function extracts a netcdf file into a structure which is
%   returned to the command line.  You need only specify a filename of a
%   netCDF file.  This version uses the java class to read data directly
%   from the netCDF file.  All '+' are converted to 'p' and '-' to 'n' in
%   attribute, dimension, and variable names.  All '.' are converted to
%   underscores.
%
%   Options:
%   READ_NETCDF(filename,'flipdim') This will automatically permute the
%   arrays read from the file in reverse order.
%
%   READ_NETCDF(filename,'strip') This will strip non-alphanumerica
%   characters from element names to allow them to be used as structure
%   elements.
%
%   Note:
%   Some versions of MATLAB will require a local copy of the netcdfAll.jar
%   files.  To temporarily add these use:
%   javaaddpath('ftp://ftp.unidata.ucar.edu/pub/netcdf-java/v4.1/netcdfAll-4.1.jar');
%   Or you can download the files locally and add them to your classpath.txt
%   file.  You can find your classpath.txt file using the which command:
%   which classpath.txt
%   Then you may add the jar file permanently by using the command:
%   edit classpath.txt.  Also note you will need the netcdfAll java class
%   to use this routine.
%
%   Example
%       data=read_netcdf('input.nc');
%
%   Version 1.5
%   Maintained by: Samuel Lazerson (lazerson@pppl.gov)
%   Date  08/01/2011

% 02/15/11-SAL Updated to support multiple attribute types.  Thanks to 
% Zelalem Engida for bringing this to my attention.
%
% 04/12/11-SAL All '+' converted to 'p' in names.  All '-' converted to 'n'
% in names.  Thanks to Thomas Leslie Leith for bringing this to my
% attention.
%
% 08/01/11-SAL All '.' converted to '_' in names.  Code now checks to
% ensure that the maximum name length is not exceeded and that the
% variables names do not begin with non-alphabet characters.

% Allow the user to pass some variables.
flipdim=0;
strip=0;
netcdffound=[];
maxlength=namelengthmax;
if nargin>1
    for i=2:nargin
        switch varargin{i-1}
            case 'flipdim'
                flipdim=1;
            case 'strip'
                strip=1;
        end
    end
end
% Test for netCDF-File
loadedclasses=javaclasspath('-all');
for i=1:length(loadedclasses)
    netcdffound=strfind(loadedclasses{i},'netcdfAll');
    if ~isempty(netcdffound)
        break
    end
end
if ~isempty(netcdffound)
else
    disp('ERROR:  Could not find netcdfALL JAVA class in your loaded classpaths');
    disp('        Please try:');
    disp('           help read_netcdf');
    disp('        for instructions on loading the proper JAVA class.');
    return
end
import ucar.nc2.*;
% Open the File
try
    netcdfdata=NetcdfFile.open(filename);
catch read_netcdf_error
    data=-1;
    disp(['ERROR: Opening netCDF File: ' filename]);
    disp(['  -identifier: ' read_netcdf_error.identifier]);
    disp(['  -message:    ' read_netcdf_error.message]);
    disp('      For information type:  help read_netcdf');
    return
end
% Get information on number of elements
ndimen=netcdfdata.getDimensions.size;
nvars=netcdfdata.getVariables.size;
ngatts=netcdfdata.getGlobalAttributes.size;
unlimidimid=netcdfdata.getUnlimitedDimension;
% Get global attributes
if ngatts >0
    for i=0:ngatts-1
        att=netcdfdata.getGlobalAttributes.get(i);
        attname=char(att.getName);
        % Do some checking
        if ~isletter(attname(1))
            disp(['Error Attribute not read: ' attname '  Attributes must begin with a character.']);
            break
        end
        if length(attname) > maxlength
            disp(['Warning Variable: ' attname ' too long.  Max structure name length:' num2str(maxlength)]');
        end
        attname(attname=='.')='_';
        attname(attname=='-')='n';
        attname(attname=='+')='p';
        if strip, attname(~isstrprop(attname,'alphanum'))=''; end;
        % We handle multiple data types
        switch char(att.getDataType)
            case 'int'
                data.(attname)=att.getNumericValue.intValue;
            case 'float'
                data.(attname)=att.getNumericValue.floatValue;
            case 'double'
                data.(attname)=att.getNumericValue.doubleValue;
            case 'short'
                data.(attname)=att.getNumericValue.shortValue;
            case 'long'
                data.(attname)=att.getNumericValue.longValue;
            case 'byte'
                data.(attname)=att.getNumericValue.byteValue;
            case 'String'
                data.(attname)=char(att.getStringValue);
            otherwise
                disp(['-----Error reading: ' attname]);
        end
    end
end
% Get dimensions
if ndimen >0
    for i=0:ndimen-1
        dim=netcdfdata.getDimensions.get(i);
        dimname=char(dim.getName);
        % Do some checking
        if ~isletter(dimname(1))
            disp(['Error Dimension not read: ' dimname '  Dimensions must begin with a character.']);
            break
        end
        if length(dimname) > maxlength
            disp(['Warning Dimension Name: ' dimname ' too long.  Max structure name length:' num2str(maxlength)]');
        end
        attname(dimname=='.')='_';
        dimname(dimname=='-')='n';
        dimname(dimname=='+')='p';
        if strip, dimname(~isstrprop(dimname,'alphanum'))=''; end;
        data.(dimname)=dim.getLength;
    end
end
% Get Variables
if nvars >0
    for i=0:nvars-1
        var=netcdfdata.getVariables.get(i);
        varname=char(var.getName);
        % Do some checking
        if ~isletter(varname(1))
            disp(['Error Variable not read: ' varname '  Dimensions must begin with a character.']);
            break
        end
        if length(varname) > maxlength
            disp(['Warning Variable Name: ' varname ' too long.  Max structure name length:' num2str(maxlength)]');
        end
        attname(varname=='.')='_';
        varname(varname=='-')='n';
        varname(varname=='+')='p';
        if strip, varname(~isstrprop(varname,'alphanum'))=''; end;
        natts=var.getAttributes.size;
        if (var.getDataType.isString)
            try
                data.(varname)=char(var.getStringValue);
            catch exception
                try
                    data.(varname)=char(var.readScalarString);
                catch exception
                    try
                        data.(varname)=char(var.read);
                    catch exception
                        disp(['-----Error reading: ' varname]);
                    end
                end
            end
        elseif (var.getDataType.isNumeric)
            if var.getSize == 1
                data.(varname)=var.read.copyTo1DJavaArray;
            else
                if flipdim
                    temp=var.read.copyToNDJavaArray;
                    temp=permute(temp,ndims(temp):-1:1);
                    data.(varname)=temp;
                else
                    data.(varname)=var.read.copyToNDJavaArray;
                end
            end
        end
        if natts > 0
            for j=0:natts-1
                att=var.getAttributes.get(j);
                attname=char(att.getName);
                if strip, attname(~isstrprop(attname,'alphanum'))=''; end;
                varattname=strcat(varname,'_',attname);
                switch char(att.getDataType)
                    case 'int'
                        data.(varattname)=att.getNumericValue.intValue;
                    case 'float'
                        data.(varattname)=att.getNumericValue.floatValue;
                    case 'double'
                        data.(varattname)=att.getNumericValue.doubleValue;
                    case 'short'
                        data.(varattname)=att.getNumericValue.shortValue;
                    case 'long'
                        data.(varattname)=att.getNumericValue.longValue;
                    case 'byte'
                        data.(varattname)=att.getNumericValue.byteValue;
                    case 'String'
                        data.(varattname)=char(att.getStringValue);
                    otherwise
                        disp(['-----Error reading: ' varattname]);
                end
            end
        end
    end
end
% Close the file
netcdfdata.close;
end

