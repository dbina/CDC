function pdb2mre_chl
% reads text data from *.pdb files
% the residues are identified by the variable 'reslabel'
% output as of 2021-04-28 consists of x-y-z coordinate letter label (for
% chain) and number index of atom from original *.pdb file
% ver. 1.2 2022-Oct-14 (c)DB
% last major modification: 2023-II-06

clear all; 
%close all; 

% file to analyze:
fid = fopen( 'C:\Data\Progs\dip_models\dip_path.txt' );
% fid = fopen( 'C:\DataVB\Progs\dip_models\dip_path.txt' );
fp0 = fgetl(fid);
fclose(fid);     

disp( 'starting tres_reader, selecting the tres data folder')
% get the directory and read data
directory_name = uigetdir( fp0 )
% make the list of files
d = dir( directory_name );
dsize = max( size( {d.name} ) );
str0 = {d.name};

jj = 1;
for ii = 1:dsize
     
   if ~isempty( findstr( str0{ii}, '.pdb' ) )
       str{jj} = str0{ii};
       jj = jj + 1;
   end

end
[s,v] = listdlg('PromptString','Select a file:',...
                 'SelectionMode','multiple',...
                 'ListString',str);

reslabel = 'CLA';    

% now plot all atoms of the molecule
IFPath2 = directory_name;
IFPath20 = directory_name;

if strcmp( reslabel, 'HETATM') | strcmp( reslabel, 'ATOM' )
    
    if strcmp( reslabel, 'HETATM')
        colscale = 1;
        symb = '.';
    end
    if strcmp( reslabel, 'ATOM')
        colscale = 2;
        symb = '.';
    end
    
    for ii = 1:length( s ) 
        dfile = str{s(ii)};

        C = read_PDB_coord_atoms( [IFPath2,'\',dfile] , 1, reslabel );    
        C(:, 1:3) = C(:, 1:3) / 10;  % from A to nm
                       
        for jj = 1:size( C , 1 )
            
            switch C(jj, 5)
                case 1
                    poco = [0, 1, 0]; % C
                case 2
                    poco = [1, 0, 0]; % O
                case 3
                    poco = [0, 0, 1]; % N
                case 4
                    poco = [1, 1, 0]; % S
            end
            plot3( C(jj,1), C(jj,2), C(jj,3), symb, 'Color', poco / colscale,...
                                                    'MarkerFaceColor', poco / colscale,...
                                                    'MarkerEdgeColor', poco / colscale )
            hold on
        end        
    end
    set(gca,'DataAspectRatio', [1 1 1])
    
    return 
end

% this part is for dipole extraction
IFPath2 = directory_name;

% my_color = [ [0.32465, 0.96923, 0.60653]*0.8;
%              [0.75484, 0.8825, 0.21627]*0.8;
%              [1, 0.45783, 0.043937]; 
%              [0.07956     0.60653     0.96923];
%              [1, 0, 0] ;
%              [1, 0.1, 0.4]; 
%              [1, 0.3, 0.4];
%              [0.5, 1, 0.4]; 
%              [0.1, 0.3, 0.4];                      
%          ];        
colorcount = 8;
for ii = 1:colorcount
    
    my_color(ii,:) = myrainbow( ii, 0, colorcount+1, 1 ) ;
    
end

for ii = 1:length( s ) 
    
    CM = [];  
    dfile = str{s(ii)};
    dfile0 = dfile;
    C = read_PDB_coord_02( [IFPath2,'\',dfile] , 1, reslabel )
    C(:, 1:3) = C(:, 1:3) / 10; % angstrom to nm
%     disp( C )
%     size(C)          
    mRe = zeros( size(C, 1) / 5, 7 );
    mRex = zeros( size(C, 1) / 5, 7 );
   
%     figure(1)
    
    for jj = 1:size( C , 1) / 5
        
        block = (jj-1)*5+(1:5);
        u = C( block, 4 );        
        num = u(1);                
        
        index(jj,1) = 1;
% parametrizations        
% %.....................Chla-only LHCII ..............
%         switch num
%             case 602
%                  index(jj,1) = 1;
%             case 603
%                  index(jj,1) = 2;
%             case 604
%                  index(jj,1) = 3;
%             case 610
%                  index(jj,1) = 4;
%             case 611
%                  index(jj,1) = 5;
%             case 612
%                  index(jj,1) = 6;
%             case 613
%                  index(jj,1) = 7;
%             case 614
%                  index(jj,1) = 8;
%          end            
%...................................................        
% cyanidioschyzon 5zgh
%         if num < 800
%             index( jj, 1) = 4;
%         end
%         if jj >= 124 & jj <= 135
%             index( jj, 1) = 1;
%         end        
%          switch jj
%             case 1
%                  index(jj,1) = 3;
%             case 45
%                  index(jj,1) = 3;
%             case 40
%                  index(jj,1) = 2;
%             case 41
%                  index(jj,1) = 2; 
%             case 84
%                  index(jj,1) = 3;
%             case 38
%                  index(jj,1) = 3;
%             case 80
%                  index(jj,1) = 2;
%             case 77
%                  index(jj,1) = 2; 
%             case 119
%                  index(jj,1) = 6;
%             case 76
%                  index(jj,1) = 6;
%             case 52
%                  index(jj,1) = 7;
%             case 50
%                  index(jj,1) = 7; 
%             case 82
%                  index(jj,1) = 7;
%             case 79
%                  index(jj,1) = 7; 
%          end    
                 
% synechococcus elongatus 1JB0:
%        switch jj
%             case 1
%                  index(jj,1) = 3;
%             case 47
%                  index(jj,1) = 3;
%             case 3
%                  index(jj,1) = 2;
%             case 49
%                  index(jj,1) = 2;                 
%             case 80
%                  index(jj,1) = 4;                 
%             case 81
%                  index(jj,1) = 4;                 
%             case 82
%                  index(jj,1) = 4;                 
%             case 41
%                  index(jj,1) = 5;                 
%             case 42
%                  index(jj,1) = 5;                 
%             case 56
%                  index(jj,1) = 5;                 
%             case 35
%                  index(jj,1) = 5;                 
%             case 34
%                  index(jj,1) = 5; 
%         end
%....................................................
         
        disp( [index(jj,1) , num] )
        
        nax = find( C( block, 5 ) == 1 );
        nbx = find( C( block, 5 ) == 2 );
        ncx = find( C( block, 5 ) == 3 );
        ndx = find( C( block, 5 ) == 4 );    
        
        mgx = find( C( block, 5 ) == 6 )   
        
        xyz = C( block, 1:3 );
        na = xyz( nax, : );
        nb = xyz( nbx, : );
        nc = xyz( ncx, : );        
        nd = xyz( ndx, : );
        
        mg = xyz( mgx, : );
        
%        center = mean( [xyz; mg], 1 );
        center = mean( xyz, 1 ) % place center to midpoint of NA-NB-NC-ND quartet
%         center = mg;
%         norm( mg - center )
        
        qy = nb - nd;
        qx = na - nc;        

        qy = qy / norm( qy );
        qx = qx / norm( qx );
% Chla dipole
        mRe(jj,1) = 4;
        mRex(jj,1) = 1.6;        

% % Chlb dipole
%         mRe(jj,1) = 3.4;
%         mRex(jj,1) = 1.4;        
%         
% Chlc dipole
%         mRe(jj,1) = 2.3;
%         mRex(jj,1) = 2;        

        mRe(jj,2:4) = center;        
        mRex(jj,2:4) = center;              
        mRe(jj,5:7) = qy;        
        mRex(jj,5:7) = qx;        
        
        plot3( na(1), na(2), na(3), 'ro', 'MarkerFaceColor', 'r' )
        hold on
        plot3( nb(1), nb(2), nb(3), 'bo', 'MarkerFaceColor', 'b' )
        plot3( nc(1), nc(2), nc(3), 'ko', 'MarkerFaceColor', 'k' )
        plot3( nd(1), nd(2), nd(3), 'co', 'MarkerFaceColor', 'c' )    
        
        plot3( [mRe(jj,2)-mRe(jj,5)/5, mRe(jj,2)+mRe(jj,5)/5],...
               [mRe(jj,3)-mRe(jj,6)/5, mRe(jj,3)+mRe(jj,6)/5],...
               [mRe(jj,4)-mRe(jj,7)/5, mRe(jj,4)+mRe(jj,7)/5], 'm-', 'LineWidth', 2 ) 
        text( mRe(jj,2)+0.05, mRe(jj,3)+0.05, mRe(jj,4)+0.05, [num2str( num ),'...', num2str(jj) ] ) 
        hold on
        na = 2.125*(na - center);
        nb = 2.125*(nb - center);
        nc = 2.125*(nc - center);
        nd = 2.125*(nd - center);        
        
        plot3( [mRe(jj,2)+na(1), mRe(jj,2)+nb(1)],...
               [mRe(jj,3)+na(2), mRe(jj,3)+nb(2)],...
               [mRe(jj,4)+na(3), mRe(jj,4)+nb(3)], '-', 'LineWidth', 2, 'Color', my_color( index(jj) , :))

        plot3( [mRe(jj,2)+nb(1), mRe(jj,2)+nc(1)],...
               [mRe(jj,3)+nb(2), mRe(jj,3)+nc(2)],...
               [mRe(jj,4)+nb(3), mRe(jj,4)+nc(3)], '-', 'LineWidth', 2, 'Color', my_color( index(jj) , :))
           
        plot3( [mRe(jj,2)+nc(1), mRe(jj,2)+nd(1)],...
               [mRe(jj,3)+nc(2), mRe(jj,3)+nd(2)],...
               [mRe(jj,4)+nc(3), mRe(jj,4)+nd(3)], '-', 'LineWidth', 2, 'Color', my_color( index(jj) , :))
           
        plot3( [mRe(jj,2)+nd(1), mRe(jj,2)+na(1)],...
               [mRe(jj,3)+nd(2), mRe(jj,3)+na(2)],...
               [mRe(jj,4)+nd(3), mRe(jj,4)+na(3)], '-', 'LineWidth', 2, 'Color', my_color( index(jj) , :))           
        
        mRe_out(jj,:) = [jj, num, index(jj), mRe(jj,:)]; % the output is [number, residue label, index(=pigment class), dipole, coordinates, dipole vector components]
        mRex_out(jj,:) = [jj, num, index(jj), mRex(jj,:)]; % the output is [number, residue label, index(=pigment class), dipole, coordinates, dipole vector components]
        
      
    end
    
%     plot3( C(:,1), C(:,2), C(:,3), '.' )
    hold on
    plot3( mRe(:,2), mRe(:,3), mRe(:,4), 'rsq' )
    set(gca,'DataAspectRatio', [1 1 1])
    grid on
    
    dfile0 = dfile;
    IFPath20 = IFPath2;
    [dfile,IFPath2] = uiputfile([IFPath2,'\',dfile(1:end-4),'_', reslabel,'_vecdata_qy', '.txt'], 'vector file') 
    if dfile ~= 0
        save( [IFPath2,'\',dfile], 'mRe_out', '-ascii') 
        save( [IFPath2,'\',dfile(1:end-6), 'qx.txt'], 'mRex_out', '-ascii')        
    end
  
    [dfile,IFPath2] = uiputfile([IFPath20,'\index_',dfile0(1:end-4),'_', reslabel,'_vecdata','.txt'], 'index file') 
    if dfile ~= 0
        save( [IFPath2,'\',dfile], 'index', '-ascii') 
    end
    
%     if ~isempty( C )
%         na = length(C(:,1))
%         if length( s) == 1
%            [dfile,IFPath2] = uiputfile([IFPath2,'\',dfile(1:end-4),'_', reslabel, '.xyz']) 
%            if dfile == 0
%                break
%            else
%                save( [IFPath2,'\',dfile], 'C', '-ascii') 
%            end
%         else
%            save( [IFPath2,'\',dfile(1:end-4),'_', reslabel, '.xyz'], 'C', '-ascii')
%         end
%     else
%         disp( ['this file does not contain residue ',reslabel,', nothing saved' ]  )
%     end
    
end

% C = read_PDB_coord_atoms( [IFPath20,'\',dfile0] , 1, 'HETATM' );    
% C(:, 1:3) = C(:, 1:3) / 10;  % from A to nm
% colscale = 1;
% for jj = 1:size( C , 1 )
%     
%     switch C(jj, 5)
%         case 1
%             poco = [0, 1, 0]; % C
%         case 2
%             poco = [1, 0, 0]; % O
%         case 3
%             poco = [0, 0, 1]; % N
%         case 4
%             poco = [1, 1, 0]; % S
%     end
%     plot3( C(jj,1), C(jj,2), C(jj,3), '.', 'Color', poco / colscale,...
%         'MarkerFaceColor', poco / colscale,...
%         'MarkerEdgeColor', poco / colscale )
%     hold on
% end        
% end
% set(gca,'DataAspectRatio', [1 1 1])


disp('end OK' )
%==========================================================
function XYZ = read_PDB_coord_02(FName, startL, reslabel)
% opens and analyses a .PDB file 
% version as of Oct 2015
%
fid = fopen(FName);
disp(['extracting coordinates from file ', FName, ', please wait...'])
a = textread( FName,'%s','delimiter', '\n');
fclose(fid);

XYZ = [];
ii = 1;
jj = 1;
~strcmp('all', reslabel )

if (~isempty( findstr( a{ii} , 'HETATM') )) & (isempty( findstr( a{ii} , reslabel) ))
    disp(a{ii} )
end

if ~strcmp('all', reslabel )
    for ii = startL:size(a,1)-1
        if (~isempty( findstr( a{ii} , 'HETATM') )) & ( (~isempty( findstr( a{ii} , reslabel) ) ) | (~isempty( findstr( a{ii} , 'CL0') ) ))
% ATOM      1  N   SER C  14     -45.638 122.405  15.425  1.00 46.80           N                    
% HETATM41476  NA  CLA A 824     256.600 255.991 267.607  1.00 33.47           N  
% HETATM28407  O1D CHL 4 315      36.507 105.935 -33.192  1.00 96.36           O  

            chl_at_label = a{ii}( 13 : 17 );
            Nnum = 0;
            if ~isempty( findstr( chl_at_label , 'NA' ) )
                Nnum = 1;                
            end
            if ~isempty( findstr( chl_at_label , 'NB' ) )
                Nnum = 2;                
            end
            if ~isempty( findstr( chl_at_label , 'NC' ) )
                Nnum = 3;                
            end
            if ~isempty( findstr( chl_at_label , 'ND' ) )
                Nnum = 4;                
            end
            
            if ~isempty( findstr( chl_at_label , 'MG' ) )
                Nnum = 6;                
            end
            
%             disp( a{ii} )
            reslab = a{ii}( 23 : 29 );
            resnumlab = str2num( char( strread(reslab,'%s') ) );
%             disp( resnumlab );
%             chainlab = a{ii}(poslab+4);
%             atomindexinchain = a{ii}(poslab+5:30);
            
            XYZpart1 = a{ii}(31:38);
            XYZpart2 = a{ii}(39:47);
            XYZpart3 = a{ii}(48:56);        
            Coord1 = strread(XYZpart1,'%s');
            Coord2 = strread(XYZpart2,'%s');        
            Coord3 = strread(XYZpart3,'%s');        
            if ~isempty(Coord1) & Nnum ~= 0
                XYZ(jj,1) = str2num(char(Coord1));
                XYZ(jj,2) = str2num(char(Coord2));
                XYZ(jj,3) = str2num(char(Coord3));    
                XYZ(jj,4) = resnumlab; %double( chainlab);        
                XYZ(jj,5) = Nnum;%str2num( atomindexinchain );        
                
                jj = jj + 1;
            end
        end
        
        if ~isempty( findstr( a{ii} , 'CONNECT') )
            break
        end
        
    end
else
    for ii = startL:size(a,1)-1
        %     disp( a{ii} )
%         if ~isempty( findstr( a{ii}(1:7) , 'ATOM') ) | ~isempty( findstr( a{ii}(1:7), 'HETATM') )
         if ~isempty( findstr( a{ii} , 'HETATM') ) & ~isempty( findstr( a{ii}, 'CLA' ) )    
%               disp( a{ii} )
%             poslab = findstr( a{ii}, reslabel );
            chainlab = a{ii}(22);
            
            XYZpart1 = a{ii}(31:38);
            XYZpart2 = a{ii}(39:47);
            XYZpart3 = a{ii}(48:56);        
            Coord1 = strread(XYZpart1,'%s');
            Coord2 = strread(XYZpart2,'%s');        
            Coord3 = strread(XYZpart3,'%s');        
            if ~isempty(Coord1)
                XYZ(jj,1) = str2num(char(Coord1));
                XYZ(jj,2) = str2num(char(Coord2));
                XYZ(jj,3) = str2num(char(Coord3));    
                XYZ(jj,4) = double( chainlab);        
                
                jj = jj + 1;
            end
        end
        
        if ~isempty( findstr( a{ii} , 'CONNECT') )
            break
        end
    end    
end


disp('...done')

%------------------------------------------------------------
function XYZ = read_PDB_coord_atoms(FName, startL, reslabel)
% opens and analyses a .PDB file 
% version as of Oct 2015
%
fid = fopen(FName);
disp(['extracting coordinates from file ', FName, ', please wait...'])
a = textread( FName,'%s','delimiter', '\n');
fclose(fid);

XYZ = [];
ii = 1;
jj = 1;

if (~isempty( findstr( a{ii} , 'HETATM') )) & (isempty( findstr( a{ii} , reslabel) ))
    disp(a{ii} )
end

for ii = startL:size(a,1)-1
    if (~isempty( findstr( a{ii} , reslabel) )) 
% ATOM      1  N   SER C  14     -45.638 122.405  15.425  1.00 46.80           N                    
% HETATM41476  NA  CLA A 824     256.600 255.991 267.607  1.00 33.47           N  
% HETATM28407  O1D CHL 4 315      36.507 105.935 -33.192  1.00 96.36           O  

            atom_at_label = a{ii}( 77 : end );
            Nnum = 0;
            if ~isempty( findstr( atom_at_label , 'C' ) )
                Nnum = 1;                
            end
            if ~isempty( findstr( atom_at_label , 'O' ) )
                Nnum = 2;                
            end
            if ~isempty( findstr( atom_at_label , 'N' ) )
                Nnum = 3;                
            end
            if ~isempty( findstr( atom_at_label , 'S' ) )
                Nnum = 4;                
            end
            
%             disp( a{ii} )
            reslab = a{ii}( 23 : 29 );
            resnumlab = str2num( char( strread(reslab,'%s') ) );
%             disp( resnumlab );
%             chainlab = a{ii}(poslab+4);
%             atomindexinchain = a{ii}(poslab+5:30);
            
            XYZpart1 = a{ii}(31:38);
            XYZpart2 = a{ii}(39:47);
            XYZpart3 = a{ii}(48:56);        
            Coord1 = strread(XYZpart1,'%s');
            Coord2 = strread(XYZpart2,'%s');        
            Coord3 = strread(XYZpart3,'%s');        
            if ~isempty(Coord1) & Nnum ~= 0
                XYZ(jj,1) = str2num(char(Coord1));
                XYZ(jj,2) = str2num(char(Coord2));
                XYZ(jj,3) = str2num(char(Coord3));    
                XYZ(jj,4) = resnumlab; %double( chainlab);        
                XYZ(jj,5) = Nnum;%str2num( atomindexinchain );        
                
                jj = jj + 1;
            end
        end
        
        if ~isempty( findstr( a{ii} , 'CONNECT') )
            break
        end        
   
end


disp('...done')

%------------------------------------------------------------
function I = faread( FileName )
%
%
%
a = textread( FileName,'%s','delimiter', '\n');

jj = 1;
kk = 0;
ii = 1;
while 1

        if jj > size(a,1)
            disp('...EOF')
            break
        end     

    if length( a{jj}) ~= 0        
        if a{jj}(1) == '>' 
            kk = kk + 1;
            ii = 1;
        end
        
        if kk > 0
            I{kk,ii} = a{jj};
            ii = ii + 1;
        end
    end
    jj = jj + 1;
    
end
disp(FileName)
disp([num2str(size(a,1)), ' lines, ', num2str(kk), ' sequences read'])       

%-----------------------------------------------------
function v = myrainbow( x, min, max, a )
%
%
%

p = ( x - min ) / ( max - min ) * 100;
if x < min
    v = [0,0,0];
end
if x > max
    v = [0.9,0.9,0.9];
end
if x >= min & x <= max
    v(1) = exp( -0.5 *  ( (p - 25) / 25)^2 );
    v(2) = exp( -0.5 *  ( (p - 50) / 25)^2 );
    v(3) = exp( -0.5 *  ( (p - 75) / 25)^2 );
end
v = v * a;

