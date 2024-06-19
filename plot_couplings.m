function plot_couplings
%
%
%
clear all;
% close all;

% display mode:
showV = 1; % couplings in 1/cm will be shown 
showindex = 1.5 ; % 1 - pigment index will be shown
               % 2 - pig. number will be shown
               % 1.5 - combines options 1 and 2
               % 3 - show coupled pigs
               % 0 - show only marker 
               % -1 - not even the marker

%% *** inputs from file(s)
% fid = fopen( 'C:\DataVB\Progs\dip_models\dip_path.txt' );
fid = fopen( 'C:\Data\Progs\dip_models\dip_path.txt' );
fp0 = fgetl(fid);
fclose(fid);     

parpath = fp0; % master directory
savepath = fp0;

[fn,fp] = uigetfile( [fp0, '*.txt' ], 'dipoles' )
if fn == 0
    disp('no file, end OK')
    return
end

data = load( [fp,fn] ) 

[fn,fp] = uigetfile( [fp0, 'index*.txt' ], 'index?' )
if fn ~= 0
    index = load( [fp,fn] );
else
    index = ones( size(data,1),1)+0;
end

% colorcount = length(unique(index));
colorcount = 5;
for ii = 1:colorcount
    
    my_color(ii,:) = myrainbow( ii, 0, colorcount , 1, 0 ) ;
    
end

% there are two types of mRe files
if size( data ,2 ) > 7
   labeldip = num2str( data(:, 2));       
   mRe = data( :, 4:end);
else
    labeldip = num2str((1:size( data, 1 ))');    
    mRe = data;
end
% mRe(6,4) = mRe(6,4) +0.2
mRe(:,1) = 5.8;
w0 = 800;
coupcutoff = 60;

% mRe2 = load( ['C:\Data\2023_08\LHC_comp_08\6a2w_to_1rwt_1mer_allCLA_vecdata_qy.txt'] );
% mRe2 = mRe2(:, 4:end);
% mv = cross( mRe(9,5:7), mRe2(9,5:7)  )
% 
% mRe(9,2:4) = mRe(9,2:4) + mv*0.2; 
% mRe([2,6], 1) = 2.2;

% pick = [ 1,6,8 ]';
% mRe = mRe( pick, : );
% index = index(pick);
% 
% labeldip = labeldip(pick,:);
DC = 1;
n = 1;

disp( [(1:size( mRe ,1) )', mRe ])

% spectra computation
SE = 1e7./w0 * ones( size(mRe,1),1 );
HB = 300 * ones( size(mRe,1),1 );
SS = CalcSpc( mRe  , SE, DC, HB, (400:0.1:750)', 1); % stick spectra
disp('transition dipoles (effective) [D]:')
m = [1e7./SS(:,1),SS(:,2).^0.5 / mRe(1,1)];
disp(m)
%...

% return
S = Vij_calc_nofile(mRe, DC);
disp(S)

% return
% Rx = mRe(:,2) * ones( 1,size( mRe ,1) ) - ( mRe(:,2) * ones( 1,size( mRe ,1) ) )';
% Ry = mRe(:,3) * ones( 1,size( mRe ,1) ) - ( mRe(:,3) * ones( 1,size( mRe ,1) ) )';
% Rz = mRe(:,4) * ones( 1,size( mRe ,1) ) - ( mRe(:,4) * ones( 1,size( mRe ,1) ) )';
% Rnorm = ( Rx.^2 +  Ry.^2 +  Rz.^2 ).^0.5 

f2 = 1;
fi = 1;

% figure(1)
% title( fn, 'Interpreter', 'none' )
% for ii = 1:size( S, 1)
%     
%     for jj = 1+ii:size( S, 2)
%         
%         if abs( S(ii,jj) ) > 100
%             disp(['[',num2str(ii),', ',num2str(jj), '] = [', labeldip(ii,:),', ',labeldip(jj,:), '], ', num2str( S(ii,jj) ), ' [1/cm]'] )
%         end
%         intscale = abs( S(ii,jj) ) / max(max( abs(S) ));
%         if S(ii,jj) ~= 0;        
%             mc = myrainbow( abs( S(ii,jj) ), 0, max(max( abs(S) )), 1+(1-intscale)*5, 1 ) ;
%         else
%             mc = [0.6, 0.6, 0.6 ] ;
%         end
%         rectangle('Position',[ ii-0.5, jj-0.5, 1 , 1], 'FaceColor', mc);
%         text( ii-0.2, jj, num2str( round( S(ii,jj)*10 )/10 ) ) 
%         fi = fi + 1;        
%         
%     end
%     
%     if rem( ii, 20 ) == 0
%         f2 = f2 + 1;
%         figure(f2)  
%         pause(eps)
%     end    
%  
% end
% set(gca,'ytick', 1:size(S, 1) )
% set(gca,'xtick', 1:size(S, 1) )
% set(gca,'XTickLabel', labeldip( 1:size(S, 1),: ) )
% set(gca,'YTickLabel', labeldip( 1:size(S, 1),: ) )
% print([savepath,'300dpi.png'], '-dpng', '-r300');

% nonzeros(triu(S))

% figure(2)
% hist( abs(nonzeros(triu(S))) );
% return

figure(2)
zz = 1;
for ii = 1:size( S, 1)
    
    for jj = ii+1:size( S, 2)

        intscale = abs( S(ii,jj) ) / max(max( abs(S) ));
               
        if S(ii,jj) ~= 0        
            mc = myrainbow( abs( S(ii,jj) ), 0 , max(max( abs(S) )), 1+(1-intscale)*4, 1 ) ;
        else
            mc = [0.6, 0.6, 0.6 ] ;
        end
               
        if abs( S(ii,jj) ) > coupcutoff % cutoff on coupling
                r = mRe( jj,2:4 ) - mRe( ii,2:4 );
                stored(zz,:) = [ii,jj,S(ii,jj),norm(r)];
                zz = zz + 1;
                vx = [mRe( ii,2 )+0.1*r(1), mRe( jj,2 )-0.1*r(1) ] ;
                vy = [mRe( ii,3 )+0.1*r(2), mRe( jj,3 )-0.1*r(2) ] ;
                vz = [mRe( ii,4 )+0.1*r(3), mRe( jj,4 )-0.1*r(3) ] ;
                plot3( vx, vy ,vz, '-', 'Color', mc , 'LineWidth', 4 )
                if showV
                    text( 0.5*(mRe( ii,2 )+mRe( jj,2 )),...
                        0.5*(mRe( ii,3 )+mRe( jj,3 )),...
                        0.5*(mRe( ii,4 )+mRe( jj,4 )), ['   ', num2str( round(S(ii,jj))),' '],...
                        'Color', mc*0.4)
                end
                hold on
                disp( [num2str(ii),'->',num2str(jj),',',num2str( S(ii,jj) )] )
        end        
    end
    
end
for ii = 1:size( mRe, 1)
    switch showindex
        case 1
            text( mRe( ii,2 ), mRe( ii,3 ), mRe( ii,4 ), ['  ', labeldip( ii, :) ] )
        case 2
            text( mRe( ii,2 ), mRe( ii,3 ), mRe( ii,4 ), ['  ', num2str( ii )] )
        case 1.5
            text( mRe( ii,2 ), mRe( ii,3 ), mRe( ii,4 ), [num2str( ii ),'-', labeldip( ii, :) ] )
    end
    if showindex > -1
        plot3( mRe( ii,2 ), mRe( ii,3 ), mRe( ii,4 ), 'sq',...
                                                  'MarkerFaceColor', my_color(index(ii),:),...
                                                  'MarkerEdgeColor', my_color(index(ii),:))
    end
end

% for ii =1:size( mRe ,1)
%     text( mRe( ii,2 )-0.5, mRe( ii,3 )-0.5, mRe( ii,4 ), ['  ', num2str( ii )] )
% end

set(gca,'DataAspectRatio', [1 1 1])
% drawdips( mRe, [0,0,0], 0.5, [] );
grid on


% figure(3)
% plotsticks( SS, [0,0,1] );
% hold off

% figure(1)
% plot3( mRe(9,2)+[0;mv(1)],  mRe(9,3)+[0;mv(2)],  mRe(9,4)+[0;mv(3)], 'k-'  )

disp( stored )

[hac,hbc] = hist( stored(:,3), (250:10:400) );
[har,hbr] = hist( stored(:,4), (0.8:0.025:1.2) );
figure(4)

disp( [hbr(:),har(:)] )
disp( [hbc(:),hac(:)] )

subplot(1,2,1)
stairs(hbr,har,'LineWidth', 2);
hold on
subplot(1,2,2)
stairs(hbc,hac,'LineWidth', 2);
hold on

outr = [hbr(:),har(:)];
outc = [hbc(:),hac(:)];

save('C:\Data\2024_05\gemma_2024\hist_dist.txt', 'outr', '-ascii')
save('C:\Data\2024_05\gemma_2024\hist_coup.txt', 'outc', '-ascii')
save('C:\Data\2024_05\gemma_2024\stored.txt', 'stored', '-ascii')

disp( [num2str( mean( stored(:,3) ) ), ' +/- ', num2str(std( stored(:,3) )) ])
disp( [num2str( mean( stored(:,4) ) ), ' +/- ', num2str(std( stored(:,4) )) ])
disp('end OK')




%==========================================================================
function S = CalcSpc_var(vecdata, siteenergies, Vij, DC, HB, Xax, sticks_only)
% a stand-alone computation of absorbance and CD spectra
% the input data format is a mRe-like dataset 
% this version houses the subroutines StickSpec and DressedSticks that were
% elliminated from 'private' folder
Xwn = 1e7 ./ Xax;
R = vecdata(:,2:4);
e = vecdata(:,5:7);
m = vecdata(:,1);
%Vij = Vij_calc_nofile(vecdata, DC);
Ham = make_Ham(Vij, siteenergies);
[E, U] = ESE(Ham) ;

[A, CD] = StickSpecL( siteenergies, U, vecdata);

DS = DressedSticksL(Xwn, E, A, CD, HB);

S = [Xax, DS(:,1), DS(:,2) ];
if sticks_only == 1
    S = [E,A ,CD ]; % outputs dipole strength (mju^2) and rotational strength ( A/E and CD/E )
end

%-------------------------------------------------------------------------
function DS = DressedSticksL(Xaxis, E, A, CD, FWHM) 
% ver. 2022-10-26
% dressed sticks: gaussian function

FWHMi = FWHM;
if max(size(FWHM)) == 1
    FWHMi = E * 0 + 1 * FWHM;
end
h = Xaxis * 0 + 1;
h2 = E * 0 + 1;
Xm = Xaxis * h2';
Am = h * A';
CDm = h * CD';
Em = h * E';
fwhm = h * FWHMi'; %FWHM to standard deviation conversion

cA = 9.19e-3;
cCD = 7.5e-5;

q = ( sqrt(log(2)/pi) * 2 );
f = exp( -0.5 * ( ( Xm - Em )./(fwhm/2.355) ).^2);

cm = fwhm / 2.355;
A  = Am ./ (9.19e-3 * cm * 2.355).* exp( -0.5 * ( ( Xm - Em )./cm ).^2 )  * q  .* Xm;
CD = CDm   ./ (7.5e-5  * cm * 2.355).* exp( -0.5 * ( ( Xm - Em )./cm ).^2 )  * q    .*Xm;
% 
% A = Am .* Em ./ ( cA * fwhm ) .* f * q;
% CD = CDm .* Em ./ ( cCD * fwhm ) .* f * q;

% A = Am .* Em ./ ( cA * fwhm ) .* f * q;
% CD = CDm ./ ( cCD * fwhm ) .* f * q;

% A = Am   ./ ( cA * fwhm ) .* f * q;
% CD = CDm   ./ ( cCD * fwhm ) .* f * q;
% 

DS = [sum(A, 2) , sum(CD, 2) ] ;

%------------------------------------------------------
function [A, CD] = StickSpecL( E, U, mRe )
% ver. 2022-10-26
% calculates stick spectra
% [absorption, circular dichroism] =
% = f( energies, eigenvectors, mRe )
% again, mRe = [ debye, [positionX,Y,Z], [unit_directionX,Y,Z ] ], dimension: n x 7
% new, compact version, up to 20-25% faster; eliminated all 'for' cycles
% but one (the one goes over exciton levels)

N = size(E,1);
mvec = [];
Rvec = [];
mju = mRe(:,1) * [1,1,1] .* mRe( :, 5:7);

MU1 = mju(:,1) * ones( 1, length(mju(:,1)) ) ;
MU2 = mju(:,2) * ones( 1, length(mju(:,2)) );
MU3 = mju(:,3) * ones( 1, length(mju(:,3)) ) ;

R1 = mRe(:,2) * ones( 1, size(mRe,1) ) - ( mRe(:,2)  * ones( 1, size(mRe,1) ) )';
R2 = mRe(:,3) * ones( 1, size(mRe,1) ) - ( mRe(:,3)  * ones( 1, size(mRe,1) ) )';
R3 = mRe(:,4) * ones( 1, size(mRe,1) ) - ( mRe(:,4)  * ones( 1, size(mRe,1) ) )';

C1 = MU2.*MU3' - MU3.*MU2';
C2 = MU3.*MU1' - MU1.*MU3';
C3 = MU1.*MU2' - MU2.*MU1';

for kk = 1:N

    mjusq = mju(:,1) * mju(:,1)' + mju(:,2) * mju(:,2)' + mju(:,3) * mju(:,3)';
    A(kk,1) = sum( sum( mjusq .* (U(:,kk) * U(:,kk)') ) );       
    CD(kk,1)  = 1.7e-5 *  sum( E .* (sum( (R1.*C1 + R2.*C2 + R3.*C3) .* (U(:,kk)*U(:,kk)') ))' ) ;   
%     CD(kk,1)  = 1.7e-5 * sum(   (sum( (R1.*C1 + R2.*C2 + R3.*C3) .* (U(:,kk)*U(:,kk)') ))' ) ;   
end

% disp( 'A new')
% disp( A ) 
% disp( 'CD new')
% disp( CD ) 

% for kk = 1:N
%     
%     a(kk,1) = 0;
%     for ii = 1:N
% 
%         for jj = 1:N
% 
%             a(kk,1) = a(kk,1) + mju(ii,1)*U(ii,kk)*mju(jj,1)*U(jj,kk) + mju(ii,2)*U(ii,kk)*mju(jj,2)*U(jj,kk) +mju(ii,3)*U(ii,kk)*mju(jj,3)*U(jj,kk);  
%  
%         end
%     end
% end
% disp( 'current')
% disp( A )
% disp( 'silly')
% disp( a )

