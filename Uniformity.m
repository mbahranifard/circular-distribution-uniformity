clear all
close all
eyename = ['test'];
BFNAME = ['BF.png'];
GRNAME = ['GR.png'];
% cd('C:\Users\mb   ahranifard3\Downloads\fiji-win64\Fiji.app\scripts')
% ImageJ;
cd ('C:\Users\mbahranifard3\Desktop\Uniformity\')
addpath ('C:\Users\mbahranifard3\Downloads\sc-toolbox-master')


% figure; imshow(transpose(modMB),[]);
imgbinary=imread(GRNAME);
% figure; imshow(imgbinary);
% axis xy
img=imread(BFNAME);
figure(1), imshow((img),[])
axis xy

bwmask=zeros(size(img,1),size(img,2));
bwmask=logical(bwmask);

for i=1
%     i=1:4
figure(1)
sprintf('Big countour Q %d',i)
roi = drawpolygon('Color','g','MarkerSize',1);
bw = createMask(roi);

% mask_boundary = boundarymask(bw);
% temp_MB = (mask_boundary);
% figure; imshow(temp_MB)
sprintf('vertices Q %d',i)
roiverx = drawpolygon('Color','b','MarkerSize',1);

ovpol=roi.Position;
overx = roiverx.Position;
% % figure()
% % plot(ovpol(:,1),ovpol(:,2))
r1=pdist2(ovpol(end,:),ovpol(end-1,:));
r2=pdist2(ovpol(end,:),ovpol(1,:));
r1trap=pdist2(overx(end-1,:),overx(end-2,:));
r2trap=pdist2(overx(end,:),overx(1,:));

rlimb = max(r1,r2);
rtrap = max(r1trap,r2trap);
difarc = ovpol(2:end-1,:)-ovpol(1:end-2,:);
eucarc = sqrt(sum(difarc.^2,2));
arclength = sum(eucarc);
r1vec = ovpol(end-1,:)-ovpol(end,:);
r2vec = ovpol(1,:)-ovpol(end,:);
theta(i) = acos(dot(r1vec,r2vec)/(r1*r2));
tetaref(i) = acos(dot(r1vec,[1,0])/r1);
if r1vec(2)<0
    tetaref(i) = 2*pi-tetaref(i);
end
minr=min(r1,r2)-rtrap;
r(i,:) = minr:(rlimb-minr)/10:rlimb;
ovpolygon{i} = ovpol;
bwmask(bw)=1;
end

% plot(xtot,ytot)

%apply the mask on the grayscale image

imgbinary(~bwmask)=0;


% temp_MB = uint8(temp_MB);
% img = imread('67ODQASl1N3_ch01.png');
% IJM.show('temp_MB');
% IJM.getDatasetAs('modMB');

% figure; imshow((imgbinary));
% axis xy
imwrite(imgbinary, 'CUTBOUND.png');
imshit = imbinarize(im2gray(imgbinary));
% figure()
% imshow(imshit)
% modMBMAT=imbinarize(imread('cutproc.png'));
modMBMAT = imshit;
if numel(modMBMAT(modMBMAT==1))>numel(modMBMAT)/2
    modMBMAT = abs(modMBMAT-1);
end
%plot location of bins


[X,Y]  = meshgrid([1:size(modMBMAT,1)],[1:size(modMBMAT,2)]);
X=X(:);
Y=Y(:);
illummat = [];
% modax = figure(6);
% imshow(modMBMAT)
% axis ij
for k=1:4
% angleconversion = sum(theta)/(2*pi);
%% change back1!
angleconversion = sum(theta)/(pi/2);
%%
% angleconversion=1;
tetasave = tetatot;
tetatot = theta(k)*angleconversion;
teta = [0:pi/12*angleconversion:theta(k)+0.05]+tetaref(k);
% tetasave = teta;
% teta=teta*.9998;
rez1 = r(1,2);
rezo = r(1,end);
rnew = linspace(rez1,rezo,10)
[Xnew,Ynew] = meshgrid(rnew(k,:),1i*sin(teta)+cos(teta));
circgrid = Xnew.*Ynew;
xc = real(circgrid);
yc = imag(circgrid); 
% hold on
% figure
xtot = xc+ovpolygon{k}(end,1);
ytot = yc+ovpolygon{k}(end,2);
ytotlin = reshape(ytot,[],1);
xtotlin = reshape(xtot,[],1);
figure(20)
if k==1
% imshow(modMBMAT)
imshow(img)
axis xy
end
hold on 
% scatter(xtotlin,ytotlin, 'gO','linestyle','-'); 
plot(xtot,ytot, 'g-O','linestyle','-','linewidth',1)
plot(xtot',ytot','g-','linewidth',1)
plot(xtot(:,end-2:end),ytot(:,end-2:end), 'b-O','linestyle','-','linewidth',1.5)
plot(xtot(:,end-2:end)',ytot(:,end-2:end)','b-','linewidth',1.5)

    for i = 1:size(xtot,2)-1
        for j = 1:size(xtot,1)-1
            boundpolX = xtot(j:j+1,i);
            boundpolX = [boundpolX;xtot(j+1:-1:j,i+1)];
            boundpolY = ytot(j:j+1,i);
            boundpolY = [boundpolY;ytot(j+1:-1:j,i+1)];
            in = inpolygon(X(:),Y(:),boundpolX,boundpolY);
            totpix = numel(X(in));
            idxmod = sub2ind(size(modMBMAT),Y(in),X(in)); %% for some reason the fucking modMBMAT matrix is the transpose of what being plotted
            litpix = numel(modMBMAT(modMBMAT(idxmod)==1));
            illummat{k}(i,j)=litpix/totpix;
    %         figure(6)
    %         hold on
    %         scatter(X(in),Y(in),10,'bO')
            clear boundpolX boundpolY
        end    
    end
end

illumatcat = cat(2,illummat{:})
savename = [eyename,'.csv']; 
writematrix(illumatcat,savename)



%% analysis

% hold all
% [xi,yi,~] = impixel() 
% r2= xi(2);
% c2= yi(2);
% % contourCCW = bwtraceboundary(modMBMAT,[r2 c2],'N',8,[],'clockwise');
% contourCCW = bwtraceboundary((modMBMAT),[r2 c2],'S');
% lctour = size(contourCCW,1);
% contourCCW = bwtraceboundary(modMBMAT,[r2 c2],'S',8,round(lctour/2)-50,'counterclockwise');
% % v= contourCCW(:,2);
% % contourCCW(:,2) = contourCCW(:,1);
% % contourCCW(:,1) = v;
% %%% the first two colums are swapped for contourCCW
% % contourCCW = bwtraceboundary(modMBMAT,[c2 r2],'N',4);
% plot(contourCCW(:,1),contourCCW(:,2),'b','LineWidth',2)
% plot(contourCCW(10,1),contourCCW(10,2),'go','markersize',5)

%% transformation
%reverse transformation
%sector meshgrid
% teta = [0:pi/30:pi/3];
% r = 0:0.1:1;
% [X,Y] = meshgrid(r,1i*sin(teta)+cos(teta));
% circgrid = X.*Y;
% xc = real(circgrid);
% yc = imag(circgrid);
% figure
% plot(xc,yc)
% p1=polygon([circgrid(1);circgrid(:,end)]);
% figure()
% plot(p1)
% f1 = diskmap(p1);
% f1=center(f1,circgrid(6,6))
% figure()
% plot(evalinv(f1,circgrid))
% 
% 
% %
% 
% p = roi.Position;
% pimg= p(:,1)+1i*p(:,2);
% ppol= polygon(pimg);
% % p = polygon(mask_boundary);
% % figure; plot(roi);
% f = diskmap(ppol);
% fcomp = composite(inv(f1),f);
% figure()
% plot(eval(fcomp,circgrid))
% hold on
% plot(pimg)
% figure()
% plot(eval(fcomp,p))
% 
% 
% figure()
% plot(f)
% figure()
% size(bw)
% figure()
% plot(bw)
% [X,Y]  = meshgrid([1:size(bw,1)],[1:size(bw,2)]);
% plot(evalinv(f,X+i*Y),'k')

% plot(exp(i*linspace(0,2*pi,180)));
% hold on, axis equal
% [X,Y] = meshgrid((-4:4)/5,(-100:100)/100);
% plot(evalinv(f,X+i*Y),'k')
% [X,Y] = meshgrid((-100:100)/100,(-4:4)/5);
% plot(evalinv(f,X'+i*Y'),'k')
% plot(p(:,1),p(:,2))
%%
% diffs = diff(contourCCW,1,1);
% ddiffs = diffs.^2;
% diffdist = sqrt(ddiffs(:,1)+ddiffs(:,2));
% diffdist (end+1) = 0;
% distmat = zeros(size(imgbinary));
% [Ximbin,Yimbin] = find(imgbinary);
% imbinelem = [Ximbin,Yimbin];
% 
% [k,eudist] = dsearchn(contourCCW,imbinelem);
% 
% figure
% plot(contourCCW(:,1),contourCCW(:,2)*-1,'ko')
% hold on
% plot(imbinelem(:,1),imbinelem(:,2)*-1,'*g')
% hold on
% plot(contourCCW(k,1),contourCCW(k,2)*-1,'*r')
% 
% cumdiffdist = cumsum(diffdist);
% totaldist = cumdiffdist(k) + eudist;
% 
% for i=1:length(k)
%     distmat(Ximbin(i),Yimbin(i))=totaldist(i);
% end

%% normalizing to iris length
% % normtotaldist = totaldist / cumdiffdist(end)/numel(k);
% normdistmat = distmat / cumdiffdist(end);
% 
% % distpar = sum(normtotaldist);
% distpar = sum(normdistmat(:))/numel(k);
% 
% 
% %% getting the intensity values
% img=imread(GRNAME);
% % img = rgb2gray(img);
% img = transpose(img);
% % figure, imshow(transpose(img),[])
% maximg = max(img(:));
% imgnorm = img*255/maximg;
% imbintemp = uint8(imgbinary);
% imgmasked = img .* imbintemp;
% imgmasked = double(imgmasked);
% 
% totint = sum(imgmasked(:));
% imgcorr = imgmasked / totint;
% % imgcorr = (img/255);
% 
% intnormdist = normdistmat .* imgcorr;
% intdistpar = sum(intnormdist(:))
% 
% imgnon = imgmasked(imgmasked~=0);