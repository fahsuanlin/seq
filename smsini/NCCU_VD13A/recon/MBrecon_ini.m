function [ EPIsos,EPIrec ] = MBrecon_ini( filename1,filename2,lambda )

file_raw=filename1;
fprintf('raw data file = [%s]\n',file_raw);


clear global ice_obj;
clear global ice_m_data;


file_in = fopen(file_raw,'r','l','US-ASCII');
fseek(file_in,0,'eof');
fileSize = ftell(file_in);

fseek(file_in,0,'bof');

meas_ID  = fread(file_in,1,'uint32');
n_measraw = fread(file_in,1,'uint32');

meas_ID=fread(file_in,1,'uint32');
file_ID=fread(file_in,1,'uint32');
measOffset = fread(file_in,1,'uint64');
measLength = fread(file_in,1,'uint64');
patientName = fread(file_in,1,'uint64');
protocolName = fread(file_in,1,'uint64');


fseek(file_in,measOffset(1),'bof');

hdrLength  = fread(file_in,1,'uint32');

fprintf('measurement has a header of %d (bytes)\n', hdrLength);

buffer=fread(file_in,hdrLength,'uchar');

fp=fopen(sprintf('vd11_meas.asc'),'w');

fprintf(fp,'%c',char(buffer));
fclose(fp);
mrprot=ice_read_prot_mbsirepi(sprintf('vd11_meas.asc'));

fclose(file_in);

EPI=readVBVD(filename1);
EPI=EPI.sqzData;
EPIacc=readVBVD(filename2);
EPIacc=EPIacc.sqzData;
for ii=2:2:size(EPI,3)
    EPI(:,:,ii,:,:)=EPI(end:-1:1,:,ii,:,:);
    EPIacc(:,:,ii,:)=EPIacc(end:-1:1,:,ii,:);
end

EPI1_ref=EPI(1:size(EPI,1)/2,:,:,:,1);
ref1=squeeze(EPI1_ref(:,:,:,round(size(EPI1_ref,4)/2-1)+1,1));
ref1=fftshift(fft(fftshift(ref1,1),[],1),1);
newref1=zeros(size(ref1,1)*29,size(ref1,2),size(ref1,3));
newref1(size(ref1,1)*14+1:size(ref1,1)*15,:,:)=ref1;
newref1=ifftshift(ifft(ifftshift(newref1,1),[],1),1);
for ii=1:size(newref1,2)
    for jj=1:size(newref1,3)
        [~,phase1(ii,jj)]=max(abs(newref1(:,ii,jj)));
        phase1_cor(ii,jj)=angle(newref1(phase1(ii,jj),ii,jj));
        phase1(ii,jj)=(phase1(ii,jj)-((size(ref1,1)/2-1)*29+1))/29;
    end
end

phase1=mean(phase1(:,:),1);
phase=mean(phase1(2:2:end))-mean(phase1(1:2:end));

EPI1_ref=squeeze(mean(EPI(size(EPI,1)/2+1:size(EPI,1),:,:,:,4:end),5));
EPI2_ref=squeeze(mean(EPI(1:size(EPI,1)/2,:,:,:,4:end),5));
EPI1_acc=squeeze(EPIacc(size(EPI,1)/2+1:size(EPI,1),:,:,:));
EPI2_acc=squeeze(EPIacc(1:size(EPI,1)/2,:,:,:));


dSlice=sqrt((mrprot.dCor(5)-mrprot.dCor(1)).^2+(mrprot.dTra(5)-mrprot.dTra(1)).^2);
center=[mrprot.dCor(1)-(mrprot.dCor(end)-mrprot.dCor(1))*((mrprot.dCor(1)*(mrprot.dCor(end)-mrprot.dCor(1))+mrprot.dTra(1)*(mrprot.dTra(end)-mrprot.dTra(1)))/((mrprot.dCor(end)-mrprot.dCor(1))^2+(mrprot.dTra(end)-mrprot.dTra(1))^2)),mrprot.dTra(1)-(mrprot.dTra(end)-mrprot.dTra(1))*((mrprot.dCor(1)*(mrprot.dCor(end)-mrprot.dCor(1))+mrprot.dTra(1)*(mrprot.dTra(end)-mrprot.dTra(1)))/((mrprot.dCor(end)-mrprot.dCor(1))^2+(mrprot.dTra(end)-mrprot.dTra(1))^2))];
d_s1=sign([mrprot.dCor(5)-mrprot.dCor(1),mrprot.dTra(5)-mrprot.dTra(1)]*[mrprot.dCor(1+2*(round(mrprot.lSlices/2/2-1)))-center(1);mrprot.dTra(1+2*(round(mrprot.lSlices/2/2-1)))-center(2)])*sqrt((mrprot.dCor(1+2*(round(mrprot.lSlices/2/2-1)))-center(1)).^2+(mrprot.dTra(1+2*(round(mrprot.lSlices/2/2-1)))-center(2)).^2);
d_s2=sign([mrprot.dCor(5)-mrprot.dCor(1),mrprot.dTra(5)-mrprot.dTra(1)]*[mrprot.dCor(2+2*(round(mrprot.lSlices/2/2-1)))-center(1);mrprot.dTra(2+2*(round(mrprot.lSlices/2/2-1)))-center(2)])*sqrt((mrprot.dCor(2+2*(round(mrprot.lSlices/2/2-1)))-center(1)).^2+(mrprot.dTra(2+2*(round(mrprot.lSlices/2/2-1)))-center(2)).^2);
d_center=sqrt(center(1).^2+center(2).^2);

StartStep=floor((1/3/2)/(1/10));

for jj=1:size(EPI,4)
    for ii=1:size(EPI,3)
        if mod(ii,3)==2
            if floor((jj-1)/10*3) <2
                EPI1_ref(:,:,ii,jj)=EPI1_ref(:,:,ii,jj)*exp(1i*2*pi/3*d_s1/dSlice);
                EPI2_ref(:,:,ii,jj)=EPI2_ref(:,:,ii,jj)*exp(1i*2*pi/3*d_s2/dSlice);
            else
                EPI1_ref(:,:,ii,jj)=EPI1_ref(:,:,ii,jj)*exp(-1i*4*pi/3*d_s1/dSlice);
                EPI2_ref(:,:,ii,jj)=EPI2_ref(:,:,ii,jj)*exp(-1i*4*pi/3*d_s2/dSlice);
            end
        elseif mod(ii,3)==0
            if floor((jj-1)/10*3) <1
                EPI1_ref(:,:,ii,jj)=EPI1_ref(:,:,ii,jj)*exp(1i*4*pi/3*d_s1/dSlice);
                EPI2_ref(:,:,ii,jj)=EPI2_ref(:,:,ii,jj)*exp(1i*4*pi/3*d_s2/dSlice);
            else
                EPI1_ref(:,:,ii,jj)=EPI1_ref(:,:,ii,jj)*exp(-1i*2*pi/3*d_s1/dSlice);
                EPI2_ref(:,:,ii,jj)=EPI2_ref(:,:,ii,jj)*exp(-1i*2*pi/3*d_s2/dSlice);
            end
        end
    end
end


for ii=1:size(EPI,3)
    if mod(ii,3)==2
        if floor(StartStep/10*3) <2
            EPI1_acc(:,:,ii,:)=EPI1_acc(:,:,ii,:)*exp(1i*2*pi/3*d_s1/dSlice);
            EPI2_acc(:,:,ii,:)=EPI2_acc(:,:,ii,:)*exp(1i*2*pi/3*d_s2/dSlice);
        else
            EPI1_acc(:,:,ii,:)=EPI1_acc(:,:,ii,:)*exp(-1i*4*pi/3*d_s1/dSlice);
            EPI2_acc(:,:,ii,:)=EPI2_acc(:,:,ii,:)*exp(-1i*4*pi/3*d_s2/dSlice);
        end
    elseif mod(ii,3)==0
        if floor(StartStep/10*3) <1
            EPI1_acc(:,:,ii,:)=EPI1_acc(:,:,ii,:)*exp(1i*4*pi/3*d_s1/dSlice);
            EPI2_acc(:,:,ii,:)=EPI2_acc(:,:,ii,:)*exp(1i*4*pi/3*d_s2/dSlice);
        else
            EPI1_acc(:,:,ii,:)=EPI1_acc(:,:,ii,:)*exp(-1i*2*pi/3*d_s1/dSlice);
            EPI2_acc(:,:,ii,:)=EPI2_acc(:,:,ii,:)*exp(-1i*2*pi/3*d_s2/dSlice);
        end
    end
end

for ii=1:size(EPI,4)
    EPI1_ref(:,:,:,ii)=EPI1_ref(:,:,:,ii)*exp(1i*(-1/2+(ii-1)/size(EPI,4))*2*pi*d_s1/dSlice);
    EPI2_ref(:,:,:,ii)=EPI2_ref(:,:,:,ii)*exp(1i*(-1/2+(ii-1)/size(EPI,4))*2*pi*d_s2/dSlice);
end


EPI1_ref=fftshift(fft(fftshift(EPI1_ref,1),[],1),1);
EPI2_ref=fftshift(fft(fftshift(EPI2_ref,1),[],1),1);

EPI1_acc=fftshift(fft(fftshift(EPI1_acc,1),[],1),1);
EPI2_acc=fftshift(fft(fftshift(EPI2_acc,1),[],1),1);



for ii=2:2:size(EPI1_ref,3)
    for jj=1:size(EPI1_ref,1)
        EPI1_ref(jj,:,ii,:)=EPI1_ref(jj,:,ii,:)*exp(1i*phase*jj*2*pi/size(EPI1_ref,1)-1i*phase*(size(EPI1_ref,1)/2+1)*2*pi/size(EPI1_ref,1));
        EPI2_ref(jj,:,ii,:)=EPI2_ref(jj,:,ii,:)*exp(1i*phase*jj*2*pi/size(EPI1_ref,1)-1i*phase*(size(EPI1_ref,1)/2+1)*2*pi/size(EPI1_ref,1));
        EPI1_acc(jj,:,ii,:)=EPI1_acc(jj,:,ii,:)*exp(1i*phase*jj*2*pi/size(EPI1_ref,1)-1i*phase*(size(EPI1_ref,1)/2+1)*2*pi/size(EPI1_ref,1));
        EPI2_acc(jj,:,ii,:)=EPI2_acc(jj,:,ii,:)*exp(1i*phase*jj*2*pi/size(EPI1_ref,1)-1i*phase*(size(EPI1_ref,1)/2+1)*2*pi/size(EPI1_ref,1));
        
    end
end

EPI1_ref=fftshift(fft(fftshift(EPI1_ref,3),[],3),3);
EPI1_ref=fftshift(fft(fftshift(EPI1_ref,4),[],4),4);
EPI2_ref=fftshift(fft(fftshift(EPI2_ref,3),[],3),3);
EPI2_ref=fftshift(fft(fftshift(EPI2_ref,4),[],4),4);


EPI1_acc=fftshift(fft(fftshift(EPI1_acc,3),[],3),3);
EPI2_acc=fftshift(fft(fftshift(EPI2_acc,3),[],3),3);


ReadRes=mrprot.lReadFoV/mrprot.lBaseResolution;
PhaseRes=mrprot.lPhaseFoV/mrprot.lPhaseEncodingLines;
PixelShift=round(d_center/ReadRes);

EPI1_ref=permute(EPI1_ref(round(mrprot.lBaseResolution/2)+1+PixelShift:round(mrprot.lBaseResolution/2)+mrprot.lBaseResolution+PixelShift,:,:,:),[1,3,2,4]);
EPI2_ref=permute(EPI2_ref(round(mrprot.lBaseResolution/2)+1+PixelShift:round(mrprot.lBaseResolution/2)+mrprot.lBaseResolution+PixelShift,:,:,:),[1,3,2,4]);
EPI1_acc=permute(EPI1_acc(round(mrprot.lBaseResolution/2)+1+PixelShift:round(mrprot.lBaseResolution/2)+mrprot.lBaseResolution+PixelShift,:,:,:),[1,3,2,4]);
EPI2_acc=permute(EPI2_acc(round(mrprot.lBaseResolution/2)+1+PixelShift:round(mrprot.lBaseResolution/2)+mrprot.lBaseResolution+PixelShift,:,:,:),[1,3,2,4]);

EPI1_rec=zeros(mrprot.lBaseResolution,mrprot.lPhaseEncodingLines,size(EPI1_ref,4),size(EPI1_acc,4));
EPI2_rec=zeros(mrprot.lBaseResolution,mrprot.lPhaseEncodingLines,size(EPI2_ref,4),size(EPI2_acc,4));

for jj=1:mrprot.lBaseResolution
    for kk=1:mrprot.lPhaseEncodingLines
        temp1=squeeze(EPI1_ref(jj,kk,:,:));
        EPI1_rec(jj,kk,:,:)=(temp1'*temp1+lambda*eye(size(EPI1_rec,3)))\(temp1'*squeeze(EPI1_acc(jj,kk,:,:)));
        temp2=squeeze(EPI2_ref(jj,kk,:,:));
        EPI2_rec(jj,kk,:,:)=(temp2'*temp2+lambda*eye(size(EPI2_rec,3)))\(temp2'*squeeze(EPI2_acc(jj,kk,:,:)));
    end
end

EPI1_sos=sqrt(sum(abs(EPI1_ref).^2,3));
EPI2_sos=sqrt(sum(abs(EPI2_ref).^2,3));


EPI1_rec=shiftback_rec(EPI1_rec);
EPI2_rec=shiftback_rec(EPI2_rec);


EPI1_sos=shiftback_ref(EPI1_sos);
EPI2_sos=shiftback_ref(EPI2_sos);


EPIsos=zeros(mrprot.lBaseResolution,mrprot.lPhaseEncodingLines,1,size(EPI1_sos,4)*2);
EPIsos(:,:,:,1:2:end)=EPI2_sos;
EPIsos(:,:,:,2:2:end)=EPI1_sos;
EPIrec=zeros(mrprot.lBaseResolution,mrprot.lPhaseEncodingLines,size(EPI1_rec,3)*2,size(EPI1_rec,4));
EPIrec(:,:,1:2:end,:)=EPI2_rec;
EPIrec(:,:,2:2:end,:)=EPI1_rec;



    function Image=shiftback_ref(Image)
        Imagetemp=repmat(Image,[1,2,1,1]);
        for ii=1:size(Image,4)
            Image(:,:,:,ii)=Imagetemp(:,mod(ii-round(size(Image,4)/2+1/2),3)*size(Image,2)/3+1:mod(ii-round(size(Image,4)/2+1/2),3)*size(Image,2)/3+size(Image,2),:,ii);
        end
    end
    function Image=shiftback_rec(Image)
        Imagetemp=repmat(Image,[1,2,1,1]);
        for ii=1:size(Image,3)
            Image(:,:,ii,:)=Imagetemp(:,mod(ii-round(size(Image,3)/2+1/2),3)*size(Image,2)/3+1:mod(ii-round(size(Image,3)/2+1/2),3)*size(Image,2)/3+size(Image,2),ii,:);
        end
    end
end

