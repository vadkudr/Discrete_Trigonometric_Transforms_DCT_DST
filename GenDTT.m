%% GenDTT.m
% TBD:
% 1. check roots of modulo polynomials
% 2. check alpha(k) and beta(k) (in some cases both are beta(k) and beta(k))

%% Table of roots
dtt_roots=cell(4);

dtt_roots{1,1}.name='DCTI';
dtt_roots{1,1}.var ='DCT1';
dtt_roots{1,1}.rootlist='[roots(TschebyshevU(N-2)); 1; -1]';

dtt_roots{1,2}.name='DCTIII';
dtt_roots{1,2}.var ='DCT3';
dtt_roots{1,2}.rootlist='[roots(TschebyshevT(N))]';

dtt_roots{1,3}.name='DCTV';
dtt_roots{1,3}.var ='DCT5';
dtt_roots{1,3}.rootlist='[roots(TschebyshevW(N-1)); 1]';

dtt_roots{1,4}.name='DCTVII';
dtt_roots{1,4}.var ='DCT7';
dtt_roots{1,4}.rootlist='[roots(TschebyshevV(N-1)); -1]';

dtt_roots{2,1}.name='DSTIII';
dtt_roots{2,1}.var ='DST3';
dtt_roots{2,1}.rootlist='[roots(TschebyshevT(N))]';

dtt_roots{2,2}.name='DSTI';
dtt_roots{2,2}.var ='DST1';
dtt_roots{2,2}.rootlist='[roots(TschebyshevU(N))]';

dtt_roots{2,3}.name='DSTVII';
dtt_roots{2,3}.var ='DST7';
dtt_roots{2,3}.rootlist='[roots(TschebyshevV(N))]';

dtt_roots{2,4}.name='DSTV';
dtt_roots{2,4}.var ='DST5';
dtt_roots{2,4}.rootlist='[roots(TschebyshevW(N))]';

dtt_roots{3,1}.name='DCTVI';
dtt_roots{3,1}.var ='DCT6';
dtt_roots{3,1}.rootlist='[roots(TschebyshevW(N-1)); 1]';

dtt_roots{3,2}.name='DCTVIII';
dtt_roots{3,2}.var ='DCT8';
dtt_roots{3,2}.rootlist='[roots(TschebyshevV(N))]';

dtt_roots{3,3}.name='DCTII';
dtt_roots{3,3}.var ='DCT2';
dtt_roots{3,3}.rootlist='[roots(TschebyshevU(N-1)); 1]';

dtt_roots{3,4}.name='DCTIV';
dtt_roots{3,4}.var ='DCT4';
dtt_roots{3,4}.rootlist='[roots(TschebyshevV(N))]';

dtt_roots{4,1}.name='DSTVIII';
dtt_roots{4,1}.var ='DST8';
dtt_roots{4,1}.rootlist='[roots(TschebyshevV(N-1)); -1]';

dtt_roots{4,2}.name='DSTVI';
dtt_roots{4,2}.var ='DST6';
dtt_roots{4,2}.rootlist='[roots(TschebyshevW(N))]';

dtt_roots{4,3}.name='DSTIV';
dtt_roots{4,3}.var ='DST4';
dtt_roots{4,3}.rootlist='[roots(TschebyshevT(N))]';

dtt_roots{4,4}.name='DSTII';
dtt_roots{4,4}.var ='DST2';
dtt_roots{4,4}.rootlist='[roots(TschebyshevU(N)); -1]';

%% polynomial type table
polytype={'TschebyshevT','TschebyshevU','TschebyshevV','TschebyshevW'};

%% scaling table
scaling={'cos(0*', 'sin(', 'cos(1/2*', 'sin(1/2*'};

%% DTT classic definition table
% cosine table
DCT_def=cell(1,8);

DCT_def{1}.teta='pi/(N-1)*k';
DCT_def{1}.freqidx='l';

DCT_def{2}.teta='pi/N*k';
DCT_def{2}.freqidx='(l+1/2)';

DCT_def{3}.teta='pi/N*(k+1/2)';
DCT_def{3}.freqidx='l';

DCT_def{4}.teta='pi/N*(k+1/2)';
DCT_def{4}.freqidx='(l+1/2)';

DCT_def{5}.teta='pi/(N-1/2)*k';
DCT_def{5}.freqidx='l';

DCT_def{6}.teta='pi/(N-1/2)*k';
DCT_def{6}.freqidx='(l+1/2)';

DCT_def{7}.teta='pi/(N-1/2)*(k+1/2)';
DCT_def{7}.freqidx='l';

DCT_def{8}.teta='pi/(N+1/2)*(k+1/2)';
DCT_def{8}.freqidx='(l+1/2)';


% sine table
DST_def=cell(1,8);

DST_def{1}.teta='pi/(N+1)*(k+1)';
DST_def{1}.freqidx='(l+1)';

DST_def{2}.teta='pi/N*(k+1)';
DST_def{2}.freqidx='(l+1/2)';

DST_def{3}.teta='pi/N*(k+1/2)';
DST_def{3}.freqidx='(l+1)';

DST_def{4}.teta='pi/N*(k+1/2)';
DST_def{4}.freqidx='(l+1/2)';

DST_def{5}.teta='pi/(N+1/2)*(k+1)';
DST_def{5}.freqidx='(l+1)';

DST_def{6}.teta='pi/(N+1/2)*(k+1)';
DST_def{6}.freqidx='(l+1/2)';

DST_def{7}.teta='pi/(N+1/2)*(k+1/2)';
DST_def{7}.freqidx='(l+1)';

DST_def{8}.teta='pi/(N-1/2)*(k+1/2)';
DST_def{8}.freqidx='(l+1/2)';

%% Tschebyshev polynomials relation table
Tsch_poly{1}.type1=1;    Tsch_poly{1}.type2=2;    Tsch_poly{1}.rel=1/2*[1 0 -1];
Tsch_poly{2}.type1=3;    Tsch_poly{2}.type2=2;    Tsch_poly{2}.rel=    [1 -1];
Tsch_poly{3}.type1=4;    Tsch_poly{3}.type2=2;    Tsch_poly{3}.rel=    [1  1];
Tsch_poly{4}.type1=1;    Tsch_poly{4}.type2=3;    Tsch_poly{4}.rel=1/2*[1  1];
Tsch_poly{5}.type1=1;    Tsch_poly{5}.type2=4;    Tsch_poly{5}.rel=1/2*[1 -1];

%% type of translation
T1='DSTVII'; T2='DSTVIII';

%% generate DTT relation m-file

% % % choosing row and column % % %
% row defines left boundary condition and scaling
% column defines right boundary condition
for i=1:numel(dtt_roots),
    if strcmp(dtt_roots{i}.name,T1), break; end
end
nT1=mod((i-1),4);   mT1=floor((i-1)/4);
for i=1:numel(dtt_roots),
    if strcmp(dtt_roots{i}.name,T2), break; end
end
nT2=mod((i-1),4);   mT2=floor((i-1)/4);
nT1=nT1+1;mT1=mT1+1;
nT2=nT2+1;mT2=mT2+1;

% % % set variables
T1var=dtt_roots{nT1,mT1}.var;
T2var=dtt_roots{nT2,mT2}.var;
idx1=str2num(T1var(end));
idx2=str2num(T2var(end));
type1=T1var(1:end-1);
type2=T2var(1:end-1);

%% preparing polynomial relation matrix

% whether table have relation for required types of polynomials?
idx=0;
for i=1:length(Tsch_poly),
    if  (Tsch_poly{i}.type1==nT1 && Tsch_poly{i}.type2==nT2) || ...
        (Tsch_poly{i}.type1==nT2 && Tsch_poly{i}.type2==nT1)
        idx=i;
        break;
    end
end

if idx~=0,

%     fidtmp=fopen('tmp,txt','w');
    bstr=[];
    
    %% determining transform direction
    if Tsch_poly{idx}.type1~=nT1,
        % exchange transform data
        temp=T1;T1=T2;T2=temp;
        temp=nT1;nT1=nT2;nT2=temp;
        temp=mT1;mT1=mT2;mT2=temp;
        temp=T1var;T1var=T2var;T2var=temp;
        temp=idx1;idx1=idx2;idx2=temp;
        temp=type1;type1=type2;type2=temp;
    end

    % prepare relation matrix
    rel=Tsch_poly{idx}.rel;
    % basic relation matrix
    N=8;
    B=toeplitz([rel(1) zeros(1,N-1)]',[rel zeros(1,N-length(rel))]);
%         fprintf(['B=toeplitz([' num2str(rel(1)) ' zeros(1,N-1)]'',[' num2str(rel) ' zeros(1,N-length(rel))]);\n']);
        bstr=[bstr 'B=toeplitz([' num2str(rel(1)) ' zeros(1,N-1)]'',[' num2str(rel) ' zeros(1,N-length([' num2str(rel) ']))]);\n'];
        B1=B;
    % left boundary condition
    B=left_bdry_condition(B,nT2);
    [T2size,T1size]=size(B);
        [row,col]=find(B~=B1);
        for i=1:length(row),
%             fprintf(['B(' num2str(row(i)) ',' num2str(col(i)) ')=' num2str(B(row(i),col(i))) ';\n']);
            bstr=[bstr 'B(' num2str(row(i)) ',' num2str(col(i)) ')=' num2str(B(row(i),col(i))) ';\n'];
        end

%     [n,~]=size(B);
%         bstr=[bstr '[n,~]=size(B);\n'];

    % increasing dimension (if necessary)
    eval(['rts='  dtt_roots{nT1,mT1}.rootlist ';']);
    eval(['rts2=' dtt_roots{nT2,mT2}.rootlist ';']);
    
    dimflag=0;  % default - no dimension increasing
    if length(rts)==length(rts2),
        if ~all(rts==rts2),
            if any(rts==1) && any(rts==-1),
                % x^2-1 polynomial exist in ring characteristic
                % add one dimension at the start of the matrix
                % and one dimension at the end of the matrix
%                 [n,m]=size(B);
                B=[B [zeros(2); B(1:end-2,end-1:end)]];
                        bstr=[bstr 'B=[B [zeros(2); B(1:end-2,end-1:end)]];\n'];
                        B1=B;
                    [T2size,T1size]=size(B);
                B=rot90(B,2); B=left_bdry_condition(B,mT2); B=rot90(B,2);
                        [row,col]=find(B~=B1);
                        for i=1:length(row),
                            bstr=[bstr 'B(' num2str(row(i)) ',' num2str(col(i)) ')=' num2str(B(row(i),col(i))) ';\n'];
                        end
                B=[ones(1,size(B,2)); B; (-1).^(0:size(B,2)-1)];
                            bstr=[bstr 'B=[ones(1,size(B,2)); B; (-1).^(0:size(B,2)-1)];\n'];
                dimflag=-11; % dimension is extended at both sides
            elseif any(rts==1),
                    % x-1 polynomial exist in ring characteristic
                    % add one dimension at the start of the matrix
%                     [n,m]=size(B);
                    B=[B [0; B(1:end-1,end)]];
                            bstr=[bstr 'B=[B [0; B(1:end-1,end)]];\n'];
                            B1=B;
                        [T2size,T1size]=size(B);
                    B=rot90(B,2); B=left_bdry_condition(B,mT2); B=rot90(B,2);
                            [row,col]=find(B~=B1);
                            for i=1:length(row),
                                bstr=[bstr 'B(' num2str(row(i)) ',' num2str(col(i)) ')=' num2str(B(row(i),col(i))) ';\n'];
                            end
                    B=[ones(1,size(B,2)); B];
                            bstr=[bstr 'B=[ones(1,size(B,2)); B];\n'];
                    dimflag=1; % dimension is extended at the start
            elseif any(rts==-1),
                    % x+1 polynomial exist in ring characteristic
                    % add one dimension at the end of the matrix
%                     [n,m]=size(B);
                    B=[B [0; B(1:end-1,end)]];
                            bstr=[bstr 'B=[B [0; B(1:end-1,end)]];\n'];
                            B1=B;
                        [T2size,T1size]=size(B);
                    B=rot90(B,2); B=left_bdry_condition(B,mT2); B=rot90(B,2);
                            [row,col]=find(B~=B1);
                            for i=1:length(row),
                                bstr=[bstr 'B(' num2str(row(i)) ',' num2str(col(i)) ')=' num2str(B(row(i),col(i))) ';\n'];
                            end
                    B=[B; (-1).^(0:size(B,2)-1)];
                            bstr=[bstr 'B=[B; (-1).^(0:size(B,2)-1)];\n'];
                    dimflag=-1; % dimension is extended at the end
            end
        end
    end
    
%     fclose(fidtmp);
else
%     B=eye(n);
    disp('No simple relation between basises. Stopping...');
    return;
end

%% open m-file for writing
fname = [T1 '_' T2 '.m'];
fid=fopen(fname,'w');

%% Writing header
fprintf(fid,['%%%% Relations between ' T1 ' and ' T2 '\n']);
fprintf(fid,['\n']);
fprintf(fid,['%%%% Definitions\n']);
fprintf(fid,['%% Transform matrix is defined for operating on column-vectors\n']);
fprintf(fid,['%% y=T*x,\n']);
fprintf(fid,['%% where y, x are column-vectors,\n']);
fprintf(fid,['%%       T    is transform matrix\n']);
fprintf(fid,['\n']);

%% analyzing transform sizes
% if nt1==nT1, [T1size,T2size]=size(B);
% else         [T2size,T1size]=size(B);
% end
% if nt1~=nT1, tmp=T1size;T1size=T2size;T2size=tmp; end
    
%% writing T1
fprintf(fid,['%%%% ' T1 ' matrix definition\n']);
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,['N1=' num2str(T1size) '; N=N1;\n']);
fprintf(fid,['k=0:N1-1;  l=0:N1-1;\n']);
if strcmp(type1,'DCT'), st1='cos'; teta1=DCT_def{idx1}.teta; freqidx1=DCT_def{idx1}.freqidx;
else                    st1='sin'; teta1=DST_def{idx1}.teta; freqidx1=DST_def{idx1}.freqidx;
end
fprintf(fid,[T1var '=' st1 '(' teta1 '''*' freqidx1 ')       %% display DCTV matrix\n']);
fprintf(fid,'\n');

ptype1=polytype{nT1};
fprintf(fid,['%%%% ' T1 ' in terms of Tschebyshev polynomials\n']);
fprintf(fid,['%% The ' T1 ' matrix can be expressed in terms of Tschebyshev polynomials [1]\n']);
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,['%% where\n']);
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,['%% are roots of polynomial\n']);
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
T1vart=[T1var 't']; rts=dtt_roots{nT1,mT1}.rootlist;
fprintf(fid,['alpha=sort(' rts ',''descend'');\n']);
fprintf(fid,[T1vart '=zeros(N1);\n']);
fprintf(fid,['for l=0:N1-1,\n']);
fprintf(fid,['    ' T1vart '(:,l+1)=polyval(' ptype1 '(l),alpha)'';\n']);
fprintf(fid,['end\n']);
D1var=['Da' num2str(idx1)]; sc1=scaling{nT1};
fprintf(fid,[D1var '=diag(' sc1 teta1 '));\n']);
fprintf(fid,[T1vart '=' D1var '*' T1vart '                    %% display DCTV matrix\n']);
fprintf(fid,'\n');

fprintf(fid,['%% compare ' T1var ' and ' T1vart ' matrices (show that both definitions above are equivalent)\n']);
fprintf(fid,['max(max(abs(' T1var '-' T1vart ')))\n']);
fprintf(fid,'\n');

%% writing T2
fprintf(fid,['%%%% ' T2 ' matrix definition\n']);
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,['N2=' num2str(T2size) '; N=N2;\n']);
fprintf(fid,['k=0:N2-1;  l=0:N2-1;\n']);
if strcmp(type2,'DCT'), st2='cos'; teta2=DCT_def{idx2}.teta; freqidx2=DCT_def{idx2}.freqidx;
else                    st2='sin'; teta2=DST_def{idx2}.teta; freqidx2=DST_def{idx2}.freqidx;
end
fprintf(fid,[T2var '=' st2 '(' teta2 '''*' freqidx2 ')       %% display DCTV matrix\n']);
fprintf(fid,'\n');

ptype2=polytype{nT2};
fprintf(fid,['%%%% ' T2 ' in terms of Tschebyshev polynomials\n']);
fprintf(fid,['%% The ' T2 ' matrix can be expressed in terms of Tschebyshev polynomials [1]\n']);
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,['%% where\n']);
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,['%% are roots of polynomial\n']);
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
T2vart=[T2var 't']; rts=dtt_roots{nT2,mT2}.rootlist;
fprintf(fid,['beta=sort(' rts ',''descend'');\n']);
fprintf(fid,[T2vart '=zeros(N);\n']);
fprintf(fid,['for l=0:N2-1,\n']);
fprintf(fid,['    ' T2vart '(:,l+1)=polyval(' ptype2 '(l),beta)'';\n']);
fprintf(fid,['end\n']);
D2var=['Db' num2str(idx2)]; sc2=scaling{nT2};
fprintf(fid,[D2var '=diag(' sc2 teta2 '));\n']);
fprintf(fid,[T2vart '=' D2var '*' T2vart '                    %% display DCTV matrix\n']);
fprintf(fid,'\n');

fprintf(fid,['%% compare ' T2var ' and ' T2vart ' matrices (show that both definitions above are equivalent)\n']);
fprintf(fid,['max(max(abs(' T2var '-' T2vart ')))\n']);

fprintf(fid,['\n']);

%% writing transform relation

fprintf(fid,['%%%% Finding relations\n']);
fprintf(fid,['%% Because there exist relation\n']);
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,['%% and\n']);
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,['%% we can express ' T1 ' through ' T2 '\n']);
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,['%% where\n']);
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');
fprintf(fid,'%%\n');

% fprintf(fid,['B=[ ...\n']);
% for i=1:size(B,2),
%     fprintf(fid,['%f,'],B(i,1:end-1));
%     fprintf(fid,['%f; ...\n'],B(i,end));
% end
% fprintf(fid,['];\n']);

% B matrix generation writing
fprintf(fid,bstr);

fprintf(fid,['\n']);

%% writing transform relation expressions

% forward relation
b='B';
% if nt1==nT1,
    d1=['Da' num2str(idx1)]; d2=['Db'  num2str(idx2)];
    t1var=T1var; t2var=T2var;
    t1=T1; t2=T2;
% else
%     d2=['Da' num2str(idx1)]; d1=['Db'  num2str(idx2)];
%     t2=T1; t1=T2;
%     t2var=T1var; t1var=T2var;
% end
fprintf(fid,['%%%% Check expression of ' t1 ' through ' t2 '\n']);
fprintf(fid,['%%\n']);
fprintf(fid,['%%\n']);
fprintf(fid,['%%\n']);
switch dimflag,
    case 1,
        im=['blkdiag(1,inv(' d2 ')*' t2var ')'];
    case -1,
        im=['blkdiag(inv(' d2 ')*' t2var ',1)'];
    case -11,
        im=['blkdiag(1,inv(' d2 ')*' t2var ',1)'];
    otherwise
        im=['inv(' d2 ')*' t2var];
end
t1expr=[d1 '*' im '*' b];
fprintf(fid,[t1var 'a=' t1expr '\n']);


fprintf(fid,['%% compare ' t1var ' and ' t1var 'a matrices (show correctness of representation of ' t1 ' through ' t2 ')\n']);
fprintf(fid,['max(max(abs(' t1var '-' t1var 'a)))\n']);

fprintf(fid,['\n']);

% inverse relation
b='inv(B)';
% if nt1~=nT1,
%     d1=['Da' num2str(idx1)]; d2=['Db'  num2str(idx2)];
%     t1var=T1var; t2var=T2var;
%     t1=T1; t2=T2;
% else
    d2=['Da' num2str(idx1)]; d1=['Db'  num2str(idx2)];
    t2=T1; t1=T2;
    t2var=T1var; t1var=T2var;
% end
fprintf(fid,['%%%% Check expression of ' t1 ' through ' t2 '\n']);
fprintf(fid,['%%\n']);
fprintf(fid,['%%\n']);
fprintf(fid,['%%\n']);
switch dimflag,
    case 1,
        im=['blkdiag(1,' d1 ')'];
        idx='2:end,2:end';
    case -1,
        im=['blkdiag(' d1 ',1)'];
        idx='1:end-1,1:end-1';
    case -11,
        im=['blkdiag(1,' d1 '1,)'];
        idx='2:end-1,2:end-1';
    otherwise
        im=d1;
        idx=':,:';
end
t2expr=[im '*inv(' d2 ')*' t2var '*' b];
fprintf(fid,[t1var 'a=' t2expr ';\n']);
fprintf(fid,[t1var 'a=' t1var 'a(' idx ')\n']);

fprintf(fid,['%% compare ' t1var ' and ' t1var 'a matrices (show correctness of representation of ' t1 ' through ' t2 ')\n']);
fprintf(fid,['max(max(abs(' t1var '-' t1var 'a)))\n']);

fprintf(fid,['\n']);

%% writing transform relation computation
fprintf(fid,['%%%% Check computation of ' T1 ' transform\n']);

fprintf(fid,['x=randn(N1,1);\n']);
fprintf(fid,['disp(''x''''='');disp(x'');\n']);
fprintf(fid,['y=' T1var '*x;                      %% true result\n']);
fprintf(fid,['disp(''y''''='');disp(y'');\n']);
% switch dimflag,
%     case 1,
%         xx='[0;x]';
%     case -1,
%         xx='[x;0]';
%     case -11,
%         xx='[0;x;0]';
%     otherwise
        xx='x';
% end
fprintf(fid,['y1=' t1expr '*' xx ';            %% compute ' T1 ' using ' T2 ' transform\n']);
fprintf(fid,['disp(''y1''''='');disp(y1'');\n']);
fprintf(fid,['\n']);


fprintf(fid,['%%%% Check computation of ' T2 ' transform\n']);

fprintf(fid,['x=randn(N2,1);\n']);
fprintf(fid,['disp(''x''''='');disp(x'');\n']);
fprintf(fid,['y=' T2var '*x;                      %% true result\n']);
fprintf(fid,['disp(''y''''='');disp(y'');\n']);
switch dimflag,
    case 1,
        xx='[0;x]';
    case -1,
        xx='[x;0]';
    case -11,
        xx='[0;x;0]';
    otherwise
        xx='x';
end
fprintf(fid,['y1=' t2expr '*' xx ';            %% compute ' T2 ' using ' T1 ' transform\n']);
fprintf(fid,['disp(''y1''''='');disp(y1'');\n']);

%% write reference
fprintf(fid,['%% Reference\n']);
fprintf(fid,['% [1] Markus Pueschel, Jose M.F. Moura. The Algebraic Approach to the Discrete\n']);
fprintf(fid,['% Cosine and Sine Transforms and their Fast Algorithms SIAM Journal of\n']);
fprintf(fid,['% Computing 2003, Vol. 32, No. 5, pp. 1280-1316.\n']);


%% close file
fclose(fid);
